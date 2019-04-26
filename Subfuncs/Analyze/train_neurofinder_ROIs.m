function trainedModel = train_neurofinder_ROIs( opt )
%TRAIN_NEUROFINDER_ROIS Train a model that subselects and predicts shape
%for neurofinder ROIs

% Load the training data
orig_regions = py.neurofinder.load([fileparts(fileparts(fileparts(get_path(opt)))) filesep 'regions' filesep 'regions.json']);
orig_ROIs = cellfun(@matpy.nparray2mat, cell(orig_regions.coordinates), 'UniformOutput', false);

% Get the centers for the training ROIs
orig_ROI_centers = cellfun(@(X)round(mean((X-1)*opt.spatial_scale+1))', orig_ROIs, 'UniformOutput', false);
orig_H = [cell2mat(orig_ROI_centers)' ones(length(orig_ROI_centers),1)];

% Load our predicted ROIs
use_iter = opt.niter;

% Load model for use_iter
load(get_path(opt,'output_iter',use_iter),'model');

[H, W, X, y_orig, y, L] = model.get_fields( 'H', 'W', 'X', 'y_orig','y','L');
szWY = size(model.y);

% Load the data
load(get_path(opt),'inp');
data = inp.data;

% if ~exist('use_cells', 'var')
%   use_cells = 1:size(H,1);
% end
      
% [ROI_mask, ROIs] = getROIs(model.opt);



% Get the pairwise ROI distance matrix between true and current ROIs
[knn_inds, knn_dists] = knnsearch(orig_H(:,1:2), H(:,1:2), 'K', 1);

% Create a classifier from X to distance to nearest actual ROI (so
% basically whether we can learn from reconstructuion weigths whether its a
% true ROI or not)
% classifTargets = (exp(-knn_dists/(opt.m/5)))>0.5

%classifTargets = knn_dists <= 5*opt.spatial_scale;
classifTargets = exp(-knn_dists ./ (5*opt.spatial_scale));
classifTargetsBin = classifTargets>=(exp(-0.8));


data_moments = cell(opt.mom,1);
for mom1 = 1:opt.mom
   %data_moments{mom1,1} = mean((inp.data.proc_stack.Y).^mom1,3);
   data_moments{mom1,1} = mean((data.proc_stack.Y).^mom1,3);
end
data_cumulants = raw2cum(data_moments);

full_feature_im = cat(3, data_cumulants{:});

% Get the cumulant images as classifFeatures
classifFeatures = zeros(size(H,1), opt.m^2 * size(full_feature_im,3));

for i1 = 1:size(H,1)
  featureIm = nan(opt.m, opt.m, size(full_feature_im,3));
  
  % Get feature coordinates
  row = H(i1, 1); col = H(i1, 2);
  [inds, cut] = mat_boundary(szWY(1:2),row-floor(opt.m/2):row+floor(opt.m/2),col-floor(opt.m/2):col+floor(opt.m/2));
  featureIm(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2), :) = full_feature_im(inds{1},inds{2}, :);
  
  classifFeatures(i1, :) = featureIm(:);
end

%classifFeatures = model.X;

classifModel = TreeBagger(30, classifFeatures, classifTargets, 'Method', 'regression');
trainedModel{1} = classifModel;


%% Given the selected ROIs, learn how to map the raw data around suggested center to the binary ROI

% Collect the features and targets
data_moments = cell(opt.mom,1);
for mom1 = 1:opt.mom
   %data_moments{mom1,1} = mean((inp.data.proc_stack.Y).^mom1,3);
   data_moments{mom1,1} = mean((data.proc_stack.Y).^mom1,3);
end
data_cumulants = raw2cum(data_moments);

full_feature_im = cat(3, data_cumulants{:});

% Get the reconstruction around matched locations
regrTargetIndices = knn_inds(classifTargetsBin);
tmp = 1:size(H,1);
regrFeatureIndices = tmp(classifTargetsBin);
regrTargetSizeOrig = int16(((opt.m/opt.spatial_scale)) + 1 - mod((opt.m/opt.spatial_scale),2));
regrTargetSize = 7;
regrTargets = zeros(length(regrTargetIndices), regrTargetSize^2);
regrFeatures = zeros(length(regrTargetIndices), opt.m^2 * size(full_feature_im,3)); 
regrModels = cell(regrTargetSize,1);



for i1 = 1:length(regrTargetIndices)
  featureIm = nan(opt.m, opt.m, size(full_feature_im,3));
  targetIm = zeros(regrTargetSizeOrig, regrTargetSizeOrig);
  
  % Get feature coordinates
  row = H(regrFeatureIndices(i1), 1); col = H(regrFeatureIndices(i1), 2);
  [inds, cut] = mat_boundary(szWY(1:2),row-floor(opt.m/2):row+floor(opt.m/2),col-floor(opt.m/2):col+floor(opt.m/2));
  featureIm(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2), :) = full_feature_im(inds{1},inds{2}, :);
  
  regrFeatures(i1, :) = featureIm(:);
  
  % Get target coordinates
  row_origsize = round((row-1)./opt.spatial_scale)+1;
  col_origsize = round((col-1)./opt.spatial_scale)+1;
  targetCoords = orig_ROIs{regrTargetIndices(i1)} + (double(regrTargetSizeOrig)+1)/2 - repmat([row_origsize col_origsize], size(orig_ROIs{regrTargetIndices(i1)},1),1);
  j1 = 1;
  while j1 <= size(targetCoords, 1)
    if (targetCoords(j1,1) < 1) || (targetCoords(j1,2) < 1) || (targetCoords(j1,1) > regrTargetSizeOrig) || (targetCoords(j1,2) > regrTargetSizeOrig)
      targetCoords(j1,:) = [];
    else
      j1 = j1+1;
    end
  end
  targetIm(sub2ind([regrTargetSizeOrig, regrTargetSizeOrig], targetCoords(:,1), targetCoords(:,2))) = 1;

  targetImSmall = imresize(targetIm, [regrTargetSize, regrTargetSize], 'nearest');
  
  regrTargets(i1,:) = targetImSmall(:);
  
%   if opt.fig > 3
%     figure(117); 
%     subplot(1,2,1); imagesc(featureIm); axis image;
%     subplot(1,2,2); imagesc(targetIm); axis image;
%     pause;
%   end
end
  


parfor m1 = 1:(regrTargetSize^2)
  %regrModels{m1} = TreeBagger(opt.m^2, regrFeatures, regrTargets(:,m1));
  regrModels{m1} = TreeBagger(10, regrFeatures, regrTargets(:, m1), 'Method', 'regression');
end

trainedModel{2} = regrModels;

end

