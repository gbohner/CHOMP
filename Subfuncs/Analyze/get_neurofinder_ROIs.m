function [ROI_image, ROIs] = get_neurofinder_ROIs( opt, trainedModel, varargin )
%SUBSELECT_ROIS Given the neurofinder training data, subselect amongst the
%proposed location (based on training data)

if nargin <= 2
  classifThreshold = (exp(-1)*0.8);
else
  classifThreshold = varargin{1};
end

% Load our predicted ROIs
use_iter = opt.niter;

% Load model for use_iter
load(get_path(opt,'output_iter',use_iter),'model');

[H, W, X, y_orig, y, L] = model.get_fields( 'H', 'W', 'X', 'y_orig','y','L');
szWY = size(model.y);

% Load the data
load(get_path(opt),'inp');
data = inp.data;

szRaw = size(y_orig);
ROI_image = zeros([szRaw(1:2), 3]);



% Collect the features
data_moments = cell(opt.mom,1);
for mom1 = 1:opt.mom
   %data_moments{mom1,1} = mean((inp.data.proc_stack.Y).^mom1,3);
   data_moments{mom1,1} = mean((data.proc_stack.Y).^mom1,3);
end
data_cumulants = raw2cum(data_moments);

full_feature_im = cat(3, data_cumulants{:});

regrFeatures = zeros(size(H,1), opt.m^2 * size(full_feature_im,3)); 

for i1 = 1:size(H,1)
  featureIm = nan(opt.m, opt.m, size(full_feature_im,3));
  
  % Get feature coordinates
  row = H(i1, 1); col = H(i1, 2);
  [inds, cut] = mat_boundary(szWY(1:2),row-floor(opt.m/2):row+floor(opt.m/2),col-floor(opt.m/2):col+floor(opt.m/2));
  featureIm(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2), :) = full_feature_im(inds{1},inds{2}, :);
  
  regrFeatures(i1, :) = featureIm(:);
end

% Get the valid ROIs based on the classifier
% Subselect ROIs based on the trained classifier
classifPred = predict(trainedModel{1}, regrFeatures);
% classifPred = predict(trainedModel{1}, model.X);

%classifPred = str2num(cell2mat(classifPred));
use_cell = (classifPred>classifThreshold);
H = H(use_cell,:);

% Get the predicted ROI images
regrTargetSize = sqrt(length(trainedModel{2}));
predIms = zeros(size(H,1), regrTargetSize^2);
parfor m1 = 1:regrTargetSize^2
  tmp = predict(trainedModel{2}{m1}, regrFeatures(use_cell,:));
  %predIms(:,m1) = str2num(cell2mat(tmp));
  predIms(:,m1) = tmp;
end

% Convert the images back to original ROI size and store in the usual ROI
% struct
for i1 = 1:size(H,1)
  row = H(i1, 1); col = H(i1, 2); type = H(i1, 3);
  row = round((row-1)./opt.spatial_scale)+1;
  col = round((col-1)./opt.spatial_scale)+1;
  cutsize = round(opt.m./opt.spatial_scale);
  cutsize = cutsize + (1-mod(cutsize,2));
  
  reconst = reshape(predIms(i1,:),regrTargetSize,regrTargetSize);
  
  reconst = imresize(reconst,[cutsize,cutsize], 'bicubic');
  
  reconst = reconst>0.5;

%   
%   % Close the image
%   reconst = imopen(reconst,strel('rectangle',[3,3]));
%   reconst = imclose(reconst,strel('rectangle',[3,3]));
  
  if sum(reconst(:)) == 0
    reconst((numel(reconst)+1)/2) = 1;
  end
  
  
  
  [inds, cut] = mat_boundary(szRaw(1:2),row-floor(cutsize/2):row+floor(cutsize/2),col-floor(cutsize/2):col+floor(cutsize/2));
  ROI_image(inds{1},inds{2}, :) = ROI_image(inds{1},inds{2}, :) + ...
    mply(reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2)), 0.7*randn(1,3));
  
  
  %Store results in cells now, later add option for json output %TODO
  ROIs{i1} = struct('col', col, 'row', row, 'type', type, 'mask', reconst, 'reconst', zeros(size(reconst)));
end

% if opt.fig > 3
%   for i1 = 1:size(predIms,1)
%     figure(118); 
%     subplot(1,3,1); imagesc(reshape(regrFeatures(2*i1,:),opt.m, opt.m)); axis image;
%     subplot(1,3,2); imagesc(reshape(predIms(i1,:), regrTargetSize, regrTargetSize)); axis image;
%     subplot(1,3,3); imagesc(reshape(regrTargets(2*i1,:), regrTargetSize, regrTargetSize)); axis image;
%     pause;
%   end
% end


end

