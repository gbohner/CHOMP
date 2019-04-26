if ~exist('use_iter', 'var')
  use_iter = 1;
end

% Load model for use_iter
load(get_path(opt,'output_iter',use_iter),'model');

[H, W, X, y_orig, y, L] = model.get_fields( 'H', 'W', 'X', 'y_orig','y','L');
szWY = size(model.y);

if ~exist('use_cells', 'var')
  use_cells = 1:size(H,1);
end

update_visualize(model.y,model.H(use_cells,:), ...
  reshape(model.W,model.opt.m,model.opt.m,size(model.W,ndims(model.W))),...
  model.opt,1);

% Set temporary options
opt_tmp = model.opt;
opt_tmp.niter = use_iter;
      
[ROI_mask, ROIs] = getROIs(opt_tmp, use_cells);

% Show original mean image
figure(2); clf; imagesc(model.y); axis image; colorbar;

% Show reconstruction of mean
[rh,rl,Wfull] = reconstruct_cell( opt, W, X(use_cells,:));

full_reconst = zeros(szWY(1:2));
for c1 = 1:size(rl{1},3)
  row_hat = H(c1,1); col_hat = H(c1,2);
  [inds, cut] = mat_boundary(size(full_reconst),row_hat-floor(opt.m/2):row_hat+floor(opt.m/2),col_hat-floor(opt.m/2):col_hat+floor(opt.m/2));
  full_reconst(inds{1},inds{2}) = full_reconst(inds{1},inds{2}) + rl{1}(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2),c1);
end

figure(4); imagesc(full_reconst); axis image; colorbar;

% Show maximum moment image

if opt_tmp.mom >= 2
  load(get_path(opt_tmp)); % Load the input data
  data = inp.data;
  % Get the data cumulants
  data_moments = cell(opt_tmp.mom,1);
  for mom1 = 1:opt_tmp.mom
     %data_moments{mom1,1} = mean((inp.data.proc_stack.Y).^mom1,3);
     data_moments{mom1,1} = mean((data.proc_stack.Y).^mom1,3);
  end
  data_cumulants = raw2cum(data_moments);
  
%   data_moments_raw = cell(opt_tmp.mom,1);
%   raw_data = single(inp.data.raw_stack.Y(:,:,:));
%   for mom1 = 1:opt_tmp.mom
%      data_moments_raw{mom1,1} = mean(raw_data.^mom1,3);
%   end
%   data_cumulants_raw = raw2cum(data_moments_raw);
%   
%   % Whiten raw data and get cumulants
%   opt_tmp2 = opt_tmp;
%   opt_tmp2.m = round(opt_tmp2.m ./ opt.spatial_scale);
%   opt_tmp2.smooth_filter_mean = opt_tmp2.smooth_filter_mean ./ opt.spatial_scale;
%   a_tmp.proc_stack.Y = raw_data;
%   raw_data_white = whiten_proc(opt_tmp2, a_tmp);
%   raw_data_white = raw_data_white.proc_stack.Y;
%   
%   data_moments_raw_white = cell(opt_tmp.mom,1);
%   for mom1 = 1:opt_tmp.mom
%      data_moments_raw_white{mom1,1} = mean(raw_data_white.^mom1,3);
%   end
%   data_cumulants_raw_white = raw2cum(data_moments_raw_white);
  
  
  %(alternative, exactly the same - as tested)
  % W_point = zeros(opt.m);
  % W_point((opt.m+1)/2, (opt.m+1)/2) = 1;
  % [data_cumulantsWY, GW, WnormInv] = compute_filters(inp.data, repmat(W_point(:),[1, opt_tmp.KS]), opt_tmp );
  
  figure(5); imagesc(data_cumulants{opt_tmp.mom}); axis image; colorbar;

  
  
  full_reconst_max = zeros(szWY(1:2));
  for c1 = 1:size(rl{1},3)
    row_hat = H(c1,1); col_hat = H(c1,2);
    [inds, cut] = mat_boundary(size(full_reconst_max),row_hat-floor(opt.m/2):row_hat+floor(opt.m/2),col_hat-floor(opt.m/2):col_hat+floor(opt.m/2));
    full_reconst_max(inds{1},inds{2}) = full_reconst_max(inds{1},inds{2}) + rl{opt_tmp.mom}(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2),c1);
  end

  figure(6); imagesc(full_reconst_max); axis image; colorbar;
end

% Show ROIs extracted
figure(7); imagesc(ROI_mask); axis image; colorbar;

% Get neurofinder results
orig_regions = py.neurofinder.load([fileparts(fileparts(fileparts(get_path(opt)))) filesep 'regions' filesep 'regions.json']);
orig_ROIs = cellfun(@matpy.nparray2mat, cell(orig_regions.coordinates), 'UniformOutput', false);

orig_ROI_mask = zeros([size(ROI_mask,1), size(ROI_mask,2), 3]);
for i1 = 1:size(orig_ROIs,2)
  cur_color = 0.5*rand(1,3);
  try
    for j1 = 1:size(orig_ROIs{i1},1)
      orig_ROI_mask(orig_ROIs{i1}(j1,1), orig_ROIs{i1}(j1,2),:) = cur_color;
    end
  end
end

figure(8); imagesc(orig_ROI_mask); axis image; colorbar;
figure(9); imagesc(model.y_orig); axis image; colorbar;
