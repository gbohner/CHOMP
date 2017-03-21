function [ ROI_image, ROIs ] = getROIs( opt, varargin )
%GETROIS Summary of this function goes here
%   Detailed explanation goes here

load(get_path(opt, 'output_iter', opt.niter), 'model');
 [H, W, X, y_orig, y] = model.get_fields( 'H', 'W', 'X', 'y_orig','y');

%update_visualize( y_orig,H,reshape(W,opt.m,opt.m,size(W,2)),opt,1);

if nargin>1
  to_reconst = varargin{1};
  num_reconst = max(varargin{1});
else
  to_reconst = 1:size(H,1);
  num_reconst = size(H,1);
end

if nargin>2
  if varargin{2}
    %Get random ROIs
    H = floor(1+rand(size(H)).*numel(y));
    X = randn(size(X));
  end
end

sz = size(y);

ROIs = cell(num_reconst,1);
switch opt.ROI_type
    case 'quantile'
      ROI_image = zeros(sz(1:2));
    case 'quantile_dynamic'
      ROI_image = zeros(sz(1:2));
    case 'quantile_origsize'
      szRaw = size(y_orig);
      ROI_image = zeros(szRaw(1:2));
    case 'quantile_dynamic_origsize'
      szRaw = size(y_orig);
      ROI_image = zeros(szRaw(1:2));
    case 'mean_origsize'
      szRaw = size(y_orig);
      ROI_image = zeros(szRaw(1:2));
end


[all_reconst,all_reconst_lowdim] = reconstruct_cell( opt, W, X);

for i1 = to_reconst
  row = H(i1, 1); col = H(i1, 2); type = H(i1, 3);
%   if opt.mom>=2 %Then reconstruct variance image
%     reconst = all_reconst_lowdim{2}(:,:,i1);
%   else
    reconst = all_reconst_lowdim{1}(:,:,i1);
  
  reconst_orig = reconst;
  
  if opt.fig>=3
    figure; imagesc(reconst); pause(0.2);
  end
  
  switch opt.ROI_type
    case 'quantile'
      reconst = reconst > quantile(reconst(:),opt.ROI_params(1));
      reconst = bwconvhull(reconst);
      [inds, cut] = mat_boundary(sz(1:2),row-floor(opt.m/2):row+floor(opt.m/2),col-floor(opt.m/2):col+floor(opt.m/2));
      ROI_image(inds{1},inds{2}) = reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2));
    case 'quantile_origsize'
      cutsize = round(opt.m./opt.spatial_scale);
      cutsize = cutsize + (1-mod(cutsize,2));
      reconst = imresize(reconst,[cutsize,cutsize]);
      reconst = reconst > quantile(reconst(:),opt.ROI_params(1));
      reconst = imerode(reconst,strel('rectangle',[3,3]));
      reconst = imdilate(reconst,strel('rectangle',[3,3]));
      reconst = imdilate(reconst,strel('rectangle',[3,3]));
      reconst = imerode(reconst,strel('rectangle',[3,3]));
      reconst = imerode(reconst,strel('rectangle',[3,3]));
      %reconst = bwconvhull(reconst); %looks better for visualizations
      row = round(row./opt.spatial_scale);
      col = round(col./opt.spatial_scale);
      [inds, cut] = mat_boundary(szRaw(1:2),row-floor(cutsize/2):row+floor(cutsize/2),col-floor(cutsize/2):col+floor(cutsize/2));
      ROI_image(inds{1},inds{2}) = reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2));
    case 'mean_origsize'
      cutsize = round(opt.m./opt.spatial_scale);
      cutsize = cutsize + (1-mod(cutsize,2));
      reconst = imresize(reconst,[cutsize,cutsize]);
      %reconst = reconst > quantile(reconst(:),opt.ROI_params(1));
      %reconst = bwconvhull(reconst); %looks better for visualizations
      row = round(row./opt.spatial_scale);
      col = round(col./opt.spatial_scale);
      [inds, cut] = mat_boundary(szRaw(1:2),row-floor(cutsize/2):row+floor(cutsize/2),col-floor(cutsize/2):col+floor(cutsize/2));
%       %get the covariance between pixels of reconstruction and original
%       %mean
%       reconst = reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2)) .* y_orig(inds{1},inds{2});
      ROI_image(inds{1},inds{2}) = reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2));
    case 'quantile_dynamic_origsize'
      cutsize = round(opt.m./opt.spatial_scale);
      cutsize = cutsize + (1-mod(cutsize,2));
      reconst = imresize(reconst,[cutsize,cutsize]);
      row = round(row./opt.spatial_scale);
      col = round(col./opt.spatial_scale);
      [inds, cut] = mat_boundary(szRaw(1:2),row-floor(cutsize/2):row+floor(cutsize/2),col-floor(cutsize/2):col+floor(cutsize/2));
      tmp = opt.ROI_params(1)*(max(reshape(y_orig(inds{1},inds{2}),1,[])')-min(reshape(y_orig(inds{1},inds{2}),1,[])'))+min(reshape(y_orig(inds{1},inds{2}),1,[])');
      tmp = sum(sum(y_orig(inds{1},inds{2})>tmp))./numel(y_orig(inds{1},inds{2})); %ratio of pixels to pick
      [reconst_vals, tmp_ind] = sort(reconst(:),'descend');
      last_val = reconst_vals(floor(tmp*numel(reconst)));
      last_ind = find(reconst_vals == last_val,1,'last');
      reconst(:) = 0;
      reconst(tmp_ind(1:last_ind))=1;
      ROI_image(inds{1},inds{2}) = reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2));
    otherwise
      error('CHOMP:roi:method',  'Region of interest option string (opt.ROI_type) does not correspond to implemented options.')
  end
  
  %Store results in cells now, later add option for json output %TODO
  ROIs{i1} = struct('col', col, 'row', row, 'type', type, 'mask', reconst, 'reconst', reconst_orig);
  
  
  
  
end

  if opt.fig >1
    figure(2);
    imshow(ROI_image); axis square;
  end

end

