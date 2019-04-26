function [ output_args ] = wave_propagation_roi( opt )
%WAVE_PROPAGATION_ROI Summary of this function goes here
%   Detailed explanation goes here

%Load the input file we wanna use (mostly manual for now)
load(get_path(opt),'inp');
data = inp.data;
load(get_path(inp.opt,'output_iter',inp.opt.niter),'model');
opt=inp.opt;

mean_img = mean(single(data.raw_stack.Y(:,:,:)),3);
var_img = var(single(data.raw_stack.Y(:,:,:)),[],3);


%% Getting ROIs
[H, W, X, y_orig, y, L, V] = model.get_fields( 'H', 'W', 'X', 'y_orig','y','L', 'V');
if opt.fig
  update_visualize( y,H,reshape(W,opt.m,opt.m,size(W,2)),opt,1,0);
end

%% Start from a 3x3 area around the center points H and add pixels that are sufficiently similar to the included pixels

for h1 = 1:size(H,1)
  Mask = zeros(size(y));  
  row = H(i1, 1); col = H(i1, 2); type = H(i1, 3);
  tmp1, tmp2 = meshgrid(row-1:row+1, col-1:col+1);
  
  pix_lin_coord = sub2ind(size(y), row, col);
  
  Mask(pix_lin_coord)=1;
  
  neighbors = get_4neigh_idxs(pix_lin_coord, size(y));
  
  % Remove pixels that are already 1
  neighbors(Mask(neighbors)) = [];
  
  % Check which neighbors to keep:
  
  
  

end


function out = get_4neigh_idxs(pix_lin_coord, sz)
  % Return all valid 4-neighborhood linear indices in a 2D sz-sized array
  % for all coordinates in pix_lin_coord
 
  neighbor_offsets = [-1, sz(2), 1, -sz(2)];
  
  out = bsxfun(@plus, pix_lin_coord, neighbor_offsets);
  out = out(:);
  
  % TODO: Check for boundaries (atm not an issue due to padding mask)
end
  
  
function neighbor_cost_func(pix_to_check, row, col, mean_img, var_img, Mask)
  % Check if the current pixel is sufficiently different from the average
  % properties of the included pixels
  
  % Mean difference cost (how many sigmas away)
  mdc = abs(mean_img(pix_to_check) - mean(mean_img(Mask(:))))./std(mean_img(Mask(:)));
  
  % Var difference cost (how many sigmas away)
  mdc = abs(var_img(pix_to_check) - mean(var_img(Mask(:))))./std(var_img(Mask(:)));
  
  % Distance from center cost
  
  %TODO
end
    

