function [W, Sv] = update_dict(datas,Hs,W,opts,k,type)
%UPDATE_DICT Summary of this function goes here
%   Detailed explanation goes here
% Use the mean image to find a reasonable backprojection of W to the
% original space, via the locations

if iscell(opts), opt=opts{1}; else opt = opts; end
  
if strcmp(opt.learn_decomp,'NMF')
  [patches, num_cells, col_count] = pick_patches(datas,Hs,opts,type,3);
elseif strcmp(opt.learn_decomp,'COV_RAW')
  [patch_cov, num_cells, col_count] = pick_patches(datas,Hs,opts,type, 2);
elseif strcmp(opt.learn_decomp,'COV')
  [patch_cov, num_cells, col_count] = pick_patches(datas,Hs,opts,type, 1);
else
  [patches, num_cells, col_count] = pick_patches(datas,Hs,opts,type,0);
end

if num_cells < k
  warning(sprintf('\nOnly %d of type %d objects have been found, returning the previous basis functions without learning\n', num_cells, type));
  W = W;
  Sv = zeros(size(W,2),1);
  return;
end

switch opt.learn_decomp
  case 'COV_RAW'
    [U,Sv,explained] = pcacov(patch_cov);
    Sv = sqrt(Sv);%*(col_count-num_cells));
  case 'COV'
    [U,Sv,explained] = pcacov(patch_cov);
    Sv = sqrt(Sv);%*(col_count-num_cells));
  case 'NMF'
    U = nnmf(reshape(patches-min(patches(:)), size(patches,1),[]), k);
    Sv = ones(size(U,2),1);
  case 'MTF'
    %TODO Maybe sometime later
  case 'HOSVD'
    [U, Sv] = svds(patches,k);
  otherwise %raise error
    error('CHOMP:learning:dict_update', 'Dictionary update option string (opt.learn_decomp) does not correspond to implemented options.');
    
end
    

% if opt.fig >1
%   figure(12);
%   plot(diag(Sv)); title('singular values of HOSVD')
% end

W(:,1:min(size(W,2),k)) = U(:,1:min(size(W,2),k));

% Transform singular values to be "usage weights"
Sv = 1./Sv(1:size(W,2));
% %Sv((min(size(W,2),k)+1):end) = 0;
% Sv = Sv.^2/sum(Sv.^2); % Get fraction of variance explained
% Sv = (1-Sv);


% % Make sure that basis functions are zerosum norms are ~1. 
% for i1 = 1:size(W,2)
%   W(:,i1) = W(:,i1) - mean(W(:,i1));
%   W(:,i1) = W(:,i1)./norm(W(:,i1)+1e-6);
%   W(:,i1) = W(:,i1) + 0.5/size(W,1);
%   W(:,i1) = W(:,i1)./norm(W(:,i1)+1e-6);
% end

%Translation step (such that center of mass is in the middle overall)
m = opt.m;
d = floor(m/2);

W = reshape(W,m,m,[]);
xs  = repmat(-d:d, m, 1);
ys  = xs';

absW = abs(W);
absW = absW/mean(absW(:));
x0 =  mean2(mean(absW,3) .* xs);
y0 =  mean2(mean(absW,3) .* ys);

xform = [1 0 0; 0 1 0; -x0 -y0 1];
tform_translate = maketform('affine',xform);

for k = 1:size(W,3)
    W(:,:,k) = imtransform(W(:,:,k), tform_translate,...
        'XData', [1 m], 'YData',   [1 m]);
%     W(:,:,k) = W(:,:,k) - mean2(W(:,:,k));
%     if std2(W(:,:,k))~=0
%       W(:,:,k) = W(:,:,k)./std2(W(:,:,k));
%     end
end

W = reshape(W,m^2,[]);

%Make sure all Ws are positive near the center
[~, mask] = transform_inds_circ(0,0,150,opt.m,max((opt.m-5)/2,1),0);
mask = logical(mask);
for k = 1:size(W,2)
  if sum(W(mask(:),k)) < sum(W(~mask(:),k))
    W(:,k) = -1.*W(:,k);
  end
end

end

