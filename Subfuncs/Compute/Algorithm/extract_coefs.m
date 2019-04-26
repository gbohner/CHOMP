function [ H, X, L] = extract_coefs( WY, GW, WnormInv, W, opt, varargin)
%EXTRACT_COV_COEFS Summary of this function goes here
%   Detailed explanation goes here

% opt.fig = 1; %TMP
if opt.fig > 1
%   h_dl = figure(7);  
%   h_dl2 = figure(8);  
%   h_dl3 = figure(9);
  if opt.fig > 2
    h_comb = figure(107);
    set(h_comb,'Units','Normalized')
    set(h_comb,'Position',[0.1,0.1,0.7,0.8]);
  end
  if opt.fig >3
%     Video_dl = VideoWriter('likelihood_map.avi');
%     Video_dl.FrameRate = 2;
%     open(Video_dl);
%     Video_dl2 = VideoWriter('mean_score_map.avi');
%     Video_dl2.FrameRate = 2;
%     open(Video_dl2);
%     Video_dl3 = VideoWriter('var_score_map.avi');
%     Video_dl3.FrameRate = 2;
%     open(Video_dl3);
    Video_comb = VideoWriter('inference_video.avi');
    Video_comb.FrameRate = 2;
    open(Video_comb);
  end
end

Ntypes = opt.NSS;
szWY = size(WY{1}{1}); % DataWidth * DataHeight * n_filter

n_regressors = zeros(opt.mom,1);
for mom1 = 1:opt.mom
  n_regressors(mom1) = size(WY{1}{mom1},3);
end

H = zeros(opt.cells_per_image,3); %Location (row, col, type)
X = zeros(opt.cells_per_image, sum(n_regressors)); % Basis function coefficients
L = zeros(opt.cells_per_image, opt.mom); % Likelihood gains


if opt.mask
  Mask = varargin{1};
else
  Mask = ones(szWY(1:2)); % Possible cell placements (no overlap / nearby cells);
end

Mask(1:floor(opt.m/4),:) = 0; Mask(end-floor(opt.m/4):end,:) = 0; Mask(:, 1:floor(opt.m/4)) = 0; Mask(:, end-floor(opt.m/4):end) = 0;%Don't allow cells near edges
%Mask(1:floor(opt.m),:) = 0; Mask(end-floor(opt.m):end,:) = 0; Mask(:, 1:floor(opt.m)) = 0; Mask(:, end-floor(opt.m):end) = 0;%Don't allow cells near edges
Mask = double(Mask);
%TODO: Maybe modify it such that mask is per object type

%% Initialize coefficients and likelihood maps

WYorig = WY; % TODO - Store only for testing purposes

dL_mom = zeros([szWY(1:2),Ntypes,opt.mom]);
dL = zeros([szWY(1:2), Ntypes]); % delta log likelihood
xk = cell(opt.NSS,1); % Coefficients for image filter reconstruction
for obj_type=1:opt.NSS
  xk{obj_type} = cell(opt.mom,1);
%   for mom1 = 1:opt.mom
%     xk{obj_type}{mom1} = zeros(size(WY{obj_type}{mom1}));
%   end
end

%% Locate cells 

for j = 1:opt.cells_per_image
  
  % Compute delta log-likelihoods
  if j == 1
    % Compute filter coefficients (MAP estimate)  - we'll update it in small areas
    for obj_type=1:opt.NSS
      for mom1 = 1:opt.mom
        xk{obj_type}{mom1} = ...
          mply(WY{obj_type}{mom1}, WnormInv{obj_type}{mom1}, 1);
        
        % xk{obj_type}{mom1}(xk{obj_type}{mom1}<0) = 0; %TOTHINK - NMF style constraint
      end
    end
  end
  
  % xk(xk<0) = 0; %TOTHINK
  
  %Compute delta log-likelihood
  for obj_type=1:opt.NSS
    for mom1 = 1:opt.mom
      
      if isempty(opt.W_weights)
        dL_mom(:,:,obj_type,mom1) = ...
          -sum(WY{obj_type}{mom1} .* xk{obj_type}{mom1},3);
      else % Compute penalised (by basis function weights) likelihood
        Wnorm = inv(WnormInv{obj_type}{mom1});
        dL_mom(:,:,obj_type,mom1) = ...
          -sum(WY{obj_type}{mom1} .* ...
            mply(WY{obj_type}{mom1}, ...
              inv( Wnorm + mean(diag(Wnorm))*diag(opt.W_weights(opt.Wblocks{obj_type}).^2)), ...
            1), 3);
      end
    end
  end
  
  % Renormalise
  for obj_type=1:opt.NSS
    mom_weights = ones(opt.mom,1);
    for mom1 = 1:opt.mom
        % TOTHINK: Give relative weight to the moments based on how many elements
        %they involve
        
%         % Use the inverse norm as weighting for the summation;
%         mom_weights(mom1,1) = 1./norm(reshape(dL_mom(:,:,obj_type,mom1),1,[])); 
        
        % TOTHINK: compute the likelihood change not just based on how much
        % the likelihood improve, but also give a penalty for using certain
        % basis functions given their singular value during learning step
        % dL_mom(:,:,obj_type,mom) = - sum(WY(:,:,opt.Wblocks{type},mom) .* xk(:,:,opt.Wblocks{type},mom),3);
  %       if mom>=2
  %         dL_mom(:,:,mom) = dL_mom(:,:,mom)./abs(mean2(dL_mom(:,:,mom))); %normalize the moment-related discrepencies
  %       end
       szdL_mom = [size(dL_mom), 1, 1];
       dL_mom = reshape(dL_mom, [prod(szdL_mom(1:2)), szdL_mom(3), szdL_mom(4)]);
       dL_mom(Mask(:)>0,obj_type,mom1) = ...
         zscore(dL_mom(Mask(:)>0,obj_type,mom1),[], 1);
       dL_mom = reshape(dL_mom, szdL_mom);
    end
    
    % If the user supplied extra weighting, apply that as well:
    if ~isempty(opt.mom_weights)
      mom_weights = mom_weights(:) .* opt.mom_weights(:);
    end
    
    %linear sum of individual zscored differences %TODO GMM version of "zscoring jointly"

    dL(:,:,obj_type) = mply(dL_mom(:,:,obj_type,:), mom_weights, 1); 
%     dL = dL_mom(:,:,obj_type,4); % Make it kurtosis pursuit
  end

  % Find maximum decrease  
  [AbsMin, ind] = min( dL(:).*repmat(Mask(:),Ntypes,1) );
  [row_hat, col_hat, type_hat] = ind2sub(size(dL),ind); 
    
  %Check if there is not enough likelihood decrease anymore
%   if AbsMin >= 0
%     H = H(1:(j-1), :);
%     break;
%   end    
  
  %Store the values
  H(j, :) = [row_hat, col_hat, type_hat]; % Estimated location and type
  for mom1 = 1:opt.mom
    X(j, sum([0; n_regressors(1:mom1-1)])+1:sum(n_regressors(1:mom1))) = ...
      squeeze(xk{type_hat}{mom1}(row_hat,col_hat,:));
  end
  L(j,:) = squeeze(dL_mom(row_hat,col_hat,type_hat,:)); % Estimated likelihood gain per moment
    
    
  if opt.fig >1
%     set(0,'CurrentFigure',h_dl); imagesc(dL_mom(:,:,1)); colorbar; pause(0.05);
%     set(0,'CurrentFigure',h_dl); imagesc(dL(:,:,1).*Mask); colorbar; axis square;  pause(0.05);
%     set(0,'CurrentFigure',h_dl2); imagesc(dL_mom(:,:,1,min(1,size(dL_mom,4))).*Mask); colorbar; axis square; pause(0.05);
%     set(0,'CurrentFigure',h_dl3); imagesc(dL_mom(:,:,1,min(2,size(dL_mom,4))).*Mask); colorbar; axis square; pause(0.05);
  
    if opt.fig >2
      
      [~,rl,~] = reconstruct_cell( opt, W, X(1:j,:));
      
      full_reconst = zeros(szWY(1:2));
      for c1 = 1:size(rl{1},3)
        row_hat = H(c1,1); col_hat = H(c1,2);
        [inds, cut] = mat_boundary(size(full_reconst),row_hat-floor(opt.m/2):row_hat+floor(opt.m/2),col_hat-floor(opt.m/2):col_hat+floor(opt.m/2));
        full_reconst(inds{1},inds{2}) = full_reconst(inds{1},inds{2}) + rl{1}(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2),c1);
      end
     %load(get_path(opt,'output_iter',4),'model')
     %set(0,'CurrentFigure',h_comb);
     %subplot(2,2,1); imagesc(model.y_orig); colormap gray; axis square; axis off; title('Mean Data');
     subplot(2,2,1); imagesc(-dL(:,:,1)); colormap default; colorbar; axis square; axis off; title('Log likelihood increase');  pause(0.0001);
     subplot(2,2,3); imagesc(full_reconst); colorbar; axis square; axis off; title('Reconstructed cells'); pause(0.0001); 
     subplot(2,2,2); imagesc(-dL_mom(:,:,min(1,size(dL_mom,3)),min(1,size(dL_mom,4)))); colorbar; axis square; axis off; title('mom=1 LL change'); pause(0.0001);
     subplot(2,2,4); imagesc(-dL_mom(:,:,min(1,size(dL_mom,3)),min(4,size(dL_mom,4)))); colorbar; axis square; axis off; title(sprintf('mom=%s LL change',size(dL_mom,4))); pause(0.0005);
%    set(0,'CurrentFigure',h_dl3); imagesc(Mask(:,:,1)); colorbar; axis square; pause(0.05);
    end
  end
  
  
  
  
  %Affected local area
  % Size(,1) : number of rows, size(,2): number of columns
 [inds, cut] = mat_boundary(szWY(1:2),(row_hat-opt.m+1):(row_hat+opt.m-1),(col_hat-opt.m+1):(col_hat+opt.m-1));
  
 
 %% Update WY
 % Compute reconstruction
 if j==1 % Generate Wfull for each object type
   load(get_path(opt, 'precomputed'), 'GPT');
   Wfull = cell(opt.NSS,1);
   for obj_type = 1:opt.NSS
     [~,~,Wfull_cur] = reconstruct_cell( opt, W(:,opt.Wblocks{type_hat}), X(j,:));
     Wfull{obj_type} = Wfull_cur;
   end
 end
 
 reconst = reconstruct_cell( opt, W(:,opt.Wblocks{type_hat}), X(j,:),'Wfull',Wfull{type_hat});
 
 % Compute change in WY (by computing the projection of order-mom
 % reconstruction onto each order-mom filter
 WYchange = cell(opt.NSS,1);
 for obj_type = 1:opt.NSS
   WYchange{obj_type} = cell(opt.mom,1);
   for mom1 = 1:opt.mom
%      reconst_cur = reshape(reconst{mom1},opt.m, opt.m, []); % m x m x REST tensor
     WYchange{obj_type}{mom1} = zeros(2*opt.m-1,2*opt.m-1,size(Wfull{obj_type}{mom1},2));
     for i1 = 1:size(Wfull{obj_type}{mom1},2) % Iterate through each filter
       for s1 = 1:2*opt.m-1
         for s2 = 1:2*opt.m-1
           WYchange{obj_type}{mom1}(s1,s2,i1) = ...
            reconst{mom1}(GPT{s1,s2,mom1}(:,2))'*Wfull{obj_type}{mom1}(GPT{s1,s2,mom1}(:,1),i1);
         end
       end
%        W_cur = reshape(Wfull{obj_type}{mom1}(:,i1),[opt.m*ones(1,2*mom1),1]); % m x m x ... x m tensor (2*mom1 dimensions)
%        for dim1 = 1:ndims(W_cur)
%          W_cur = flip(W_cur,dim1); % flip all dimensions to work correctly with convolution
%        end
%        W_cur = reshape(W_cur,opt.m, opt.m, []); % m x m x REST tensor BUT with dimensions flipped accordingly for convolution
%        for dim1 = 1:size(reconst_cur,3) % Iterate through each layer of the tensor (cause we need correlation along shifts in the first 2 dimensions, the rest stays)
%          WYchange{obj_type}{mom1}(:,:,i1) = ...
%            WYchange{obj_type}{mom1}(:,:,i1) + conv2(reconst_cur(:,:,dim1), W_cur(:,:,dim1),'full');
%        end
     end
   
   % Update WY
   WY{obj_type}{mom1}(inds{1},inds{2},:) = ...
     WY{obj_type}{mom1}(inds{1},inds{2},:) - WYchange{obj_type}{mom1}((1+cut(1,1)):(end-cut(1,2)),(1+cut(2,1)):(end-cut(2,2)),:);
  end
 end
 
   
   
   
 
%% WY UPDATE OLD (FASTER?) VERSION - SEEMS TO BE A BUG WITH GW ?
%   % Compute the changes in WY (the effect of the saved filter on nearby
%   % locations for all object types)
%   for mom1 = 1:opt.mom
%     for filt1_ind = 1:size(WY{type_hat}{mom1},3)
%       for obj_type = 1:opt.NSS
%         for filt2_ind = 1:size(WY{obj_type}{mom1},3)
%           WY{obj_type}{mom1}(inds{1},inds{2},filt2_ind) = ... 
%             WY{obj_type}{mom1}(inds{1},inds{2},filt2_ind) - ...
%             ( ...
%               GW{mom1}{filt1_ind,filt2_ind}(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2)) * ... % The large interaction tensor
%               X(j, sum([0; n_regressors(1:mom1-1)])+filt1_ind) ... % The stored filter coefficient
%             );
%         end
%       end
%     end
%   end
%  
  
%% Update the changed xk values and delta log-likelihoods
  for obj_type=1:opt.NSS
    for mom1 = 1:opt.mom
      xk{obj_type}{mom1}(inds{1},inds{2},:) = ...
        mply(WY{obj_type}{mom1}(inds{1},inds{2},:), WnormInv{obj_type}{mom1}, 1);
      
%       % TOTHINK - Remove negative values
%       xk{obj_type}{mom1}(inds{1},inds{2},:) = ...
%         xk{obj_type}{mom1}(inds{1},inds{2},:) .* (xk{obj_type}{mom1}(inds{1},inds{2},:)>0);

%       dL_mom(inds{1},inds{2},obj_type,mom1) = ...
%         -sum(WY{obj_type}{mom1}(inds{1},inds{2},:) .* xk{obj_type}{mom1}(inds{1},inds{2},:),3);
    end
  end
    
  
%    figure(4); imagesc(WY(:,:,1)); colorbar; pause(0.05);  
 
  % Update the patch around the point found
%   Mask(max(row-3,1):min(row+3,end),max(col-3,1):min(col+3,end),type) = 0; % Make it impossible to put cells to close to eachother
%   Mask(max(row-1,1):min(row+1,end),max(col-1,1):min(col+1,end),:) = 0; % Make it impossible to put cells to close to eachother
    Mask(row_hat,col_hat,:) = 0; %Make it impossible to put a cell into the exact same location
  
if ~isempty(opt.spatial_push)
  [yinds, ycut] = mat_boundary(szWY(1:2),row_hat-opt.m:row_hat+opt.m,col_hat-opt.m:col_hat+opt.m);  
  [gridx,gridy] = meshgrid((-opt.m+ycut(1,1)):(opt.m-ycut(1,2)),(-opt.m+ycut(2,1)):(opt.m-ycut(2,2))); %make sure to cut the corresponding dimensions using ycut
  gridx = gridx'; gridy = gridy'; %meshgrid(1:n,1:m) creates mxn matrices, need to transpose
  
  grid_dist = sqrt(gridx.^2+gridy.^2);
  grid_dist = opt.spatial_push(grid_dist); % Specified distance based function
  Mask(yinds{1},yinds{2},:) = Mask(yinds{1},yinds{2},:).*repmat(grid_dist,[1,1,size(Mask,3)]); % Make it impossible to put cells to close to eachother
end
  
 if opt.fig >3
%   writeVideo(Video_dl, getframe(h_dl));
%   writeVideo(Video_dl2, getframe(h_dl2));
%   writeVideo(Video_dl3, getframe(h_dl3));
  writeVideo(Video_comb, getframe(h_comb));
 end  
%   disp([num2str(j) ' cells found, current type: ' num2str(type)]);
end

  if opt.fig >3 
%     close(Video_dl);
%     close(Video_dl2);
%     close(Video_dl3);
    close(Video_comb);
  end
% close(Video_yres);
end

