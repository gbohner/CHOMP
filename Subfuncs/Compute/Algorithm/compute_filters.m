function [WY, GW, WnormInv] = compute_filters(data, W, opt )
%COMPUTE_FILTERS Computes the correlation between filter and data, plus the
% original MAP coefficients
%   Detailed explanation goes here

%Load the shift tensor GPT (we dont really want to keep it memory all the time
if ~exist(get_path(opt, 'precomputed'),'file')
  GPT = precompute_shift_tensors(opt);
else
  load(get_path(opt, 'precomputed'), 'GPT');
end

szY = chomp_size(data.proc_stack,'Y'); %size of the data tensor

%% Initialisation of outputs

WY = cell(opt.NSS,1); % Projection of data onto basis functions
WnormInv = cell(opt.NSS,1);  % Inverse Interaction between basis functions
GW = cell(opt.mom,1); % Effect of a reconstruction onto nearby projections

for obj_type = 1:opt.NSS
  WY{obj_type} = cell(opt.mom,1);
  WnormInv{obj_type} = cell(opt.mom,1);
  for mom1 = 1:opt.mom
    if opt.diag_tensors 
      filter_combs = opt.KS; % Fast version where we don't use combinations of filter responses for estimation
    else
      filter_combs = nchoosek(opt.KS+mom1-1,mom1); % For each object type consider all combinations of filters; exploiting (super)symmetricity to reduce storage and computation
    end
    WY{obj_type}{mom1} = zeros([szY(1:2),filter_combs]); %For every location store the regression wrt all filter_combs (elements of the N_filter^mom1 projection tensor for the mom1-th moment
    WnormInv{obj_type}{mom1} = zeros([filter_combs, filter_combs]);
    GW{mom1} = cell([filter_combs*opt.NSS, filter_combs*opt.NSS]); % This may be HUGE too big for large moments and large filter sizes - %TODO store on disk
  end
end


%% Computing WY - projections

% Compute the convolutions with the filters for each timepoint
for t1 = 1:szY(end) %TODO Can be done parallelly or on GPU 
  for obj_type = 1:opt.NSS
    conv_result = zeros([szY(1:2), length(opt.Wblocks{obj_type})]); % Result of the convolution with each filter belonging to the current object
    for filt = opt.Wblocks{obj_type} %Convolution with each filter for an object type
      Wconv = reshape(W(:,filt),opt.m,opt.m);
      conv_result(:,:,mod(filt-1,opt.KS)+1) = conv2(data.proc_stack.Y(:,:,t1),Wconv,'same');    
    end
    for mom1 = 1:opt.mom  %Get raw moments of the projected time course at each possible cell location %TODO - this might be wrong, because it assumes equal weighting??? but it's filter by filter, so maybe the linear combination of filters is still linear in the higher moments %TOTHINK Nah it seems correct
      [all_combs, mom_combs] = all_filter_combs(opt.KS, mom1, opt.diag_tensors); % Get the required tensors
      
      for i1 = 1:size(mom_combs,1) % Iterate over the possible combinations (should be opt.KS^mom1)
        % Turn the rows of all_combs into a vector that counts the occurance
        % of the k-th number for easy use with existing data_structure
        WY{obj_type}{mom1}(:,:,i1) = WY{obj_type}{mom1}(:,:,i1) + prod(bsxfun(@power, conv_result, shiftdim(mom_combs(i1,:),-1)),3);
      end
    end
  end
end

% Devide by number of timesteps to get the actual moment tensors
for obj_type = 1:opt.NSS
  for mom1 = 1:opt.mom
    WY{obj_type}{mom1} = WY{obj_type}{mom1}./szY(3);
  end
end

%TODO - FORGET THIS FOR NOW
%Convert the raw moment estimates into cumulant estimates
for obj_type = 1:opt.NSS
  if opt.diag_tensors
    WY{obj_type} = raw2cum(WY{obj_type});
  else
    % TODO - This is quite inefficient in terms of access, but for proof of
    % concept it's ok
    
    % Store the tensor-collapsing and tensor-inducing index vectors (takes long if within the loop)
    comb_inds_all = cell(opt.mom,2);
    for mom1 = 1:opt.mom
      [~, ~, comb_inds, comb_inds_rev] = all_filter_combs(opt.KS, mom1, opt.diag_tensors);
      comb_inds_all{mom1, 1} = comb_inds;
      comb_inds_all{mom1, 2} = comb_inds_rev;
    end
    
    for i1 = 1:szY(1)
      for i2 = 1:szY(2)
        % For every location, convert into a cell array with moments as its
        % elements as a proper tensor
        tmp = cell(opt.mom,1);
        for mom1 = 1:opt.mom
          tmp{mom1} = reshape(WY{obj_type}{mom1}(i1,i2,comb_inds_all{mom1,2}),[opt.KS*ones(1,mom1),1]);
        end
        
        % Compute the cumulants
        tmp = raw2cum_multivariate(tmp);
        
        % Convert the resulting cumulant tensors back (into the cheap storeage of only unique elements) and store
        for mom1 = 1:opt.mom
          WY{obj_type}{mom1}(i1,i2,:) = tmp{mom1}(comb_inds_all{mom1,1});
        end
      end
    end  
  end
end



%% Computing WnormInv and GW - interactions between filters and projections


for obj_type = 1:opt.NSS
  for mom1 = 1:opt.mom  %Get raw moments of the projected time course at each possible cell location %TODO - this might be wrong, because it assumes equal weighting??? but it's filter by filter, so maybe the linear combination of filters is still linear in the higher moments %TOTHINK Nah it seems correct
    [all_combs, ~, comb_inds] = all_filter_combs(opt.KS, mom1, opt.diag_tensors);
    
    all_combs = all_combs(comb_inds,:); % Use only unique combinations - everything else is (super)symmetric - be careful to symmetrize during (re)constructions 
    
    % For each row in all_combs compute the filter 
    for filt1_ind = 1:size(all_combs,1)
      % Compute the filter tensor
      filt1 = get_filter_comb(W(:,opt.Wblocks{obj_type}), all_combs(filt1_ind,:));
      for filt2_ind = filt1_ind:size(all_combs,1)
        filt2 = get_filter_comb(W(:,opt.Wblocks{obj_type}), all_combs(filt2_ind,:));
        % TOTHINK - Do we need to renormalise these objects? Think not. Old code did it though
        
        % Computing WnormInv - reconstruction update
        % ----------------------------------------
        WnormInv{obj_type}{mom1}(filt1_ind,filt2_ind) = filt1(:)'*filt2(:);
        WnormInv{obj_type}{mom1}(filt2_ind,filt1_ind) = WnormInv{obj_type}{mom1}(filt2_ind,filt1_ind); %Use symmetricity
      end 
    end
    WnormInv{obj_type}{mom1} = inv(WnormInv{obj_type}{mom1} + 1e-3 * eye(size(WnormInv{obj_type}{mom1}))); % Regularised inverse
  end
end


%% Computing GW - reconstruction update
% ----------------------------------------
% Each cell is going to be a cell of (2*m-1)^2 shifts and at each shift and
% each moment we'll have a vector of features^moment to describe how much
% the corresponding WY entry is modified if we set the coeffecient of active
% filt1 at moment mom to 1.
%
% We need to compute the interaction between all filter responses,
% regardless of objecet type, so we need an extra loop for filt2 in
% here, and index accordingly

for mom1 = 1:opt.mom
  [all_combs, ~, comb_inds] = all_filter_combs(opt.KS, mom1, opt.diag_tensors);

  all_combs = all_combs(comb_inds,:); % Use only unique combinations - everything else is (super)symmetric - be careful to symmetrize during (re)constructions 
  szAC = size(all_combs);
      
  for filt1_ind = 1:size(GW{mom1},1)
    obj_type1 = floor((filt1_ind-1)/szAC(1))+1;
    % Compute the filter tensor
    filt1 = get_filter_comb(W(:,opt.Wblocks{obj_type1}), all_combs(mod(filt1_ind-1,szAC(1))+1,:));
    for filt2_ind = filt1_ind:size(GW{mom1},1)
      obj_type2 = floor((filt2_ind-1)/szAC(1))+1;
      filt2 = get_filter_comb(W(:,opt.Wblocks{obj_type2}), all_combs(mod(filt2_ind-1,szAC(1))+1,:));
        
      if exist(['./Subfuncs/Compute/Mex/computeGW.' mexext],'file') %Quicker c for loop
        GW{mom1}{filt1_ind,filt2_ind} = computeGW(GPT,filt1,filt2,opt.m,opt.mom,mom1);
        GW{mom1}{filt2_ind,filt1_ind} = GW{mom1}{filt1_ind,filt2_ind}';
      else %Slower Matlab for loop
        GW{mom1}{filt1_ind,filt2_ind} = zeros(2*opt.m-1, 2*opt.m-1); 
        for s1 = 1:(2*opt.m-1)
          for s2 = 1:(2*opt.m-1)
            GW{mom1}{filt1_ind,filt2_ind}(s1,s2) = filt1(GPT{s1,s2,mom1}(:,2))'* filt2(GPT{s1,s2,mom1}(:,1)); %compute the shifted effect in original space via the shift tensors GPT. Because the Worigs were computed to correspond to the best inverse of the Ws
            GW{mom1}{filt2_ind,filt1_ind}(s2,s1) = GW{mom1}{filt1_ind,filt2_ind}(s1,s2); % Use symmetricity (Swapping both filter indices AND the shift
          end
        end
      end
    end
  end
end

clearvars -except WY GW WnormInv







end

