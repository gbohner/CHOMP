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
GW = cell(opt.NSS,1); % Effect of a reconstruction onto nearby projections

for obj_type = 1:opt.NSS
  WY{obj_type} = cell(opt.mom,1);
  WnormInv{obj_type} = cell(opt.mom,1);
  GW{obj_type} = cell(opt.mom,1);
  for mom1 = 1:opt.mom
    if opt.diag_tensors 
      filter_combs = opt.KS; % Fast version where we don't use combinations of filter responses for estimation
    else
      filter_combs = opt.KS^mom1; % For each object type consider all combinations of filters % TODO - exploit (super)symmetricity to reduce storage and computation
    end
    WY{obj_type}{mom1} = zeros([szY(1:2),filter_combs]); %For every location store the regression wrt all filter_combs (elements of the N_filter^mom1 projection tensor for the mom1-th moment
    WY{obj_type}{mom1} = zeros([filter_combs, filter_combs]);
    GW{obj_type}{mom1} = cell([filter_combs, filter_combs]); % This may be HUGE too big for large moments and large filter sizes - %TODO store on disk
  end
end

%% Computing WY - projections

% Compute the convolutions with the filters for each timepoint
for t1 = 1:szY(end) %TODO Can be done parallelly or on GPU 
  for obj_type = 1:opt.NSS
    conv_result = zeros([szY(1:2), length(opt.Wblocks{obj_type})]); % Result of the convolution with each filter belonging to the current object
    for filt = opt.Wblocks{obj_type} %Convolution with each filter for an object type
      Wcur = W(:,filt);
      Wcurc = Wcur(:);
      Wcurc = Wcurc./norm(Wcurc+1e-6); % make sure it has norm of 1.
      Wconv = reshape(Wcurc,opt.m,opt.m);
      conv_result(:,:,filt) = conv2(data.proc_stack.Y(:,:,t1),Wconv,'same');    
    end
    for mom1 = 1:opt.mom  %Get raw moments of the projected time course at each possible cell location %TODO - this might be wrong, because it assumes equal weighting??? but it's filter by filter, so maybe the linear combination of filters is still linear in the higher moments %TOTHINK Nah it seems correct
      all_combs = all_filter_combs(opt.KS, mom1, opt.diag_tensors); % Get the required tensors
      for i1 = 1:size(all_combs,1) % Iterate over the possible combinations (should be opt.KS^mom1)
        % Turn the rows of all_combs into a vector that counts the occurance
        % of the k-th number for easy use with existing data_structure
        mom_vec = histc(all_combs(i1,:),1:opt.KS);
        mom_vec = shiftdim(mom_vec,-1);
        WY{obj_type}{mom1}(:,:,i1) = WY{obj_type}{mom1}(:,:,i1) + prod(bsxfun(@power, conv_result, mom_vec),3);
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
    WY{obj_type} = raw2cum_multivariate(WY{obj_type});
  end
end



%% Computing WnormInv and GW - interactions between filters and projections


for obj_type = 1:opt.NSS
  for mom1 = 1:opt.mom  %Get raw moments of the projected time course at each possible cell location %TODO - this might be wrong, because it assumes equal weighting??? but it's filter by filter, so maybe the linear combination of filters is still linear in the higher moments %TOTHINK Nah it seems correct
    all_combs = all_filter_combs(opt.KS, mom1, opt.diag_tensors);
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
        

        % Computing GW - reconstruction update
        % ----------------------------------------
        % Each cell is going to be a cell of (2*m-1)^2 shifts and at each shift and
        % each moment we'll have a vector of features^moment to describe how much
        % the corresponding WY entry is modified if we set the coeffecient of active
        % filt1 at moment mom to 1.
        if exist(['./Subfuncs/Compute/Mex/computeGW.' mexext],'file') %Quicker c for loop
          GW{obj_type}{mom1}{filt1_ind,filt2_ind} = computeGW(GPT,filt1,filt2,opt.m,opt.mom,mom1);
        else %Slower Matlab for loop
          GW{obj_type}{mom1}{filt1_ind,filt2_ind} = zeros(2*opt.m-1, 2*opt.m-1); 
          for s1 = 1:(2*opt.m-1)
            for s2 = 1:(2*opt.m-1)
              GW{obj_type}{mom1}{filt1_ind,filt2_ind}(s1,s2) = filt1(GPT{s1,s2,mom1}(:,2))'* filt2(GPT{s1,s2,mom1}(:,1)); %compute the shifted effect in original space via the shift tensors GPT. Because the Worigs were computed to correspond to the best inverse of the Ws
              GW{obj_type}{mom1}{filt2_ind,filt1_ind}(s2,s1) = GW{obj_type}{mom1}{filt1_ind,filt2_ind}(s1,s2); % Use symmetricity (Swapping both filter indices AND the shift
            end
          end
        end
      end
    end
    WnormInv{obj_type}{mom1} = inv(WnormInv{obj_type}{mom1} + 1e-3 * eye(size(WnormInv{obj_type}{mom1}))); % Regularised inverse
  end
end


clearvars -except WY GW WnormInv


function all_combs = all_filter_combs(n_filt, k_choose, only_diag)
% ALL_FILTER_COMBS - Returns all required combinations given opt.KS filters
% and a mom1-combination 
  if only_diag
    all_combs = repmat((1:n_filt)',1,k_choose); % Only use [1 1 1 1], [2 2 2 2], etc [opt.KS opt.KS opt.KS opt.KS] type rows (diagonal elements of the moment tensor)
  else
    all_combs = unique(nchoosek(repmat(1:n_filt,1,n_filt), k_choose), 'rows');  
    % Returns all unique mom1-combinations of opt.KS filters (with repetition) as rows.
    % Also using this way of listing all_combs will ensure that if we
    % reshape our vector we get the correct tensor
  end
end


function filt = get_filter_comb(W, filt_comb)
  % Given a set of filters and a row of all_combs, return the requested combination
  filt = W(:, filt_comb(1));
  for i11 = 2:size(filt_comb,2)
    filt = mply(filt,  W(:, filt_comb(i11)'); % always add an extra dim by multiplying from the right with a row vector
  end
end

end

