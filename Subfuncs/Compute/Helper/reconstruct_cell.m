function [ reconst, reconst_lowdim, Wfull ] = reconstruct_cell( opt, W, X, varargin )
%RECONSTRUCT_CELL Given a basis set W and corresponding coefficients X, reconstructs the mom-th moment of the cell 

%% Define inputs
p = inputParser();
p.addRequired('opt',@(x)isa(x,'chomp_options'));
p.addRequired('W',@isnumeric);
p.addRequired('X',@isnumeric);
p.addParameter('Wfull',{}); % If given, then no need to recompute the higher order regressors
p.addParameter('mom_max',opt.mom); % Return reconstructions up to this moment
p.parse(opt,W,X,varargin{:});

Wfull = p.Results.Wfull;
mom_max = p.Results.mom_max;

%% Compute the required new set of Ws 
if isempty(Wfull)
  Wfull = cell(opt.mom,1);
  for mom1 = 1:opt.mom
    [all_combs, ~, comb_inds] = all_filter_combs(opt.KS, mom1, opt.diag_tensors);
    all_combs = all_combs(comb_inds,:); % Use only unique combinations - everything else is (super)symmetric - be careful to symmetrize during (re)constructions
    Wfull{mom1} = zeros((opt.m^2)^mom1, size(all_combs,1)); 
    
    % For each row in all_combs compute the filter
    for filt1_ind = 1:size(all_combs,1)
      % Compute the filter tensor
      filt1 = get_filter_comb(W, all_combs(filt1_ind,:));
      Wfull{mom1}(:,filt1_ind) = filt1(:);
    end
  end
else
  assert(numel(Wfull)>=mom_max, 'CHOMP:reconstruct:wrongWfullGiven');
end


%% Compute the reconstructions

reconst = cell(mom_max,1);
reconst_lowdim = cell(mom_max,1);
curXstart = 0;

for mom1 = 1:mom_max
  % Compute reconstruction as an order mom1 tensor
  reconst{mom1} = Wfull{mom1}*X(:,(curXstart+1):(curXstart+size(Wfull{mom1},2)))';
  curXstart = curXstart + size(Wfull{mom1},2);
  %reconst{mom1} = reshape(reconst{mom1},[opt.m^2*ones(1,mom1), size(reconst{mom1},2)]);
  
  reconst_lowdim{mom1} = zeros(opt.m,opt.m,size(reconst{mom1},2));
  % Reshape into 2D images of cells
  for cell_num = 1:size(reconst_lowdim{mom1},3)
    if mom1==1 
        reconst_lowdim{mom1}(:,:,cell_num) = reshape(reconst{mom1}(:,cell_num),opt.m,opt.m);
    else
      % For higher moment just take the reconstruction on the diagonal to
      % show in the original (low-dim) space
      reconst_lowdim{mom1}(:,:,cell_num) = ...
        reshape(reconst{mom1}(linspace(1,size(reconst{mom1},1),opt.m^2),cell_num),opt.m,opt.m);
    end
  end
end


end
