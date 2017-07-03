function [all_combs, mom_combs, comb_inds, comb_inds_rev] = all_filter_combs(n_filt, k_choose, only_diag)
% ALL_FILTER_COMBS - Returns all required combinations given opt.KS filters
% and a mom1-combination 
  if only_diag
    all_combs = repmat((1:n_filt)',1,k_choose); % Only use [1 1 1 1], [2 2 2 2], etc [opt.KS opt.KS opt.KS opt.KS] type rows (diagonal elements of the moment tensor)
  else
    all_combs = unique(nchoosek(repmat(1:n_filt,1,k_choose), k_choose), 'rows');  
    % Returns all unique mom1-combinations of opt.KS filters (with repetition) as rows.
    % Also using this way of listing all_combs will ensure that if we
    % reshape our vector we get the correct tensor
  end
    
  % Find the moment combinations given the filter combinations (i.e. aggregating the numbers picked multiple times);
  mom_combs = zeros(size(all_combs,1), n_filt); 
  
  for i11 = 1:size(all_combs,1)
    mom_combs(i11,:) = histc(all_combs(i11,:),1:n_filt);
  end

  % Find the unique rows of mom_combs, exploiting the known supersymmetry of
  % our tensors
  [mom_combs, comb_inds, comb_inds_rev] = unique(mom_combs,'rows','stable');
end
