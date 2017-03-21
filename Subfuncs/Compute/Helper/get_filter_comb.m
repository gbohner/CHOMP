function filt = get_filter_comb(W, filt_comb)
  % Given a set of filters and a row of all_combs, return the requested combination
  filt = W(:, filt_comb(1));
  for i11 = 2:size(filt_comb,2)
    % Make sure everything we do is (super)symmetric!
    filt = symmetrise(mply(filt,  W(:, filt_comb(i11))', 0)); % always add an extra dim by multiplying from the right with a row vector
  end
end