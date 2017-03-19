function out = raw2cum_multivariate( inp )
%RAW2CUM Converts raw moments to cumulants for multivariate distributions for cell array inp, 
% where raw moments are stored within the individual cells inp{1} = mean,
% inp{2} = second moment etc

%Convert the raw moments into cumulants
out = inp;
for moms = 1:numel(inp)
  for i1=1:moms-1
      % Symmetrise what we subtract (as we know these moment and cumulant
      % tensors have to be supersymmetric)
      out{moms} = out{moms} - nchoosek(moms-1,i1-1) .* symmetrise(mply(out{i1}, shiftdim(inp{moms-i1},-1))); 
  end
end

end

