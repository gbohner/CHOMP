function out = raw2cum( inp )
%RAW2CUM Converts raw moments to cumulants for cell array inp, 
% where raw moments are stored within the individual cells inp{1} = mean,
% inp{2} = second moment etc AND the individual entries of the matrices are
% independent

%Put the moment dimension first
out = inp;
for moms = 1:numel(out)
  for i1=1:moms-1
      out{moms,:} = out{moms} - nchoosek(moms-1,i1-1) .* out{i1}.*inp{moms-i1}; % .* instead of outer product due to independent entries 
  end
end

end

