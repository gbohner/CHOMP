function [ out ] = symmetrise( inp )
%SYMMETRISE Symmetrises an n-dimensional hypercube by rotating through it's
%dimensions, summing up the objects then dividing by the number of
%dimensions.

out = zeros(size(inp));
for j1 = 0:(ndims(inp)-1)
  out = out + shiftdim(inp,j1);
end

out = out./ndims(inp);


end

