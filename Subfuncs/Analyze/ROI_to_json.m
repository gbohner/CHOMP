function save_path = ROI_to_json( opt, ROIs, varargin )
%ROI_TO_JSON Exports the CHOMP ROI format into the one required by
%neurofinder


if nargin > 2
  num_cells = varargin{1};
else
  num_cells = numel(ROIs);
end

regions = {};
for i1 = 1:num_cells
  [rows, cols] = find(ROIs{i1}.mask);
  rows = rows + ROIs{i1}.row - (size(ROIs{i1}.mask,1)+1)/2;
  cols = cols + ROIs{i1}.col - (size(ROIs{i1}.mask,2)+1)/2;
  regions{i1} = struct('id', i1, 'coordinates', [rows(:), cols(:)]);
end

save_path =  [opt.results_folder filesep 'sources_' opt.timestamp '.json'];

savejson('', regions, save_path);

% % Include all datasets for final submission
% 
% results = [
%   struct('dataset', '00.00.test', 'regions', struct('coordinates', [[0, 1]; [2, 3]])),
%   struct('dataset', '00.01.test', 'regions', struct('coordinates', [[0, 1]; [2, 3]]))
% ]
% savejson('', results, 'results.json')
% 
% 
% % for single datasets, just collect all sources as: 
% savejson('', {results.regions}, 'results_sources.json')

end

