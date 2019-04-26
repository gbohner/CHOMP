function get_neurofinder_results( opt, ROIs )
%GET_NEUROFINDER_RESULTS Summary of this function goes here
%   Detailed explanation goes here
% Make sure that the "neurofinder" script is on Matlab's PATH
PATH = getenv('PATH');

% % Version for local laptop
% if isempty(strfind(PATH, '/Users/gergobohner/anaconda2/bin'))
%   setenv('PATH', [PATH ':/Users/gergobohner/anaconda2/bin']);
% end

% Version for office computer
if isempty(strfind(PATH, '/nfs/nhome/live/gbohner/anaconda2/bin'))
  setenv('PATH', [PATH ':/nfs/nhome/live/gbohner/anaconda2/bin']);
end
  
figure;
for num_cells = [10, 30, 50:50:numel(ROIs), numel(ROIs)]
  
    eval_command = ['neurofinder evaluate ' fileparts(fileparts(opt.data_path)) '/regions/regions.json'...
      ' ' ROI_to_json(opt, ROIs, num_cells)];
    [status,cmdout] = unix(eval_command,'-echo');

    recall = str2num(cmdout((strfind(cmdout, '"recall"') + 10):(strfind(cmdout, '"combined"') - 3)));
    precision = str2num(cmdout((strfind(cmdout, '"precision"') + 13):(strfind(cmdout, '"inclusion"') - 3)));

    scatter(num_cells, 2*recall*precision/(recall+precision)); hold on;
end

end

