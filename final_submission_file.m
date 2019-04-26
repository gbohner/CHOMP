fname = '/nfs/data3/gergo/Neurofinder_update/final_submission_file.json';


results = [];

% For each submission file, pick the "sources_" file we want to submit

% % Dataset 00.00.test
% %cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.00.00.test/CHOMP/results/sources_20170727T055447.json');
% cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.00.00.test/CHOMP/results/sources_20170727T140331.json');
% cur_regions_struct = struct('coordinates',[]);
% for i1 = 1:length(cur_regions)
%   cur_regions_struct(i1).coordinates = cur_regions{i1}.coordinates;
% end
% results = [results, struct('dataset', '00.00.test', 'regions', cur_regions_struct)];
% 
% % Dataset 00.01.test
% %cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.00.01.test/CHOMP/results/sources_20170727T055447.json');
% cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.00.01.test/CHOMP/results/sources_20170727T140331.json');
% cur_regions_struct = struct('coordinates',[]);
% for i1 = 1:length(cur_regions)
%   cur_regions_struct(i1).coordinates = cur_regions{i1}.coordinates;
% end
% results = [results, struct('dataset', '00.01.test', 'regions', cur_regions_struct)];
% 
% % Dataset 01.00.test
% %cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.01.00.test/CHOMP/results/sources_20170727T055454.json');
% cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.01.00.test/CHOMP/results/sources_20170727T144747.json');
% cur_regions_struct = struct('coordinates',[]);
% for i1 = 1:length(cur_regions)
%   cur_regions_struct(i1).coordinates = cur_regions{i1}.coordinates;
% end
% results = [results, struct('dataset', '01.00.test', 'regions', cur_regions_struct)];
% 
% 
% % Dataset 01.01.test
% %cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.01.01.test/CHOMP/results/sources_20170727T055454.json');
% cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.01.01.test/CHOMP/results/sources_20170727T144747.json');
% cur_regions_struct = struct('coordinates',[]);
% for i1 = 1:length(cur_regions)
%   cur_regions_struct(i1).coordinates = cur_regions{i1}.coordinates;
% end
% results = [results, struct('dataset', '01.01.test', 'regions', cur_regions_struct)];



% Dataset 02.00.test
%cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.02.00.test/CHOMP/results/sources_20170727T152557.json');
%cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.02.00.test/CHOMP/results/sources_20170809T104143.json');
%cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.02.00.test/CHOMP/results/sources_20170815T172238.json');
cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.02.00.test/CHOMP/results/sources_20170816T101642.json');
cur_regions_struct = struct('coordinates',[]);
for i1 = 1:length(cur_regions)
  cur_regions_struct(i1).coordinates = cur_regions{i1}.coordinates;
end
results = [results, struct('dataset', '02.00.test', 'regions', cur_regions_struct)];


% Dataset 02.01.test
%cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.02.01.test/CHOMP/results/sources_20170727T152557.json');
%cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.02.01.test/CHOMP/results/sources_20170815T172238.json');
cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.02.01.test/CHOMP/results/sources_20170816T101642.json');
cur_regions_struct = struct('coordinates',[]);
for i1 = 1:length(cur_regions)
  cur_regions_struct(i1).coordinates = cur_regions{i1}.coordinates;
end
results = [results, struct('dataset', '02.01.test', 'regions', cur_regions_struct)];


% Dataset 03.00.test
%cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.03.00.test/CHOMP/results/sources_20170727T055501.json');
cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.03.00.test/CHOMP/results/sources_20170816T130434.json');
cur_regions_struct = struct('coordinates',[]);
for i1 = 1:length(cur_regions)
  cur_regions_struct(i1).coordinates = cur_regions{i1}.coordinates;
end
results = [results, struct('dataset', '03.00.test', 'regions', cur_regions_struct)];


% Dataset 04.00.test
%cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.04.00.test/CHOMP/results/sources_20170727T153355.json');
cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.04.00.test/CHOMP/results/sources_20170816T130548.json');
cur_regions_struct = struct('coordinates',[]);
for i1 = 1:length(cur_regions)
  cur_regions_struct(i1).coordinates = cur_regions{i1}.coordinates;
end
results = [results, struct('dataset', '04.00.test', 'regions', cur_regions_struct)];


% Dataset 04.01.test
%cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.04.01.test/CHOMP/results/sources_20170727T153355.json');
cur_regions = loadjson('/nfs/data3/gergo/Neurofinder_update/neurofinder.04.01.test/CHOMP/results/sources_20170816T130548.json');
cur_regions_struct = struct('coordinates',[]);
for i1 = 1:length(cur_regions)
  cur_regions_struct(i1).coordinates = cur_regions{i1}.coordinates;
end
results = [results, struct('dataset', '04.01.test', 'regions', cur_regions_struct)];




% Fill up the other datasets with empty regions

% missing_datasets = ...
%   {'02.00.test', '02.01.test', '04.00.test','04.01.test'};

missing_datasets = ...
  {'00.00.test', '00.01.test', '01.00.test', '01.01.test'};


cur_regions_struct = struct('coordinates', []);
cur_regions_struct(1).coordinates = [[0, 1]; [2, 3]];
cur_regions_struct(2).coordinates = [[0, 1]; [2, 3]]+5;

for i1 = 1:length(missing_datasets)
  results = [results, struct('dataset', missing_datasets{i1}, 'regions', cur_regions_struct)];
end

savejson('', results, fname);
