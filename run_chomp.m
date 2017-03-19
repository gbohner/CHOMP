close all;
clear;

%Example local run

cd(fileparts(mfilename('fullpath')));
addpath(genpath('.'))

setenv('CHOMP_ROOT_FOLDER', '');

opt = chomp_options(); %Initialize default options

% % Setup folder structure for local run and local data 
%opt.data_path = '~/stanford/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001_Cycle00001_Ch2_000001.ome.tif';
opt.data_path = '/Users/gergobohner/Data/stacksForUCL/GreenChan/GreenChanB_0000.tif';
% opt.data_path = '/Users/gergobohner/Dropbox/Gatsby/Research/forUCL/stacksForUCL/RedChan/RedChanD_0000.tif';

opt.data_type = 'frames';
%opt.data_path = '/mnt/gatsby/nfs/data3/gergo/Neurofinder/neurofinder.00.00';
% opt.data_path = '~/tmp/nf/neurofinder.00.00';
% opt.timestamp = '20161117T174032';

%opt.data_path = '~/Data/Neurofinder/neurofinder.01.01/images/image00001.tiff';
%opt.data_path = '~/Data/Neurofinder/neurofinder.03.00/images/image00001.tiff';
%opt.src_string = '.tiff';

opt.src_string = '.tif';

opt.m = 5;
opt.stabilize = 0;
opt.whiten = 1;
opt.spatial_scale = 0.2;
opt.time_scale = 0.1;
opt.NSS = 1;
opt.KS = 3;
opt.mom = 4;
opt.init_model = 'filled'; %{'filled', 'donut'};
opt.niter = 3;
opt.learn = 1;

% opt.m = 7;
% opt.stabilize = 0;
% opt.spatial_scale = 0.28;
% opt.time_scale = 0.2;
% opt.NSS = 1;
% opt.KS = 4;
% opt.mom = 2;
% opt.init_model = 'filled'; %{'filled', 'donut'};
% opt.niter = 5;


% opt.input_folder = '~/stanford/input/';
% opt.output_folder = '~/stanford/output/';
% opt.precomputed_folder = '~/stanford/precomputed/';
% opt.results_folder = '~/stanford/results/';
 
opt.precomputed_folder = './precomputed/';

% Set the subfolders - neurofinder specific
opt.input_folder = [fileparts(fileparts(opt.data_path)) '/CHOMP/input/'];
opt.output_folder = [fileparts(fileparts(opt.data_path)) '/CHOMP/output/'];
opt.results_folder = [fileparts(fileparts(opt.data_path)) '/CHOMP/results/'];

[~, tmp1, tmp2] = fileparts(fileparts(fileparts(opt.data_path))); % Important if you want to have nice file prefixes (corresponding to main folder name)
opt.file_prefix = [tmp1 tmp2];

opt.fig = 2;
opt.verbose = 2;


ROI_type = 'quantile_dynamic_origsize';
ROI_params = [0.4];



%%
% %Learning with 100 cells;
% opt.cells_per_image = 100;
% %opt.cells_per_image = 1000;
% 
% % Ensure no overlaps at all
% opt.spatial_push = @(grid_dist)grid_dist>=(2*opt.m -1);
% 
% %opt.spatial_push = @(grid_dist)grid_dist>=(2*opt.m - 13);
% [opt, ROI_mask, ROIs] = chomp(opt);

opt.cells_per_image = 100;
opt.spatial_push = @(grid_dist)grid_dist>=25;

[opt, ROI_mask, ROIs] = chomp(opt);

  %%
%Inferring with 500 cells;
opt.cells_per_image = 400;
opt.init_iter = opt.niter; %Just running the inference
opt.niter = opt.niter+1; %Just running the inference

% Use a less restrictive spatial push for inference
% opt.spatial_push = @(grid_dist)logsig(0.5*grid_dist-floor(opt.m/2-1)); %@(grid_dist, sharp)logsig(sharp*grid_dist-floor(sharp*2*obj.m/2-1));
%opt.spatial_push = @(grid_dist)logsig(0.5*grid_dist-floor(opt.m/2-1)).*(grid_dist>opt.m);
opt.spatial_push = @(grid_dist)grid_dist>=20;

[opt, ROI_mask, ROIs] = chomp(opt);

% 
get_cell_timeseries(opt);

%% Show evaluation results with neurofinder script

% 
% % Make sure that the "neurofinder" script is on Matlab's PATH
% PATH = getenv('PATH');
% if isempty(strfind(PATH, '/Users/gergobohner/anaconda2/bin'))
%   setenv('PATH', [PATH ':/Users/gergobohner/anaconda2/bin']);
% end
%   
% figure;
% for num_cells = [10, 30, 50:50:1000]
% 
%   eval_command = ['neurofinder evaluate ' fileparts(fileparts(opt.data_path)) '/regions/regions.json'...
%     ' ' ROI_to_json(opt, ROIs, num_cells)];
%   [status,cmdout] = unix(eval_command,'-echo');
%   
%   recall = str2num(cmdout((strfind(cmdout, '"recall"') + 10):(strfind(cmdout, '"combined"') - 3)));
%   precision = str2num(cmdout((strfind(cmdout, '"precision"') + 13):(strfind(cmdout, '"inclusion"') - 3)));
%   
%   scatter(recall, precision); hold on;
% end