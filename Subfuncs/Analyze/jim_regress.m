%% Load Analyzer file and data
load('/mnt/data/Jim2016/k04_u000_004.mat')
load('/mnt/data/Jim2016/results/visualStim-002_Cycle00001_CurrentSettings_Ch1_000001_20160318T192722_results.mat', 'timeseries');

trials = 1:(size(timeseries,2)./13); %Each trial was 13 frame
types = 1:9;

len = 13;
type = 9;
rep = 6;

%% Create trial averaged traces
sz = size(timeseries);
trial_avg = reshape(timeseries,sz(1),len,type,rep);

%%
figure; 
for ty1 = 1:type
  plot(mean(trial_avg(2,:,ty1,:),4));
  colormap autumn
  hold on;
end

%%
figure; 
to_plot = [1:min(size(timeseries,1),10)];%[15:20];%[10:15]+30;
v = max(std(timeseries(to_plot,:),1))./5;
for i1 = to_plot
  plot(timeseries(i1,:) + numel(to_plot)*v - i1*v, 'LineWidth', 2); hold on;
  set(gca,'YTick',[])
  xlabel('Frame')
end
