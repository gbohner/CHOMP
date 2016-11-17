addpath(genpath('~/Dropbox/Gatsby/Research/CHOMP_bigdata'))


load('/mnt/data/Jim2016/k04_u000_004.mat')
load('/mnt/data/Jim2016/results/visualStim-002_Cycle00001_CurrentSettings_Ch1_000001_20160318T192722_results.mat', 'timeseries','opt', 'ROIs', 'model');


trials = 1:(size(timeseries,2)./13); %Each trial was 13 frame
types = 1:9;

len = 13;
type = 9;
rep = 6;

%% Get models for each iteration
for iters = 1:10
  model_it = get_path(opt,'output_iter',iters);
  tmp = load(model_it);
  models{iters} = tmp.model;
  
end

%% Visualise stuff
for iters = 1:10;
  W = models{iters}.W;
  W = reshape(W,opt.m,opt.m,size(W,2));
  
  NSS = opt.NSS;
KS = opt.KS;
  %Just show first obs
  NSS = 1;
  W = W(:,:,1:KS);

Nmaps = size(W,3);
isfirst = zeros(1,Nmaps);
m = size(W,1);
d = (m-1)/2;
xs  = repmat(-d:d, m, 1);
ys  = xs';
rs2 = (xs.^2+ys.^2);


    sign_center = -squeeze(sign(W(d,d,:)));
  sign_center(:) = 1;
  Wi = reshape(W, m^2, Nmaps);
  nW = max(abs(Wi), [], 1);
  %             nW = sum(Wi.^2, 1).^.5;

  Wi = Wi./repmat(sign_center' .* nW, m*m,1);

  figure(1); subplot(1,10,iters)
	visualSS(Wi, 4, KS, [-1 1]); colormap('jet')

  
end

%% Create trial averaged traces
sz = size(timeseries);
trial_avg = reshape(timeseries,sz(1),len,type,rep);