function timeseries = get_cell_timeseries(opt)

%% Loading
%Add Dan's plotting tools if on neurofast and can uncomment the commented
%plots
% addpath('/home/djoshea/code/rig1/analysis/djoshea')
% addpath('/home/djoshea/code/rig1/analysis/djoshea/utils')
% import Regress.plotTuningColorGuide


%Load the input file we wanna use (mostly manual for now)
load(get_path(opt),'inp');
data = inp.data;
load(get_path(inp.opt,'output_iter',inp.opt.niter),'model');
opt=inp.opt;



%% Getting ROIs
[H, W, X, y_orig, y, L] = model.get_fields( 'H', 'W', 'X', 'y_orig','y','L');
if opt.fig
  update_visualize( y_orig,H,reshape(W,opt.m,opt.m,size(W,2)),opt,1,0);
end
%opt.ROI_type = 'quantile_dynamic_origsize';
[ROI_mask, ROIs] = getROIs(opt);

%% Showing ROIs
h_roi_figure= figure;
if ~opt.fig, set(h_roi_figure,'Visible','off'); end 
% subplottight(2, 1, 2);
h_rois = imagesc(model.y); colormap gray; axis image;
hold on;
mycolor = 'rymg';
for i1 = [1, 80, 243, 180, 200, 262, 16]%1:400%1:min(numel(ROIs),50)
  to_draw = bwboundaries(ROIs{i1}.mask);
  for c1 = 1:numel(to_draw)
    to_draw{c1}(:,2) = to_draw{c1}(:,2)+ROIs{i1}.col-floor(opt.m/2);
    to_draw{c1}(:,1) = to_draw{c1}(:,1)+ROIs{i1}.row-floor(opt.m/2);
    line(to_draw{c1}(:,2), to_draw{c1}(:,1), 'Color', mycolor(mod(ROIs{i1}.type-1,length(mycolor))+1),'LineWidth',2);
  end
  text(ROIs{i1}.col-9, ROIs{i1}.row-9, num2str(i1), 'Color', mycolor(mod(ROIs{i1}.type-1,length(mycolor))+1),'FontSize',20,'FontWeight','bold');
end

set(gca, 'XTick', []);
set(gca, 'YTick', []);

print(gcf,get_path(opt,'results_image'),'-dpng')


%    
% subplottight(2, 1, 1);
% imagesc(out_im); axis image
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% axes('Position', [0.85 0.85 0.15 0.15]);
% plotTuningColorGuide();
 
pause(0.3);
% 
% if getRandom
%   print(gcf,'cur_ROIs_rand.eps','-depsc2')
% else
%   print(gcf,'cur_ROIs.eps','-depsc2')
%   print(gcf,'cur_ROIs.png','-dpng')
% end


%% Getting timeseries from ROIs
szY = chomp_size(data.raw_stack, 'Y');
timeseries = zeros(numel(H), szY(3));

% Uses a ton of memory
patches = get_patch(data.raw_stack,opt,H,1:szY(3),'scaled',1);

for i1 = 1:numel(H)
  % % Uses low memory but more file access
  % patches = get_patch(data.raw_stack,opt,H(i1),1:szY(3),'scaled',0);
  timeseries(i1,:) = mply(ROIs{i1}.mask, patches(:,:,:,i1),2)./sum(ROIs{i1}.mask(:));
end

save(get_path(opt,'results'),'timeseries','ROIs','ROI_mask','patches','model','opt', '-v7.3');

%% Just plotting

%load(get_path(opt,'results'),'timeseries')

if 1
  figure; 
  to_plot = [1:min(size(timeseries,1),50)];%[15:20];%[10:15]+30;
  to_plot = [1:10];
  v = max(std(timeseries(to_plot,:),1,2))*2;
  %v = 5; 
  yticks = []; yticklabels = {};
  to_plot_count = 0;
  for i1 = to_plot
    to_plot_count = to_plot_count +1;
    %plot(zscore(timeseries(i1,:)) + numel(to_plot)*v - to_plot_count*v, 'LineWidth', 0.5); hold on;
    plot(timeseries(i1,:)-mean(timeseries(i1,:)) + numel(to_plot)*v - to_plot_count*v, 'LineWidth', 1.5); hold on;
    yticks = [yticks, numel(to_plot)*v - to_plot_count*v];
    yticklabels{end+1} = sprintf('Cell %3d',i1);
    xlabel('Frame')
  end
  set(gca,'YTick',yticks(end:-1:1))
  set(gca,'YTickLabel',yticklabels(end:-1:1))
  set(gca, 'FontSize', 16)
end
% %% Plotting the PCs of the timeseries
% [coeff,score,latent] = pca(timeseries,'Algorithm', 'svd', 'Centered',true);
% %[coeff,score,latent] = pca(timeseries);
% [u,s,v]=  svd(bsxfun(@times, bsxfun(@minus, timeseries, mean(timeseries,1)),1./std(timeseries,[],1))');
% %[u,s,v]=  svd(bsxfun(@minus, timeseries, mean(timeseries,1))');
% [u,s,v] = svd(timeseries');
% %[u,s,v1] = eig( zscore(timeseries)' * zscore(timeseries));
% coeff = coeff';
% figure;
% to_plot = 1:5;
% v = max(std(coeff(to_plot,:),1))*3;
% for i1 = to_plot
%   plot(coeff(i1,:) + numel(to_plot)*v - i1*v); hold on;
%   set(gca,'YTick',[])
%   xlabel('Frame')
% end

end