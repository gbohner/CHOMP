
%% Control (red) channel picture
% Image was created in imagej by first gauss filtering (1.2 sigma) then
% averaging the stack over time, then locally contrast-normalising (10 pix,
% 1.5 std);
% model.y = imread('/Users/gergobohner/Dropbox/Gatsby/Research/forUCL/stacksForUCL/AVG_RedChan_proc.tif'); % red
% figure; imshow(model.y, [0, 22000]);
model.y = imread('/Users/gergobohner/Dropbox/Gatsby/Research/forUCL/stacksForUCL/AVG_GreenChan_proc.tif'); % green
figure; imshow(model.y);

% Show ROIs (in a zoomed in part) 1, 80, 243, 180, 200, 262, 16

%% Distributions of pixels belonging to ROIs within control (red) channel

% Run this after get_cell_timeseries, having loaded the appropriate opt and
% model files
fig_num = 132;
figure(fig_num);
close(fig_num);
h_redchan = figure(fig_num); h_orig = histogram(model.y(:), [0:500:39500, Inf], 'Normalization', 'probability', 'DisplayStyle','bar', 'LineWidth', 1);

% Get random ROIs and plot
rng(1337)
[ROI_mask, ROIs] = getROIs(opt, 1:numel(H), 1);
figure(fig_num); hold on; histogram(model.y(ROI_mask==1), h_orig.BinEdges, 'Normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 5)

% Get all ROIs and plot
[ROI_mask, ROIs] = getROIs(opt, 1:numel(H), 0);
figure(fig_num); hold on; histogram(model.y(ROI_mask==1), h_orig.BinEdges, 'Normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 3)

% Get worst 100 and plot
[ROI_mask, ROIs] = getROIs(opt, (numel(H)-99):numel(H), 0);
figure(fig_num); hold on; histogram(model.y(ROI_mask==1), h_orig.BinEdges, 'Normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 3)

% Get best 100 ROIs and plot
[ROI_mask, ROIs] = getROIs(opt, 1:100, 0);
figure(fig_num); hold on; histogram(model.y(ROI_mask==1), h_orig.BinEdges, 'Normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 5)

xlim([-1000, 4*10^4])
legend({'Raw image', 'Random ROIs', 'All ROIs', '301-400. ROIs', 'Top 100 ROIs'})

xlabel('Control channel pixel brightness');
ylabel('Probability density');



%% Spatial coefficient of variation of temporal means and variances of pixels within an ROI 

To_plot = cell(4,1); % var of means and var of vars for random rois then same for extracted ROIs

% Do it with random rois
rng(1337)
[ROI_mask, ROIs] = getROIs(opt, 1:numel(H), 1);
szRaw = size(inp.data.raw_stack.Y);
for i1 = 1:numel(ROIs)
  % Build indices
  H_cur(i1) = sub2ind(szRaw(1:2), ROIs{i1}.row, ROIs{i1}.col);
end
patches = get_patch(data.raw_stack,opt,H_cur,1:szY(3),'scaled',0);
patches = reshape(patches, [size(patches,1)*size(patches,2), size(patches,3), size(patches,4)]);

for i1=1:numel(ROIs)
  tmp = mean(patches(:,:,i1),2);
  tmp = tmp(ROIs{i1}.mask(:)==1,:);
  To_plot{1}(i1,1) = squeeze(std(tmp,[],1)./mean(tmp,1));
  tmp = std(patches(:,:,i1),[],2);
  tmp = tmp(ROIs{i1}.mask(:)==1,:);
  To_plot{2}(i1,1) = squeeze(std(tmp,[],1)./mean(tmp,1));
end


% Do it with real ROIs
[ROI_mask, ROIs] = getROIs(opt, 1:numel(H), 0);
szRaw = size(inp.data.raw_stack.Y);
for i1 = 1:numel(ROIs)
  % Build indices
  H_cur(i1) = sub2ind(szRaw(1:2), ROIs{i1}.row, ROIs{i1}.col);
end
patches = get_patch(data.raw_stack,opt,H_cur,1:szY(3),'scaled',0);
patches = reshape(patches, [size(patches,1)*size(patches,2), size(patches,3), size(patches,4)]);

for i1=1:numel(ROIs)
  tmp = mean(patches(:,:,i1),2);
  tmp = tmp(ROIs{i1}.mask(:)==1,:);
  To_plot{3}(i1,1) = squeeze(std(tmp,[],1)./mean(tmp,1));
  tmp = std(patches(:,:,i1),[],2);
  tmp = tmp(ROIs{i1}.mask(:)==1,:);
  To_plot{4}(i1,1) = squeeze(std(tmp,[],1)./mean(tmp,1));
end

mean_im = mean(data.raw_stack.Y(:,:,:),3);
tmp = mean_im(:);
mean_thr = squeeze(std(tmp,[],1)./mean(tmp,1));

std_im = std(data.raw_stack.Y(:,:,:),[],3);
tmp = std_im(:);
std_thr = squeeze(std(tmp,[],1)./mean(tmp,1));

%%
figure; h_cv_mean = axes(); 
stairs(0:0.02:1,cumsum(hist(To_plot{1}, 0:0.02:1))/numel(To_plot{1}), 'LineWidth', 3)
hold on; line([mean_thr, mean_thr], [0, 1], 'Color',[0, 0.7, 0], 'LineWidth', 3)
hold on; stairs(0:0.02:1,cumsum(hist(To_plot{3}, 0:0.02:1))/numel(To_plot{1}), 'LineWidth', 3)
legend({'Random ROIs', 'All FoV', 'True ROIs'}, 'Location', 'SouthEast')
xlabel('Coefficient of Variation, pixel temporal means')
% xlim([0,1.7])
% ylim([0,1.2])
ylabel('Cumulative density')
set(h_cv_mean, 'FontSize', 16)
figure; h_cv_std = axes();
stairs(0:0.02:1,cumsum(hist(To_plot{2}, 0:0.02:1))/numel(To_plot{1}), 'LineWidth', 3)
hold on; line([std_thr, std_thr], [0, 1], 'Color',[0, 0.7, 0], 'LineWidth', 3)
hold on; stairs(0:0.02:1,cumsum(hist(To_plot{4}, 0:0.02:1))/numel(To_plot{1}), 'LineWidth', 3)
legend({'Random ROIs', 'All FoV', 'True ROIs'}, 'Location', 'SouthEast')
xlabel('Coefficient of Variation, pixel temporal stds')
ylabel('Cumulative density')
%hold on; histogram(To_plot{4}(1:100),[0:0.02:0.95, Inf], 'Normalization', 'cdf', 'DisplayStyle','stairs', 'LineWidth', 3); 
set(h_cv_std, 'FontSize', 16)





