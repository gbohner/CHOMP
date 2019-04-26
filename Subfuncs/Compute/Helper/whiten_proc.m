function data = whiten_proc( opt, data )
%WHITEN_RAW Spatio-temporally whitens a 2D+time image stack


    data_tmp =data.proc_stack.Y(:,:,:);

    % Compute a temporal mean / median for each pixel
    % averaged over large spatial areas so cells have smaller effect
    
    conv_filter_size = opt.m*5+1-mod(opt.m*5,2);
    mean_t_filtered = conv2(padarray(...
      mean(data_tmp,3),... 
      [(conv_filter_size+1)/2-1, (conv_filter_size+1)/2-1], 'replicate', 'both'),...
      fspecial('average', [conv_filter_size, conv_filter_size]), 'valid');
    
    min_t_filtered = conv2(padarray(...
      min(data_tmp,[],3),... 
      [(conv_filter_size+1)/2-1, (conv_filter_size+1)/2-1], 'replicate', 'both'),...
      fspecial('average', [conv_filter_size, conv_filter_size]), 'valid');
    
    % As we assume this estimates some affine intensity scaling, subtract a baseline value (coming from thermal noise on photo-multiplier possibly) then divide all pixel values by
    % the filtered mean, then low-pass filter to reduce noise
    min_intensity = 0; %min(mean_t_filtered(:));
    mean_t_filtered = mean_t_filtered - min_t_filtered;

    baseline_intensity = 0; %max(quantile(mean_t_filtered(:), 0.1), mean(mean_t_filtered(:))/5);
    % Ensure we do not divide with small/negative numbers);
    mean_t_filtered(mean_t_filtered < baseline_intensity) = baseline_intensity;
    
    filt_gauss_size = opt.smooth_filter_sigma;
    
    h = fir1(128, min(1./(opt.m/5), 0.85), 'low'); % Anything smaller than a 1/5th of a cell should disappear.
    h = h'*h;
    
    parfor t1 = 1:size(data.proc_stack.Y,3)
      data_tmp(:,:,t1) = (data_tmp(:,:,t1)-min_t_filtered)./mean_t_filtered;
      %data_tmp(:,:,t1) = imgaussfilt(data_tmp(:,:,t1), filt_gauss_size);
      %data_tmp(:,:,t1) = filter2(h, data_tmp(:,:,t1), 'same');
    end
   
    % Ensure that spatio-temporal global background mean is around 0
    szData = size(data_tmp);
    data_tmp = reshape(data_tmp, [prod(szData(1:2)), szData(3)]);
    
    if isempty(opt.mask_image)
        opt.mask_image = ones(szData(1:2));
    end
    
    if isempty(opt.cell_pixel_fraction)
      data_tmp = data_tmp - mean(reshape(data_tmp(opt.mask_image(:)==1,:),1,[]));
    else
      data_tmp = data_tmp - quantile(reshape(data_tmp(opt.mask_image(:)==1,:), 1, []), max(1-opt.cell_pixel_fraction, 0));
    end
    
    data.proc_stack.Y = reshape(data_tmp, szData);
    
    clear data_tmp



%     % WallisFilter version
%     tmp = data.proc_stack.Y(:,:,:);
%     windowsize = opt.m;
%     smoothsize = opt.smooth_filter_mean;
%     parfor t1 = 1:size(data.proc_stack.Y,3)
%       tmp(:,:,t1) = WallisFilter(tmp(:,:,t1), ...
%         windowsize, 0, 1, 1e10, 1, 0);
%       tmp(:,:,t1) = WallisFilter(tmp(:,:,t1), ...
%         round(windowsize/4), 0, 1, 2, 1, smoothsize);
%       tmp(:,:,t1) = WallisFilter(tmp(:,:,t1), ...
%         windowsize, 0, 1, 1e10, 1, 0);
%     end
%     data.proc_stack.Y(:,:,:) = tmp(:,:,:);


%       % Old version
%       for t1 = 1:T
%         tmp = squeeze(data.proc_stack.Y(:,:,t1)-opt.A) ./ (opt.B.^0.5);
%         data.proc_stack.Y(:,:,t1) = tmp;
%       end

end

