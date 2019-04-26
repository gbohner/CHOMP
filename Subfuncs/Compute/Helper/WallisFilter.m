function [G, M, D] = WallisFilter(F, W, Md, Dd, Amax, p, smoothing)
%WALLISFILTER Implements local contrast normalisation in an image
% F - image
% W - window size of "local" filtering
% Md - Desired local mean
% Dd - Desired local stanard deviation
% p - Proportionality between original and desired mean (set to 1 for only desired mean)
% Amax - Laplace smoothing factor for low variance data (maximum allowed standard deviation of a single pixel)
% smoothing - Binary variable, if nonzero, gaussian smoothing in window W with spatial stanard deviation "smoothing"
if mod(W,2)==0
    W = W+1;
end

F = double(F);
% Image of ones to count number of elements (to handle edge effects)
convNumEl = ones(size(F));


if smoothing > 0
  gauss_blur = fspecial('gaussian', W, smoothing)*(W^2);
  F = conv2(F, gauss_blur, 'same')./conv2(convNumEl,ones(W),'same');
end

% Get the mean in every window W
M = conv2(F,ones(W),'same')./conv2(convNumEl,ones(W),'same');

% Get the standard deviation in every window W
D = (conv2((F-M).^2,ones(W),'same')./conv2(convNumEl,ones(W),'same')).^0.5;

% Create the normalised image with Md mean and Dd standard deviation in
% every window
G = (F-M) .* Amax .* Dd ./ (Amax .* D + Dd) + p * Md + (1-p) * M;

end