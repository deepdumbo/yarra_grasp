function img_out = reconPfFilt2d(k_data, res, w)

% Perform an inverse Fourier Transform reconstruction of a 2D k-space
% acquired with Partial Fouier, using zero padding to full resolution and
% sine filtering to reduce ringing artifacts (by creating a smooth
% transition between measurement & zero-padded area)
% The phase encoding direction is assumed to be vertical, frequency
% encoding horizontal.
% k_data is the input k-data
% res is the full resolution we want to reconstruct to: default is the
% number of frequency encoding steps
% w is the width of the transition band of the sine filter

if nargin<3
    w = floor(size(k_data,2)/8);
end
if nargin<2
    res = size(k_data,2);
end

% Zero-pad the k-space
kPad = zeros(size(k_data,2),res);
kPad(1:size(k_data,1),:) = k_data;

% Create a filter function:
% Initialize to 1s and 0s
filter = ones(res,1);
transitionStart = size(k_data,1) + 1 - w;
filter(transitionStart:end) = zeros(res-transitionStart+1,1);
% Make sine-shaped transition
filter(transitionStart:transitionStart+w-1) = sin( (w:-1:1) * pi / (2*w) );


% Multiply k-space with the filter
filter = filter * ones(1,size(k_data,2));
kFilt = kPad .* filter;

% Perform inverse Fourier transform
img_out = ifft2(fftshift(kFilt));
