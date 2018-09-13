function img_out = reconPfFillConj2d(k_data, res)

% Perform an inverse Fourier Transform reconstruction of a 2D k-space
% acquired with Partial Fouier, by filling up missing lines with the
% complex conjugate and mirrored version of the acquired lines.
% The phase encoding direction is assumed to be vertical, frequency
% encoding horizontal.
% k_data is the input k-data
% res is the full resolution we want to reconstruct to: default is the
% number of frequency encoding steps

if nargin<2
    res = size(k_data,2);
end

nLinesMissing = res-size(k_data,1);
kFill = circshift(conj(flipud(fliplr(k_data(2:nLinesMissing+1,:)))), [0 1]);
kRec = [k_data; kFill];

% Perform inverse Fourier transform
img_out = ifft2(fftshift(kRec));


