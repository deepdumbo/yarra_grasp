function kRec = pfFillConj2d(k_data, res)

% Fill in missing data in 2D k-space acquired with Partial Fouier, by 
% filling up missing lines with the complex conjugate and mirrored version 
% of the acquired lines.
% The phase encoding direction is assumed to be vertical, frequency
% encoding horizontal.
% k_data is the input k-data
% res is the full resolution we want to reconstruct to: default is the
% number of frequency encoding steps
% If k_data has more than 2 dimensions, the filling is done in the first 2
% dimensions only.

if nargin<2
    res = size(k_data,2);
end

s = size(k_data);
nLinesMissing = res-s(1);

% Calculate chunk of k-space data that is to be filled in
% Note that shifting by 1,1 is necessary because of (usually) even number
% of lines and columns in k-space. This will probably mess up if odd
% dimensions are used.
kFill = circshift(conj(flip(flip(k_data(2:nLinesMissing+1,1:res,:),1),2)), [0 1]);
% Restore dimensions of filly-chunky
kFill = reshape(kFill, [size(kFill,1), s(2:end)]);
% Append it to partial k-space data
kRec = cat(1,k_data,kFill);



