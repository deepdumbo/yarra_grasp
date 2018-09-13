function raw_out = zeroPadPar_Filt(raw_in, nparNew)
% Zero-pad k-space data raw_in in the partition direction to npar
% partitions, using a sine filter to reduce ringing.

% Get size of raw data
[nx, nc, ntviews, npar, ne] = size(raw_in);

if nparNew>npar
    % Permute raw data so that npar is 1st, nx 2nd
    % (This makes npar*nx 'slices', that will be filtered in npar direction)
    raw_in = permute(raw_in, [4,1,3,2,5]);
    
    % Concatenate these slices to npar * gigantic matrix
    raw_in = reshape(raw_in, npar, nx*ntviews*nc*ne);
    
    % Create filter
    filter = pfFilter(size(raw_in,1));
    
    % Filter k-space in partition direction
    raw_in = raw_in .* (filter'*ones(1,size(raw_in,2)));
    
    % Zero-pad k-space in partition direction.
    raw_out = [zeros(floor((nparNew-npar)/2),size(raw_in,2));
        raw_in;
        zeros(ceil((nparNew-npar)/2),size(raw_in,2))];
    
    % Reshape & permute back
    raw_out = reshape(raw_out, nparNew, nx, ntviews, nc, ne);
    raw_out = permute(raw_out, [2,4,3,1,5]);
else
    disp('Warning: new number of partitions is not larger than old one. Skipping.');
    raw_out = raw_in;
end