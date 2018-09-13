function raw_out = pfGoldenAngleMe(raw_in, twix)
% Perform filling of missing k-space data in a multi-echo golden angle 
% radial stack of stars acquisition performed with Partial Fourier in the 
% partition (slice) direction.
% For now, we consider each 'view' (ie all lines at a given angle, ie
% nx*npar points) as a separate 2D k-space and loop over those. There are
% probably (much) more efficient methods.
% 
% Written by Marnix Maas, Radboud UMC (Marnix.Maas@radboudumc.nl)
% September 2018
% 
% Uses Jeff Fesslers IRT Toolbox.

% Get size of raw data
[nx, nc, ntviews, npar, ne] = size(raw_in);

% Get number of k-space lines in partition direction that we should have
% after partial fourier filling
nparPF = twix.hdr.Meas.Partitions;

if npar ~= nparPF
    % Fill missing k-space lines using POCS algorithm:
    % Permute raw data to [nx, npar, ntviews, nc, ne]
    raw_in = permute(raw_in, [1 4 3 2 5]);
    % Reshape raw data to [nx, npar, ntviews*nc*ne]
    raw_in = reshape(raw_in, [nx npar ntviews*nc*ne]);
    % Loop over all nx-npar 'blades' and perform partial Fourier filling
    raw_out = single(zeros(nx, nparPF, ntviews*nc*ne));
    for i=1:ntviews*nc*ne
        [~, raw_out(:,:,i)] = ir_mri_partial_fourier_3d(raw_in(:,:,i), [nx nparPF], 'pf_location', [0 0 0], 'niter', 5);
    end
    % Reshape back to [nx, npar, ntviews, nc, ne]
    raw_out = reshape(raw_out, [nx nparPF ntviews nc ne]);
    % Permute back to original [nx, nc, ntviews, npar, ne]
    raw_out = permute(raw_out, [1 4 3 2 5]);
    
end
