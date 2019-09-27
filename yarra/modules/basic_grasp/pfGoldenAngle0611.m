function raw_out = pfGoldenAngle0611(raw_in, nparPF)
% Perform filling of missing k-space data in a golden angle 
% radial stack of stars acquisition performed with Partial Fourier in the 
% partition (slice) direction.
% For now, we consider each 'view' (ie all lines at a given angle, ie
% nx*npar points) as a separate 2D k-space and loop over those. There are
% probably (much) more efficient methods.
% 
% Input parameters:
% - raw_in:     Raw k-space data. Dimensions should be in same order as
%               read in by mapVBVD, ie [nx,nc,ntviews,npar,...]
% - nparPF:     Number of partitions after partial Fourier filling
% 
% The function properly handles all dimensions beyond npar, in that it
% simply loops over them. It therefore correctly handles things like 
% multiple echoes, averages, etc.
% 
% Written by Marnix Maas, Radboud UMC (Marnix.Maas@radboudumc.nl)
% September 2018
% 
% Uses Jeff Fesslers IRT Toolbox.

% Get data dimensions
[nx,~,~,npar,~] = size(raw_in);
if npar ~= nparPF
    % Permute dimensions:
    % from [nx,nc,ntviews,npar,...]
    %ivom: ie samples along spoke, coils, spokes, partitions.
    % to   [nx,npar,nc,ntviews,...]
    %ivom: put spokes along 2nd dim...
    ndk     = length(size(raw_in)); 
    order   = [1 4 2 3 5:ndk];                % Permutation order
    raw_in  = permute(raw_in,order);
    % Reshape raw data to [nx, npar, nc*ntviews*...]
    sr      = size(raw_in);
    raw_in  = reshape(raw_in, [nx npar prod(sr(3:end))]);
    % Loop over all nx-npar 'blades' and perform partial Fourier filling
    raw_out = single(zeros(nx, nparPF, prod(sr(3:end))));
    for i=1:prod(sr(3:end))
        [~, raw_out(:,:,i)] = ir_mri_partial_fourier_3d(raw_in(:,:,i), [nx nparPF], 'pf_location', [0 0 0], 'niter', 5);
    end
    % Reshape back to [nx,npar,nc,ntviews,...]
    raw_out = reshape(raw_out, [nx nparPF sr(3:end)]);
    % Permute back to original [nx,nc,ntviews,npar,...]
    raw_out = permute(raw_out, [1 3 4 2 5:ndk]);
    
end
