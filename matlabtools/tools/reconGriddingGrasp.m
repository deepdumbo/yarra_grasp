function out_img = reconGriddingGrasp(kdata)
% Function for performing a reconstruction of single-echo golden-angle 
% radial sparse (GRASP) k-space data using non-cartesian Fourier 
% transformations ('gridding' reconstruction). No binning in different
% motion or temporal states is performed.
% NOTE: it is assumed that a fourier transform along the slice direction
% has already been performed!
% Input parameters:
% kdata: contains the raw k-space data. Must have dimensions
% [nx,nc,ntviews,nz], which is the same order they are read in by mapVBVD.
% Adapted from the demo file 'iterative_multicoil.m' by Tobias Block
% by Marnix Maas (Marnix.Maas@radboudumc.nl), August 2018

% Get data dimensions
[nx,nc,ntviews,nz]=size(kdata);

% Generate trajectory
Traj=Trajectory_GoldenAngle(ntviews,nx,1,[0 0]);

% Calculate density compensation (DCF)
dcf = dcfGridding(ntviews, nx);

% ## Prepare the NUFFT
fftSize=[nx,nx];
%         (k,w,phase,shift,imSize, mode)
FT = NUFFT(Traj(:,1:end), 1, 1, [0,0], fftSize, 2);

% Create image to sum the gridded channels (for sum-of-squares combination)
slices = 1:nz;
out_img=single(zeros(nx,nx,length(slices)));
for sl=slices
    fprintf('Starting Gridding recon of slice %d / %d\n', (sl-slices(1)+1), length(slices));
    % Loop over channels and calculate image for each coil
    for channel=1:nc
        % Fetch rawdata for slice and channel and multiply with density
        % compensation function
        workKData=dcf.*double(squeeze(kdata(:,channel,:,sl-(slices(1)-1))));
        
        % Run the NUFFT
        workImage=FT'*(workKData);
        
        % Add squared channel image to buffer
        out_img(:,:,sl-(slices(1)-1)) = out_img(:,:,sl-(slices(1)-1)) + abs(workImage.*workImage);
    end
end
% Calculate the root (as final part of the sum-of-squares calculation)
out_img = sqrt(out_img);
out_img = out_img./max(out_img(:));