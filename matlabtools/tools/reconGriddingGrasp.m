function out_img = reconGriddingGrasp(kdata, coilCombineMode)
% Function for performing a reconstruction of single-echo golden-angle 
% radial sparse (GRASP) k-space data using non-cartesian Fourier 
% transformations ('gridding' reconstruction). No binning in different
% motion or temporal states is performed.
% NOTE: it is assumed that a fourier transform along the slice direction
% has already been performed!
% Input parameters:
% - kdata: contains the raw k-space data. Must have dimensions
%          [nx,nc,ntviews,nz], which is the same order they are read in by 
%          mapVBVD.
% - coilCombineMode: 
%          method for combining signals from the different coil
%          elements. 0 = sum-of-squares (default), 1 = using Ricardo Otazo's
%          adapt_array_2d (unsure whether this is the same as Siemens' 'Adaptive
%          Combine' option
% Adapted from the demo file 'iterative_multicoil.m' by Tobias Block
% by Marnix Maas (Marnix.Maas@radboudumc.nl), August 2018

% Parse inputs
if nargin<2, coilCombineMode = 0; end

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
out_img=single(zeros(nx,nx,nz));
% out_img=single(zeros(nx,nx,nz,nc));
for sl=1:nz
    fprintf('Starting Gridding recon of slice %d / %d\n', sl, nz);
    % Loop over channels and calculate image for each coil
    workImage = zeros(nx,nx,nc);
    for channel=1:nc
        % Fetch rawdata for slice and channel and multiply with density
        % compensation function
        workKData=dcf.*double(squeeze(kdata(:,channel,:,sl)));
        
        % Run the NUFFT
        workImage(:,:,channel) = FT'*(workKData);
    end
    % Combine coil elements
    if coilCombineMode
        % "Adaptive combine" coil combination
        [~, out_img(:,:,sl)] = adapt_array_2d(workImage);
    else
        % Sum of squares combination
        out_img(:,:,sl) = sum(abs(workImage.*workImage),3);
    end
end

% Normalize image
out_img = out_img./max(abs(out_img(:)));