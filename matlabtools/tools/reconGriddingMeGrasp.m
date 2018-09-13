function out_img = reconGriddingMeGrasp(kdata)
% Function for performing a reconstruction of multi-echo golden-angle 
% radial sparse (GRASP) k-space data using non-cartesian Fourier 
% transformations ('gridding' reconstruction). No binning in different
% motion or temporal states is performed.
% NOTE: it is assumed that a fourier transform along the slice direction
% has already been performed!
% Input parameters:
% kdata: contains the k-space data that has been Fourier transformed along z. 
% Must have dimensions [nx,nc,ntviews,nz,ne], which is the same order they
% are read in by mapVBVD.
% Adapted from the demo file 'iterative_multicoil.m' by Tobias Block
% by Marnix Maas (Marnix.Maas@radboudumc.nl), August 2018

% Get data dimensions
[nx,nc,ntviews,nz,ne]=size(kdata);

% Generate trajectory
Traj=Trajectory_GoldenAngle_ME(ntviews, nx, ne, 2, 'normalize',1,'flip',[0 0]);

% Calculate density compensation (DCF)
dcf = dcfGridding(ntviews, nx);

% Set parameters
fftSize=[nx,nx];
% slices = 1:nz;

% Pre-allocate output image
out_img=single(zeros(nx,nx,nz,ne));

for ec = 1:ne
    % Prepare the NUFFT
    FT = NUFFT(Traj(:,:,ec), 1, 1, [0,0], fftSize, 2);
    %         (k,             w,phase,shift,imSize, mode)
    
    for sl=1:nz %slices
%         fprintf('Starting Gridding recon of echo %d, slice %d / %d\n', ec, (sl-slices(1)+1), length(slices));
        fprintf('Starting Gridding recon of echo %d, slice %d / %d\n', ec, sl, nz);
        % Loop over channels and calculate image for each coil
        for channel=1:nc
            % Fetch rawdata for slice and channel and multiply with density
            % compensation function
%             workKData=dcf.*double(squeeze(kdata(:,channel,:,sl-(slices(1)-1))));
            workKData=dcf.*double(squeeze(kdata(:,channel,:,sl,ec)));
            
            % Run the NUFFT
            workImage=FT'*(workKData);
            
            % Add squared channel image to buffer
%             out_img(:,:,sl-(slices(1)-1)) = out_img(:,:,sl-(slices(1)-1)) + abs(workImage.*workImage);
            out_img(:,:,sl,ec) = out_img(:,:,sl,ec) + abs(workImage.*workImage);
        end
    end
end
% Calculate the root (as final part of the sum-of-squares calculation)
out_img = sqrt(out_img);
out_img = out_img./max(out_img(:));