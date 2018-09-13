function out_img = reconGriddingXdGrasp(kdata, Res_Signal, nresp)
% Function for performing a respiratory motion resolved reconstruction of
% golden-angle radial sparse (GRASP) k-space data using non-cartesian 
% Fourier transformations ('gridding' reconstruction). 
% NOTE: it is assumed that a fourier transform along the slice direction
% has already been performed!
% Input parameters:
% kdata: contains the raw k-space data. Must have dimensions
% [nx,nc,ntviews,nz], which is the same order they are read in by mapVBVD.
% Res-Signal: respiratory motion signal (estimated in a previous step)
% nresp: number of respiratory motion states
% Adapted from the demo file 'iterative_multicoil.m' by Tobias Block
% by Marnix Maas (Marnix.Maas@radboudumc.nl), August 2018
    
kdata = permute(kdata, [1,3,4,2]);
[nx,ntviews,nz,nc]=size(kdata);

% Generate trajectory
Traj=Trajectory_GoldenAngle(ntviews,nx,1,[0 0]);

% Calculate density compensation
dcf = dcfGridding(ntviews, nx);

% Set NUFFT parameters
fftSize=[nx,nx];

% Calculate number of spokes per respiratory phase
nLinRes=floor(ntviews/nresp);

% Sort the trajectory and dcf according to respiratory motion position
[~,index]=sort(Res_Signal,'descend');
Traj=Traj(:,index);
dcf=dcf(:,index);

% Divide k-space trajectory and dcf into motion states
Traj_u  = zeros(nx,nLinRes,nresp);
dcf_u   = Traj_u;
% FT_u    = zeros(1,nresp);
for ii=1:nresp
    Traj_u(:,:,ii)=Traj(:,(ii-1)*nLinRes+1:ii*nLinRes);
    dcf_u(:,:,ii)=dcf(:,(ii-1)*nLinRes+1:ii*nLinRes);
    FT_u(ii) = NUFFT(squeeze(Traj_u(:,:,ii)), 1, 1, [0,0], fftSize, 2);
end

% Create empty images to sum the reconstructed channels (for sum-of-squares combination)
slices = 1:nz;
out_img   = single(zeros(nx,nx,length(slices),nresp));

%%%% Motion resolved gridding reconstruction
for sl=slices
    fprintf('Starting XD-Gridding recon of slice %d / %d\n', (sl-slices(1)+1), length(slices));
    % Select slice to be reconstructed
    kdata1=squeeze(kdata(:,:,sl-(slices(1)-1),:));
    % Sort k-space data according to respiratory motion
    kdata1=kdata1(:,index,:);
    % Divide data into motion states
    kdata_u = zeros(nx,nLinRes,nc,nresp);
    for ii=1:nresp
        kdata_u(:,:,:,ii)=kdata1(:,(ii-1)*nLinRes+1:ii*nLinRes,:);
        %             kdata_u(:,:,:,ii)=kdata1(:,(1:nLinRes)+ii-1,:);
        %         (k,w,phase,shift,imSize, mode)
        
        for channel=1:nc
            % Fetch k-data for slice and channel and multiply with density
            % compensation function
            workKData=dcf_u(:,:,ii).*double(squeeze(kdata_u(:,:,channel,ii)));
            % Run the NUFFT
            workImage=FT_u(ii)'*(workKData);
            % Add squared channel image to buffer
            out_img(:,:,sl-(slices(1)-1),ii) = out_img(:,:,sl-(slices(1)-1),ii) + abs(workImage.*workImage);
        end
    end
end

% Calculate the root (as final part of the sum-of-squares calculation) and
% normalize.
out_img = sqrt(out_img);
out_img = out_img./max(out_img(:));