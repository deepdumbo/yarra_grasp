function out_img = reconGriddingXdMeGrasp(kdata, Res_Signal, nresp)
% Function for performing a respiratory motion resolved reconstruction of
% multi-echo golden-angle radial sparse (GRASP) k-space data using 
% non-cartesian Fourier transformations ('gridding' reconstruction). 
% NOTE: it is assumed that a fourier transform along the slice direction
% has already been performed!
% Input parameters:
% kdata: contains the raw k-space data that has been Fourier transformed 
% along z. Must have dimensions [nx,nc,ntviews,nz], which is the same 
% order they are read in by mapVBVD.
% Res-Signal: respiratory motion signal (estimated in a previous step)
% nresp: number of respiratory motion states
% Adapted from the demo file 'iterative_multicoil.m' by Tobias Block
% by Camille Van Speybroeck (Camille.vanSpeybroeck@radboudumc.nl), 
% November 2018
    
%kdata = permute(kdata, [1,3,4,2,5]);           k-data hoeft niet gepermute te worden.
[nx,nc,ntviews,nz,ne]=size(kdata);

% Generate trajectory.
Traj=Trajectory_GoldenAngle_ME(ntviews, nx, ne, 2, 'normalize', 1, 'flip', [0 0]);

% Calculate density compensation.
dcf = dcfGridding(ntviews, nx);

% Set NUFFT parameters.
fftSize=[nx,nx];

% Calculate number of spokes per respiratory phase.
nLinRes=floor(ntviews/nresp);

% Sort the trajectory and dcf according to respiratory motion position.
[~,index]=sort(Res_Signal,'descend');
Traj=Traj(:,index,:);
dcf=dcf(:,index);

% Divide k-space trajectory and dcf into motion states.
Traj_u  = zeros(nx,nLinRes,nresp,ne);
dcf_u   = Traj_u;
% FT_u    = zeros(1,nresp);
for ii=1:nresp                     
    ec=ne;
    Traj_u(:,:,ii,:)=Traj(:,(ii-1)*nLinRes+1:ii*nLinRes,:);
    FT_u(ii,:) = NUFFT(squeeze(Traj_u(:,:,ii,:)), 1, 1, [0,0], fftSize, 2);
    for ec=1:ne                  
        dcf_u(:,:,ii,ec)=dcf(:,(ii-1)*nLinRes+1:ii*nLinRes);
        FT_u(ii,ec) = NUFFT(squeeze(Traj_u(:,:,ii,ec)), 1, 1, [0,0], fftSize, 2);
    end
end

% Create empty images to sum the reconstructed channels (for sum-of-squares
% combination).
slices = 1:nz;
out_img   = single(zeros(nx,nx,length(slices),nresp,ne));

%%%% Multi echo Motion resolved gridding reconstruction
for sl=slices
    for ec=1:ne
        fprintf('Starting ME XD-Gridding recon of echo %d, slice %d / %d\n', ec, sl, length(slices));
        % Select slice to be reconstructed
        kdata1=squeeze(kdata(:,:,:,sl,ec)); 
        % Sort k-space data according to respiratory motion
        kdata1=kdata1(:,:,index);
        % Divide data into motion states
        kdata_u = zeros(nx,nc,nLinRes,nresp);
        for ii=1:nresp
            kdata_u(:,:,:,ii)=kdata1(:,:,(ii-1)*nLinRes+1:ii*nLinRes);
            %             kdata_u(:,:,:,ii)=kdata1(:,(1:nLinRes)+ii-1,:);
            %         (k,w,phase,shift,imSize, mode)

            for channel=1:nc
                % Fetch k-data for slice and channel and multiply with density
                % compensation function.
                workKData=dcf_u(:,:,ii,ec).*double(squeeze(kdata_u(:,channel,:,ii)));
                % Run the NUFFT.
                workImage=FT_u(ii,ec)'*(workKData);
                % Add squared channel image to buffer.
                out_img(:,:,sl,ii,ec) = out_img(:,:,sl,ii,ec) + abs(workImage.*workImage);
            end
        end
    end
end

% Calculate the root (as final part of the sum-of-squares calculation) and
% normalize.
out_img = sqrt(out_img);
out_img = out_img./max(out_img(:));