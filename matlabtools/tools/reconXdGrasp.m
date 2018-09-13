function [out_img_cs, out_img_mcgridding, t] = reconXdGrasp(kdata, Res_Signal, nresp)
% Function for performing a respiratory motion resolved reconstruction of
% golden-angle radial sparse (GRASP) k-space data using compressed sensing. 
% NOTE: it is assumed that a fourier transform along the slice direction
% has already been performed!
% Input parameters:
% kdata: contains the raw k-space data. Must have dimensions
% [nx,nc,ntviews,nz], which is the same order they are read in by mapVBVD.
% Res-Signal: respiratory motion signal (estimated in a previous step)
% nresp: number of respiratory motion states
% 
% Adapted from Demo_XDGRASP_NonContrast (NYU Demo provided by Li Feng)
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
out_img_cs   = single(zeros(nx,nx,length(slices),nresp));
out_img_mcgridding = out_img_cs;


%%% Motion resolved gridding and CS-recon using MCNUFFT 
b1 = ones(nx, nx, nc);   % temporary solution: implement b1 map estimation!
t = [];
tic;
for sl=slices
    fprintf('Starting XD-GRASP recon of slice %d / %d\n', (sl-slices(1)+1), length(slices));
    % Select slice to be reconstructed
    kdata1=kdata(:,:,sl-(slices(1)-1),:);
    
    % Estimate coil sensitivities using GROG (should this be done for every
    % respiratory state separately?)
    [Gx,Gy] = GROG.get_Gx_Gy(kdata1,Traj);
    G = GROG.init(kdata1,Traj,Gx,Gy,0);     % GROG initialization
    kref = GROG.interp(kdata1,G,1,0);       % GROG interpolation; 4th param is GPU flag, now set to 0 by default
%     kref=CropImg(kref,bas,bas);
    ref=gather(squeeze(ifft2c_mri(kref)));
    b1=adapt_array_2d(ref);clear ref
    b1=single(b1/max(abs(b1(:))));
    
    % Sort k-space data according to respiratory motion
    kdata1 = squeeze(kdata1);
    kdata1=kdata1(:,index,:);
    % Perform density compensation
    kdata1=kdata1.*repmat(sqrt(dcf),[1,1,nc]);
    % Divide data into motion states
    kdata_u = zeros(nx,nLinRes,nc,nresp);
    for ii=1:nresp
        kdata_u(:,:,:,ii)=kdata1(:,(ii-1)*nLinRes+1:ii*nLinRes,:);
    end
    
    % Perform multi-coil XD-Gridding recon as starting point for CS
    param.E=MCNUFFT(double(Traj_u),double(dcf_u),double(b1));
    param.y=double(kdata_u);
    recon_cs = param.E'*param.y;
    out_img_mcgridding(:,:,sl-(slices(1)-1),:)=recon_cs;
    
%     % Perform compressed sensing recon
%     param.TV_dim1=TV_Temp;
%     param.TVWeight_dim1=max(abs(recon_cs(:)))*0.02;
%     param.TVWeight_dim2=0;
%     param.nite = 4;
%     param.display=1;
%     
%     %         tic
%     for n=1:2
%         recon_cs = CSL1NlCg_XDGRASP(recon_cs,param);
%     end
%     %         time=toc;
    %         time=time/60
    out_img_cs(:,:,sl-(slices(1)-1),:)=recon_cs;
    t(sl-(slices(1)-1)) = toc;
end

% Normalize images
out_img_mcgridding = out_img_mcgridding/max(abs(out_img_mcgridding(:)));
out_img_cs         = out_img_cs/max(abs(out_img_cs(:)));

disp('...done.');