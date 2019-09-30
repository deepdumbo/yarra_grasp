%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This demo demonstrates the respiratory weighted iterative reconstruction

% GROG is used as a pre-interpolatation step to avoid gridding in the 
% iterative process

% The bolus signal can be used to guide the reconstruction.
% For simplicity, this is not included in this demo

% Li Feng, NYU, 12/18/2017
% Li.Feng@nyumc.org

clear
clc
addpath Data
addpath nufft_files
addpath utils
close all
load Data/kdata_Unstreaking_CoilCompression.mat
load Data/Res_Signal.mat

TA=178; % total acquisition time, this information was obtained from the raw data header file.
kdata1=kdata;clear kdata
[nx,ntviews,nz,nc]=size(kdata1);
bas=360;% the 2x full FOV is 512 and we only reconstruct 360 here.

% for zz=1:nz % loop through all slices
    zz=1 % For demonstrate purpose, reconstruct 1 slice only
    kdata=kdata1(:,:,zz,:);
    [nx,ntviews,~,nc]=size(kdata);
    
    %Generating trajectory (For GROG only, NOT for general gridding)
    Traj=Trajectory_GoldenAngle_GROG(ntviews,nx);%trajectory
    %ivom: waarom die factors?
    Traj=Traj*(bas/nx);
    
    %GROG weights calibration
    [Gx,Gy] = GROG.get_Gx_Gy(kdata,Traj);
    
    %Coil sensitivities estimation
    G = GROG.init(kdata,Traj,Gx,Gy,0); % GROG initialization 
    kref = GROG.interp(kdata,G,1,0); % GROG interpolation 
    kref=CropImg(kref,bas,bas);
    ref=gather(squeeze(ifft2c_mri(kref)));
    b1=adapt_array_2d(ref);clear ref
    b1=single(b1/max(abs(b1(:))));
    
    nline=96;nt=floor(ntviews/nline);
    Nqu=350; % Nyquist spokes, please refer to the GROG-GRASP paper for more details about this.
    
    G = GROG.init(kdata(:,end-Nqu+1:end,:,:),Traj(:,end-Nqu+1:end),Gx,Gy,0);
    D=reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1))]);
    DCF=CropImg(D,bas,bas);
    DCF=repmat(DCF,[1,1,nt,nc]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data sorting
    Traj=reshape(Traj(:,1:nt*nline),[nx,nline,nt]);
    kdata=reshape(kdata(:,1:nt*nline,:,:),[nx,nline,nt,nc]);
    Res_Signal_tmp=reshape(Res_Signal(1:nt*nline),[nline,nt]);
    [nx,ntviews,nt,nc] = size(kdata);
    
    % Sorting according to the respiratory motion signal
    for ii=1:nt
        [~,index]=sort(Res_Signal_tmp(:,ii),'descend');
        Traj(:,:,ii)=Traj(:,index,ii);
        kdata(:,:,ii,:)=kdata(:,index,ii,:);
    end
    
    [nx,ntviews,nt,nc] = size(kdata);
    G = GROG.init(kdata,Traj,Gx,Gy,0);
    DCF_U=reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1)),nt]);
    DCF_U=CropImg(DCF_U,bas,bas);
    
    ntres=4;% 4 respiratory phases
    nline2=floor(ntviews/ntres);
    kdata=reshape(kdata,[nx,nline2,ntres,nt,nc]);
    kdata=reshape(kdata,[nx,nline2,ntres*nt,nc]);
    Traj=reshape(Traj,[nx,nline2,ntres,nt]);
    Traj=reshape(Traj,[nx,nline2,ntres*nt]);
    
    % Do GROG; the preinterpolation of radial data into a Cartesian grid
    G = GROG.init(kdata,Traj,Gx,Gy,0);
    kdata = GROG.interp(kdata,G,1,0);
    kdata=CropImg(kdata,bas,bas);
    
    % calculate GROG weights, as described in the paper
    DCF_U=repmat(DCF_U,[1,1,1,nc]);
    DCF_U(find(DCF_U==0))=1;
    Weighting=DCF./DCF_U;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    kdata=reshape(kdata,[bas,bas,ntres,nt,nc]);
    kdata=permute(kdata,[1,2,4,5,3]);
    
    SoftWeight=single(zeros(bas,bas,nt,nc,ntres));
    for jj=1:ntres
        SoftWeight(:,:,:,:,jj)=(1/4)^(jj-1);
    end
    
    mask=single(kdata~=0);
    %ivom: command to send process to the gpu?
    %ivom: placing something in parentheses as gpuArray( A ) sends the
    %ivom: bracketed array to the GPU. All operations following the 
    %ivom: gpuArray( A ) command on array A are then processed on the GPU.
    %param.SG=gpuArray(SoftWeight);
    
    param.SG = SoftWeight;
    Weighting=repmat(Weighting,[1 1 1 1 ntres]);
    
    param.y=kdata.*sqrt(Weighting);
    param.E=Emat_GROG2DSG(mask,b1,Weighting);
    recon_cs=param.E'*param.y;
    
    %ivom: temporal finite difference operator TV_Temp.
    param.TV=TV_Temp;
    Weight1=0.04;
    param.TVWeight=max(abs(recon_cs(:)))*Weight1;
    param.nite = 8;param.display = 1;
    
    clc
    for n=1:3
        recon_cs = CSL1NlCg_4DDCE(recon_cs,param);
    end
    recon_cs=CropImg(recon_cs,nx/2,nx/2);
    data(:,:,:,zz)=gather(abs(single(recon_cs)));
% end
    data=data/max(data(:));
    save Data/Result.mat data 

    