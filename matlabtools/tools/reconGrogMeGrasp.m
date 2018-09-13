function [out_img, t] = reconGrogMeGrasp(kdata, doGpu)
% Function for performing a reconstruction of multi-echo golden-angle 
% radial sparse (GRASP) k-space data using compressed sensing with GROG 
% pre-interpolation. No binning into different motion or temporal states is 
% performed.
% NOTE: it is assumed that a fourier transform along the slice direction
% has already been performed!
% Input parameters:
% kdata: contains the raw k-space data. Must have dimensions
% [nx,nc,ntviews,nz,ne], which is the same order they are read in by 
% mapVBVD.
% doGpu: optional parameter specifying whether GPU should be used
% 
% Adapted from Demo4_Reconstruction.m from RACER-GRASP_GROG-GRASP demo
% package (NYU Demo provided by Li Feng)
% by Marnix Maas (Marnix.Maas@radboudumc.nl), August 2018

if nargin<2, doGpu = 0; end

% Permute dimensions of kdata to [nx,ntviews,nz,nc,ne]
kdata = permute(kdata, [1,3,4,2,5]);
[nx,ntviews,nz,nc,ne]=size(kdata);

% Set some initial parameters
% slices = 1:nz;
% nLinDyn = ntviews;
% nt=floor(ntviews/nLinDyn);              % Number of dynamic time points
Nqu=min(floor(nx/2 * pi/2), ntviews);   % Nyquist spokes, please refer to the GROG-GRASP paper for more details about this.
% MCM changed from 350
bas = nx;

% Generate trajectory (For GROG only, NOT for general gridding)
%     Traj=Trajectory_GoldenAngle(ntviews,nx,0,[0 1]);        % GROG expects no normalization and reversed spoke traversal
Traj=Trajectory_GoldenAngle_ME(ntviews, nx, ne, 2, 'normalize',0,'flip',[0 1]);
Traj=Traj*(bas/nx);

% Initialize output image & recon time measurement
out_img      = single(zeros(nx,nx,nz,ne));
t = [];
for ec=1:ne
    for sl=1:nz % loop through selected slices
        tic;
        fprintf('Starting GROG-GRASP recon of echo %d, slice %d / %d\n', ec, sl, nz);
        % Select slice to be reconstructed
        kdata1 = kdata(:,:,sl,:,ec);
%         [nx,ntviews,~,nc,ne]=size(kdata1);         % MCM: seems redundant, but we have a new number of coils since compressing!
        % Should be inside slice loop, variables are overwritten again below
        
        % Initialize GROG 
        [Gx,Gy] = GROG.get_Gx_Gy(kdata1,Traj(:,:,ec));  % weights calibration
        G = GROG.init(kdata1,Traj(:,:,ec),Gx,Gy,0);     % GROG initialization
        
        % Do GROG; the preinterpolation of radial data into a Cartesian grid
        kdata1 = GROG.interp(kdata1,G,1,doGpu);         % kdata1 is cartesian from here on
        kdata1 = CropImg(kdata1,bas,bas);
        
        % Coil sensitivities estimation
        ref = gather(squeeze(ifft2c_mri(kdata1)));
        b1  = adapt_array_2d(ref);clear ref
        b1  = single(b1/max(abs(b1(:))));
        
        % Calculate density compensation 
        D   = reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1))]);
        DCF = CropImg(D,bas,bas);
        DCF = repmat(DCF,[1,1,nc]);        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Iterative reconstruction
        
%         kdata1=reshape(kdata1,[bas,bas,nresp,nt,nc]);
        kdata1=reshape(kdata1,[bas,bas,nc]);        % Verdandert ws niets
        
        % Permute k-space data to correct order (depending on
        % algorithm to be used)
%         kdata1=permute(kdata1,[1,2,3,5,4]);         % Putting resipratory motion 3rd, dynamics 5th
        
        % Prepare iterative recon
        Weighting = DCF;
        param.y=kdata1.*sqrt(Weighting);
%         param.y=kdata1.*sqrt(DCF);
        mask=single(kdata1~=0);
        param.E=Emat_GROG2DSG(mask,b1,Weighting);
        recon_cs=param.E'*param.y;                  % Starting point for iterative recon: Fourier transform
%         Hierna heeft recon_cs 3 dimensies: 3e dimensie = coils. Dat lijkt
%         me niet de bedoeling?
% Misschien kan ik hier een 'gewone' 3D MCNUFFT operator pakken. Of
% uitvogelen hoe dit had moeten zijn.
        
        % Set parameters for optimizer
        param.TV_dim1=TV_Temp;
        param.TVWeight_dim1=max(abs(recon_cs(:)))*0.02;
        param.TVWeight_dim2=0;
        param.nite = 4;
        param.display=1;
        
        % Run iterative recon
        for n=1:2
            recon_cs = CSL1NlCg_XDGRASP_Mx(recon_cs,param);
        end
        %                 recon_cs=CropImg(recon_cs,nx/2,nx/2);
        out_img(:,:,sl,:,ec)=gather(single(recon_cs));
        
        t(ec,sl) = toc;
    end %slice loop
end %echo loop
out_img=out_img/max(abs(out_img(:)));



    