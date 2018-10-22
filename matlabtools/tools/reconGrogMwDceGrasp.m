function [out_img, t] = reconGrogMwDceGrasp(kdata, Res_Signal, pars)
% Function for performing a respiratory motion weighted reconstruction of
% golden-angle radial sparse (GRASP) k-space data using compressed sensing 
% with GROG pre-interpolation.
% NOTE: it is assumed that a fourier transform along the slice direction
% has already been performed!
% Input parameters:
% - kdata: contains the raw k-space data. Must have dimensions
%          [nx,nc,ntviews,nz], which is the same order they are read in by 
%          mapVBVD.
% - Res-Signal: respiratory motion signal (estimated in a previous step)
% - pars: structure containing reconstruction parameters, eg:
% - pars.nresp: number of respiratory motion states
% - pars.doGpu: flag specifying whether GPU should be used
% 
% Adapted from Demo4_Reconstruction.m from RACER-GRASP_GROG-GRASP demo
% package (NYU Demo provided by Li Feng)
% by Marnix Maas (Marnix.Maas@radboudumc.nl), August 2018

% Permute dimensions of kdata to [nx,ntviews,nz,nc]
kdata               = permute(kdata, [1,3,4,2]);
[nx,ntviews,nz,nc]  = size(kdata);

% Set some initial parameters
if pars.nLinDyn == 0
    nLinDyn = pars.nresp*floor(ntviews/pars.nresp);     % Ensure that final number of used spokes is divisable by number of respiratory states
else
    nLinDyn = pars.nLinDyn;
end
nt = floor(ntviews/nLinDyn);                            % Number of dynamic time points
Nqu = min(floor(nx/2 * pi/2), ntviews);                 % Nyquist spokes, please refer to the GROG-GRASP paper for more details about this.
% MCM changed from 350

% Initialize output image & recon time measurement
% out_img      = single(zeros(nx,nx,nz,pars.nresp,nt));
out_img      = single(zeros(pars.bas,pars.bas,nz,nt));
t = [];
for sl=1:nz % loop through selected slices
    tic;
    fprintf('Starting GROG MW-GRASP recon of slice %d / %d\n', sl, nz);
    % Select slice to be reconstructed
    kdata1=kdata(:,:,sl,:);                     % Kunnen we dit squeezen?
%     [nx,ntviews,~,nc]=size(kdata1);         % MCM: seems redundant, but
%     ntviews was altered below in a previous version. Consider putting
%     this line back to avoid future problems?
%     % Should be inside slice loop, variables are overwritten again below
    
    % Generate trajectory (For GROG only, NOT for general gridding)
    % This needs to stay within the slice loop; it gets modified during
    % sorting of data in respiratory phases
    Traj=Trajectory_GoldenAngle(ntviews,nx,0,[0 1]);        % GROG expects no normalization and reversed spoke traversal
    Traj=Traj*(pars.bas/nx);
    
    %GROG weights calibration
    [Gx,Gy] = GROG.get_Gx_Gy(kdata1,Traj);
    
    %Coil sensitivities estimation
    G       = GROG.init(kdata1,Traj,Gx,Gy,0);       % GROG initialization
    kref    = GROG.interp(kdata1,G,1,pars.doGpu);   % GROG interpolation; MCM changed GPU flag (4th param) from 0
    kref    = CropImg(kref,pars.bas,pars.bas);
    ref     = gather(squeeze(ifft2c_mri(kref)));
    b1      = adapt_array_2d(ref); clear ref
    b1      = single(b1/max(abs(b1(:))));
    
    % Calculate density compensation for 'Nyquist' number of spokes
    % (see paper, this is the denominator in the Weights calculation)
    % MCM: why does this use the _last_ Nqu spokes of the full
    % trajectory? What's the difference?
    G       = GROG.init(kdata1(:,end-Nqu+1:end,:,:),Traj(:,end-Nqu+1:end),Gx,Gy,0);
    D       = reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1))]);
    DCF     = CropImg(D,pars.bas,pars.bas);
    DCF     = repmat(DCF,[1,1,nt,nc]);              % size [bas,bas,nt,nc]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data sorting
    
    % Split trajectory, k-space data and respiratory signal in nt
    % dynamic time points & throw away the last few lines that don't
    % add up to a full dynamic
    Traj                = reshape(Traj(:,1:nt*nLinDyn),[nx,nLinDyn,nt]);            % Hier wordt kdata gecropt naar nt*nLinDyn spokes 
    kdata1              = reshape(kdata1(:,1:nt*nLinDyn,:,:),[nx,nLinDyn,nt,nc]);   % van [nx,ntviews,[1],nc] naar [nx,nLinDyn,nt,[1],nc]
    Res_Signal_tmp      = reshape(Res_Signal(1:nt*nLinDyn),[nLinDyn,nt]);
    [~,~,nt,~]          = size(kdata1);
    
    % Sort trajectory and k-space data according to the respiratory
    % motion signal in each dynamic time point
    for ii=1:nt
        [~,index]           = sort(Res_Signal_tmp(:,ii),'descend');
        Traj(:,:,ii)        = Traj(:,index,ii);
        kdata1(:,:,ii,:)    = kdata1(:,index,ii,:);
    end
    
    % Calculate density compensation for each dynamic time point
    % (see paper, this is the numerator in the Weights calculation)
%     [nx,ntviews,nt,nc]  = size(kdata1);          %MCM Necessary: does this change? --> Hier wordt alleen ntviews gelijk gezet aan nLinDyn, de rest verandert bmbw niet. Regel kan dus weg.
    G                   = GROG.init(kdata1,Traj,Gx,Gy,0);
    DCF_U               = reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1)),nt]);
    DCF_U               = CropImg(DCF_U,pars.bas,pars.bas);         % size [bas,bas,nt?]
    DCF_U               = repmat(DCF_U,[1,1,1,nc]);                 % size [bas,bas,nt,nc?]
    DCF_U(DCF_U==0)     = 1;
    
    % calculate GROG weights, as described in the paper
    Weighting           = DCF./DCF_U;
    Weighting           = repmat(Weighting,[1 1 1 1 pars.nresp]);   % size [bas,bas,nt,nc,nresp?]
    
    % Prepare k-space data for respiratory and dynamic time point
    % sorted GROG (cartesian interpolation)
    % Je wilt hier de data zo omvormen dat hij bestaat uit brokjes van het
    % aantal lijntjes per resp-bin binnen elk dynamisch tijdspunt. Dit om
    % elk blokje te kunnen GROG-en. De respiratory en dynamic tijdspunten
    % worden hierbij in een dimensie geveegd; kennelijk kan de GROG-routine
    % daarmee omgaan.
    nLinDynResp = floor(nLinDyn/pars.nresp);                            % Number of lines per respiratory bin within each dynamic time point
    kdata1      = reshape(kdata1,[nx,nLinDynResp,pars.nresp,nt,nc]);    % van [nx,nLinDyn,nt,nc] naar [nx,nLinDynResp,nresp,nt,nc]
    kdata1      = reshape(kdata1,[nx,nLinDynResp,pars.nresp*nt,nc]);    % naar [nx,nLinDynResp,nresp*nt,nc]
    Traj        = reshape(Traj,[nx,nLinDynResp,pars.nresp,nt]);
    Traj        = reshape(Traj,[nx,nLinDynResp,pars.nresp*nt]);
    
    % Do GROG; the preinterpolation of radial data into a Cartesian grid
    G           = GROG.init(kdata1,Traj,Gx,Gy,0);
    kdata1      = GROG.interp(kdata1,G,1,pars.doGpu);                   % Vanaf hier is kdata1 cartesisch!
    kdata1      = CropImg(kdata1,pars.bas,pars.bas);                    % ...met de 4 dimensies [bas,bas,nresp*nt,nc]
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iterative reconstruction
    
    kdata1=reshape(kdata1,[pars.bas,pars.bas,pars.nresp,nt,nc]);        % terug naar [bas,bas,nresp,nt,nc]: respiratory en tijdspunten weer in aparte dimensies
    
    % Permute k-space data to correct order (depending on algorithm to be 
    % used)
%     kdata1      = permute(kdata1,[1,2,3,5,4]);                          % [bas,bas,nresp,nc,nt] Putting resipratory motion 3rd, dynamics 5th
    kdata1      = permute(kdata1,[1,2,4,5,3]);     % Als origineel in demo. Werkt het nog als ik dit voor single time-point wil gebruiken?
    
    % Prepare iterative recon
    SoftWeight=single(zeros(pars.bas,pars.bas,nt,nc,pars.nresp));
    for jj=1:pars.nresp
        SoftWeight(:,:,:,:,jj)=(1/4)^(jj-1);
    end
    
    mask=single(kdata1~=0);
    if pars.doGpu
        param.SG = gpuArray(SoftWeight);
    else
        param.SG = SoftWeight;
    end
    
    param.y=kdata1.*sqrt(Weighting);
    param.E=Emat_GROG2DSG(mask,b1,Weighting);
    recon_cs=param.E'*param.y;
    
    param.TV=TV_Temp;
    Weight1=0.04;
    param.TVWeight=max(abs(recon_cs(:)))*Weight1;
    param.nite = pars.nIterIn;
    param.display = 1;
    
    for n=1:pars.nIterOut
        % Motion _weighted_ reconstruction
        recon_cs = CSL1NlCg_4DDCE(recon_cs,param);
    end
%     recon_cs=CropImg(recon_cs,nx/2,nx/2);
    
%     out_img(:,:,sl,:,:)       = gather(single(recon_cs));
    out_img(:,:,sl,:)       = gather(single(recon_cs));
    
%     param.y     = kdata1.*sqrt(Weighting);
%     mask        = single(kdata1~=0);
%     param.E     = Emat_GROG2DSG(mask,b1,Weighting);
%     recon_cs    = param.E'*param.y;                                     % Starting point for iterative recon: Fourier transform
%     
%     % Set parameters for optimizer
%     param.TV_dim1           = TV_Temp;
%     param.TVWeight_dim1     = max(abs(recon_cs(:)))*0.02;
%     param.TVWeight_dim2     = 0;
%     param.nite              = 4;
%     param.display           = 1;
%     
%     % Run iterative recon
% %     for n=1:2
% %         recon_cs            = CSL1NlCg_XDGRASP_Mx(recon_cs,param);
% %     end
%     for n=1:3
%         recon_cs            = CSL1NlCg_4DDCE(recon_cs,param);
%     end
%     %                 recon_cs=CropImg(recon_cs,nx/2,nx/2);
    out_img(:,:,sl,:)       = gather(single(recon_cs));
    
    t(sl) = toc;
end %slice loop

% Normalize output image
out_img=out_img/max(abs(out_img(:)));



    