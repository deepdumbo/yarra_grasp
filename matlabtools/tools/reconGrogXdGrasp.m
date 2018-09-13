function [out_img, t] = reconGrogXdGrasp(kdata, Res_Signal, nresp, doGpu)
% Function for performing a respiratory motion resolved reconstruction of
% golden-angle radial sparse (GRASP) k-space data using compressed sensing 
% with GROG pre-interpolation.
% NOTE: it is assumed that a fourier transform along the slice direction
% has already been performed!
% Input parameters:
% kdata: contains the raw k-space data. Must have dimensions
% [nx,nc,ntviews,nz], which is the same order they are read in by mapVBVD.
% Res-Signal: respiratory motion signal (estimated in a previous step)
% nresp: number of respiratory motion states
% doGpu: optional parameter specifying whether GPU should be used
% 
% Adapted from Demo4_Reconstruction.m from RACER-GRASP_GROG-GRASP demo
% package (NYU Demo provided by Li Feng)
% by Marnix Maas (Marnix.Maas@radboudumc.nl), August 2018

if nargin<4, doGpu = 0; end

% Permute dimensions of kdata to [nx,ntviews,nz,nc]
kdata = permute(kdata, [1,3,4,2]);
[nx,ntviews,nz,nc]=size(kdata);

% Set some initial parameters
slices = 1:nz;
% nLinDyn = ntviews;
nLinDyn = nresp*floor(ntviews/nresp);
nt=floor(ntviews/nLinDyn);              % Number of dynamic time points
Nqu=min(floor(nx/2 * pi/2), ntviews);   % Nyquist spokes, please refer to the GROG-GRASP paper for more details about this.
% MCM changed from 350
bas = nx;

% Initialize output image & recon time measurement
out_img      = single(zeros(nx,nx,length(slices),nresp));
t = [];
for sl=slices % loop through selected slices
    tic;
    fprintf('Starting GROG XD-GRASP recon of slice %d / %d\n', (sl-slices(1)+1), length(slices));
    % Select slice to be reconstructed
    kdata1=kdata(:,:,sl-(slices(1)-1),:);
    [nx,ntviews,~,nc]=size(kdata1);         % MCM: seems redundant, but we have a new number of coils since compressing!
    % Should be inside slice loop, variables are overwritten again below
    
    % Generate trajectory (For GROG only, NOT for general gridding)
    % This needs to stay within the slice loop; it gets modified during
    % sorting of data in respiratory phases
    Traj=Trajectory_GoldenAngle(ntviews,nx,0,[0 1]);        % GROG expects no normalization and reversed spoke traversal
    Traj=Traj*(bas/nx);
    
    %GROG weights calibration
    [Gx,Gy] = GROG.get_Gx_Gy(kdata1,Traj);
    
    %Coil sensitivities estimation
    G = GROG.init(kdata1,Traj,Gx,Gy,0);     % GROG initialization
    kref = GROG.interp(kdata1,G,1,doGpu);   % GROG interpolation; MCM added changed GPU flag (4th param) from 0
    kref=CropImg(kref,bas,bas);
    ref=gather(squeeze(ifft2c_mri(kref)));
    b1=adapt_array_2d(ref);clear ref
    b1=single(b1/max(abs(b1(:))));
    
    % Calculate density compensation for 'Nyquist' number of spokes
    % (see paper, this is the denominator in the Weights calculation)
    % MCM: why does this use the _last_ Nqu spokes of the full
    % trajectory? What's the difference?
    G = GROG.init(kdata1(:,end-Nqu+1:end,:,:),Traj(:,end-Nqu+1:end),Gx,Gy,0);
    D=reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1))]);
    DCF=CropImg(D,bas,bas);
    DCF=repmat(DCF,[1,1,nt,nc]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data sorting
    
    % Split trajectory, k-space data and respiratory signal in nt
    % dynamic time points & throw away the last few lines that don't
    % add up to a full dynamic
    Traj=reshape(Traj(:,1:nt*nLinDyn),[nx,nLinDyn,nt]);
    kdata1=reshape(kdata1(:,1:nt*nLinDyn,:,:),[nx,nLinDyn,nt,nc]);
    Res_Signal_tmp=reshape(Res_Signal(1:nt*nLinDyn),[nLinDyn,nt]);
    [nx,ntviews,nt,nc] = size(kdata1);
    
    % Sort trajectory and k-space data according to the respiratory
    % motion signal in each dynamic time poin[t
    for ii=1:nt
        [~,index]=sort(Res_Signal_tmp(:,ii),'descend');
        Traj(:,:,ii)=Traj(:,index,ii);
        kdata1(:,:,ii,:)=kdata1(:,index,ii,:);
    end
    
    % Calculate density compensation for each dynamic time point
    % (see paper, this is the numerator in the Weights calculation)
    [nx,ntviews,nt,nc] = size(kdata1);
    G = GROG.init(kdata1,Traj,Gx,Gy,0);
    DCF_U=reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1)),nt]);
    DCF_U=CropImg(DCF_U,bas,bas);
    DCF_U=repmat(DCF_U,[1,1,1,nc]);
    DCF_U(DCF_U==0)=1;
    
    % calculate GROG weights, as described in the paper
    Weighting=DCF./DCF_U;
    Weighting=repmat(Weighting,[1 1 1 1 nresp]);
    
    % Prepare k-space data for respiratory and dynamic time point
    % sorted GROG (cartesian interpolation)
%     nline2=floor(ntviews/nresp);                            % Number of lines per respiratory bin
    nline2=floor(nLinDyn/nresp);                            % Number of lines per respiratory bin
    kdata1=reshape(kdata1,[nx,nline2,nresp,nt,nc]);
    kdata1=reshape(kdata1,[nx,nline2,nresp*nt,nc]);
    Traj=reshape(Traj,[nx,nline2,nresp,nt]);
    Traj=reshape(Traj,[nx,nline2,nresp*nt]);
    
    % Do GROG; the preinterpolation of radial data into a Cartesian grid
    G = GROG.init(kdata1,Traj,Gx,Gy,0);
    kdata1 = GROG.interp(kdata1,G,1,doGpu);                 % MCM added changed GPU flag (4th param) from 0
    kdata1 = CropImg(kdata1,bas,bas);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iterative reconstruction
    
    kdata1=reshape(kdata1,[bas,bas,nresp,nt,nc]);
    
    % Permute k-space data to correct order (depending on
    % algorithm to be used)
    kdata1=permute(kdata1,[1,2,3,5,4]);         % Putting resipratory motion 3rd, dynamics 5th
    
    % Prepare iterative recon
    param.y=kdata1.*sqrt(Weighting);
    mask=single(kdata1~=0);
    param.E=Emat_GROG2DSG(mask,b1,Weighting);
    recon_cs=param.E'*param.y;                  % Starting point for iterative recon: Fourier transform
    
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
    out_img(:,:,sl-(slices(1)-1),:)=gather(single(recon_cs));
    
    t(sl-(slices(1)-1)) = toc;
end %slice loop
out_img=out_img/max(abs(out_img(:)));



    