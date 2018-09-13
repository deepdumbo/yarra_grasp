%% Script for performing GROG-GRASP reconstruction of multi-slice, single time point RAVE acquisition

clear
clc

% ## Include packages in subfolders
addpath('../mapVBVD/');                                             % for reading TWIX files
addpath('../NYU Demos/IterativeRadial_v2/');                        % for the NUFFT
addpath('../NYU Demos/IterativeRadial_v2/nufft');                   % for the NUFFT
addpath('../NYU Demos/IterativeRadial_v2/utilities');               % for the NUFFT
addpath('../NYU Demos/IterativeRadial_v2/graph');                   % for the NUFFT
addpath('../NYU Demos/IterativeRadial_v2/imagescn_R2008a/');        % for plotting images
addpath('../NYU Demos/XDGRASP_Demo/Code');
addpath('../NYU Demos/XDGRASP_Demo/Code/nufft_files');
addpath('../NYU Demos/Demo_RACER-GRASP_GROG-GRASP/utils');          % for GROG-specific tools
addpath('../NYU Demos/Demo_RACER-GRASP_GROG-GRASP/nufft_files');    % for GROG-specific tools
% addpath('../NYU Demos/IterativeRadial_v2/poblano_toolbox_1.1');  % for numerical optimization

% Set parameters:
% ...for selecting which parts of the code to run
doLoadData              = 1;            % Load data
doMotionDet             = 1;            % Perform motion detection
doFtz                   = 1;            % Perform FT in z-direction
doGridding              = 1;            % Perform gridding reconstruction
doUnstreaking           = 0;            % Perform coil unstreaking
doCoilCompression       = 1;            % Perform coil compression
doXdGrasp               = 0;            % Perform iterative recon using XD-GRASP
doGrog                  = 1;            % Perform iterative recon using GROG
doFigures               = 1;            % Display output in figures
 
% ...for reading raw data
loadTwixFile            = true;         % Read file meta data from pre-saved file 'twixData.mat'?
saveTwixFile            = false;        % Save a file twixData.mat after reading in raw data
channels                = 0;            % Select which coil channels to read, for testing purposes/reducing memory load. 0 = all.
spokes                  = 1:110;        % Select which spokes to read, for testing purposes/reducing memory load. 0 = all.
kzlines                 = 0;            % Select which kz line(s) to read from the data. 0 = all.
removeOS                = false;        % Remove readout oversampling to reduce memory load

% ...for gridding reconstruction
slices                  = 120:130;      % Select which slice(s) to reconstruct. 0 = all.
doSaveGridding          = 0;

% ...for coil unstreaking
n1=100;                                 % Number of spokes to generate artifact free images
n2=40;                                  % Number of spokes to generate images with streaking artifacts

% ...for motion detection
nSpokesMotionDet = max(spokes);         % Number of spokes used for motion detection

% ...for coil compression
ncc=8;                                  % Select to how many coil channels the data should be compressed

% ...for respiratory resolved recon
ntres=4;                                % number of respiratory phases
nLinDyn=length(spokes);                 % number of spokes per dynamic time point. Changed from 96



%% Load data
% Adapted from iterativeradial_multicoil.m (NYU Demo)

if doLoadData
    disp('Loading data...');
    % Thorax in vivo data
    % fileName = 'O:\MarMaa\Data\MR3_Prisma\180329_RaveTest_BioMR0009\meas_MID00200_FID04481_t1_rave_fs_tra_iso1_2.dat';
    % fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180329_RaveTest_BioMR0009/meas_MID00200_FID04481_t1_rave_fs_tra_iso1_2.dat';
    
    % Abdomen in vivo data
    % fileName = 'O:\MarMaa\Data\MR3_Prisma\180329_RaveTest_BioMR0009\meas_MID00214_FID04495_t1_rave_fs_tra_iso1_2.dat';
    % fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180329_RaveTest_BioMR0009/meas_MID00214_FID04495_t1_rave_fs_tra_iso1_2.dat';
    % fileName = 'O:\MarMaa\Data\MR3_Prisma\180416_MINAT05\meas_MID00156_FID10853_t1_rave_fs_tra_iso1_2.dat';
    fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180416_MINAT05/meas_MID00156_FID10853_t1_rave_fs_tra_iso1_2.dat';
    
    % Phantom data
    % fileName='O:\MarMaa\Data\MR3_Prisma\180322_RaveTest_Phantom\meas_MID00228_FID01617_RAVE';
    
    
    % ## Read metadata from the Siemens TWIX file
    if loadTwixFile
        load twixData1                               % Header info already read in to speed things up
    else
        twix=mapVBVD(fileName);
        
        % For VD software, the TWIX file may contain adjustment data. We only want
        % to look at the image data for now.
        twix = twix{2};
        if saveTwixFile
            save twixData1.mat twix
        end
    end
    
    % Set some parameters for reading the raw data
    twix.image.flagRemoveOS = removeOS;         % Remove readout oversampling to reduce memory load
    
    if channels == 0
        channels = 1:twix.image.NCha;
    end
    if spokes == 0
        spokes = 1:twix.image.NLin;
    end
    if kzlines == 0
        kzlines = 1:twix.image.NPar;
    end
    
    % Get the k-space data. Data comes in as [samples,channels,spokes,kzlines]
    rawdata = twix.image(:,channels,spokes,kzlines);     % This makes no sense in multi-slice data, does it? Should it not be FFT'd along z first?
    
    disp('...done.');
end

%% Motion detection
% Adapted from Demo2_MotionDetection.m (NYU Demo provided by Li Feng)
if doMotionDet
    disp('Performing motion detection...');
    [nx nc ntviews nz]=size(rawdata);
    
    % then, you can get your ZIP with
    ZIP = squeeze(rawdata(nx/2+2,:,:,:));
    
    % Respiratory motion detection
    ZIP=permute(ZIP,[3,2,1]);                                   % MCM changed from [2,1,3]
    ZIP=abs(fftshift(ifft(ZIP,400,1),1)); % FFT and interpolation along the kz dimension
    
    %Normalization of each projection in each coil element
    ZIP=ProjNorm(ZIP);%Normalization includes temporal smoothing
    if doFigures
        figure,imagesc(abs(ZIP(:,:,15))),axis image, axis off, colormap(gray),title('Respiratory Motion')
    end
    
    %There are 3 steps to generate a respiratory motion signal, as shown below    
    % STEP 1: find the coil elements with good representation of respiratory motion
    %         from the late enhancement spokes
    [Coil,Res_Signal_Post]=MC_Step1(ZIP,nSpokesMotionDet);
    
    %STEP 2: Estimate motion signal using PCA from the concatated coil elements
    %Those coil elements were selected in the first step
    [SI,corrm,Res_Signal,ZIP1]=MC_Step2(ZIP,Coil,nSpokesMotionDet,Res_Signal_Post);
    
    %Step 3: You noticed that the signal is not flat, due to the contrast
    %injection. So, now let's estimate the envelop of the signal and substract it
    Res_Signal=MC_Step3(Res_Signal,ZIP1);
    
    %save the estimated motion signal
    save Res_Signal.mat Res_Signal
    
    % Free up some space
    clear ZIP ZIP1 
    
    disp('...done.');
end

%% Fourier transform in Z-direction
% So that each 'star' can be reconstructed separately
% This can be done for each coil separately to save memory
% How to treat 'partial fourier' in that direction?
if doFtz
    disp('Performing Fourier Transform in z-direction...');
    [nx, nc, ntviews, nkz] = size(rawdata);
    kdata = single(zeros(size(rawdata)));
    for c = 1:nc
        % Fill up the missing partial Fourier data?
        
        % Inverse Fourier transform along z
        kdata(:,c,:,:) = fftshift(ifft(ifftshift(rawdata(:,c,:,:),4),[],4),4);          % Unsure whether order of fftshift and ifftshift is correct here (or whether that matters)
    end
    
    % delete the slices we don't need anymore
    if slices ~= 0
        kdata = kdata(:,:,:,slices);
    end
    
    % Free up some memory
    clear rawdata
    
    disp('...done.');
end


%% Coil unstreaking
% Adapted from NYU Demo Demo1_Unstreaking.m
% Should be performed before coil compression
if doUnstreaking
    disp('Performing Coil Unstreaking...');
    % MCM: Permute the k-space data first: what routine expects this?
    % Permute dimensions to prepare for Coil Unstreaking
    kdata = permute(kdata,[1,3,4,2]);
    
    [nx,ntviews,nz,nc]=size(kdata);
    % Note that the third dimension is z, NOT kz. (A FFT was performed already)
    
    % Generate trajectory (For GROG only, NOT for general gridding)
    Traj=Trajectory_GoldenAngle(ntviews,nx);
    
    % Calculate density compensation
    DensityComp = dcfGridding(ntviews, nx);
    
    Ref=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n1,nx));                 % artifact-free image
    Img=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n2,nx)*n1/n2);           % image with streaks
    
    %As described in the paper, the Diff image is calculated as the 2x FOV
    Diff=abs(Ref-Img);
    
    %The Ref and Img are then cropped to the 1x FOV
    Ref=Ref(nx/4+1:end-nx/4,nx/4+1:end-nx/4,:,:);
    Img=Img(nx/4+1:end-nx/4,nx/4+1:end-nx/4,:,:);
    
    data=sos(Ref,4);
    data(:,:,:,2)=sos(Img,4);
    
    %calculating the streak ratio
    clear StreakRatio
    for ii=1:nc
        StreakRatio(ii,1)=norm(col(Diff(:,:,:,ii)))/norm(col(Ref(:,:,:,ii)));
    end
    StreakRatio=StreakRatio/min(StreakRatio);
    if doFigures
        figure,plot(StreakRatio) % plot streak ratio for each coil element
    end
    
    %find the coil elements whose streak ratio greater than 1.3
    %unstreaking is performed only for these coils, as described in the paper
    Coil=find(StreakRatio>1.3)';
    
    % Do unstreaking
    StreakRatio=repmat(StreakRatio(Coil),[1,nx,nz,ntviews]);
    StreakRatio=permute(StreakRatio,[2,4,3,1]);
    kdata(:,:,:,Coil)=kdata(:,:,:,Coil)./StreakRatio;
    
    % reconstruct images again for comparison
    Ref=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n1,nx/2));               % artifact-free image
    Img=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n2,nx/2)*n1/n2);         % image with streaks
    
    data(:,:,:,3)=sos(Ref,4);
    data(:,:,:,4)=sos(Img,4);
    data=data./max(data(:));
    
    % Note that the images displayed below are before iterative reconstruction
    % Thus, they both have strearking artifacts.
    % However, note that the lower image has significantly less streaks
    subplot(2,1,1)
    imagesc(abs(data(:,:,8,2))),axis off, axis square, colormap(gray),title('before unstreaking')
    subplot(2,1,2)
    imagesc(abs(data(:,:,8,4))),axis off, axis square, colormap(gray),title('after unstreaking')
    
    % Permute k-space data back to its 'default': [samples,coils,spokes,partitions]
%     was: [samples,spokes,partitions,coils]
    kdata = permute(kdata, [1,4,2,3]);
    
    disp('...done.');
end

%% Coil compression
% Adapted from Demo1_Unstreaking.m (NYU Demo provided by Li Feng)
% 
% TODO: check whether this needs to be done on _all_ slices as here, or
% only on the selected slices. That would save tons of memory space.
% TODO2: check whether this step can be performed before gridding recon.
if doCoilCompression
    disp('Performing Coil Compression...');
    
    % Permute dimensions:
    % from [nx,nc,ntviews,nz]
    % to   [nx,ntviews,nz,nc]
    kdata = permute(kdata,[1,3,4,2]);
    [nx,ntviews,nz,nc]=size(kdata);

    D=reshape(kdata,nx*nz*ntviews,nc);
    [U,S,V]=svd(D,'econ');
    kdata=single(reshape(D*V(:,1:ncc),nx,ntviews,nz,ncc));
    % save Data/kdata_Unstreaking_CoilCompression.mat kdata
    
    % Free up some memory
    clear D U S V
    
    disp('...done.');
    
    % Permute dimensions back to default [nx,nc,ntviews,nz]
    kdata = permute(kdata, [1,4,2,3]);
    
    
end

%% Perform gridding reconstruction
% Based on the demo file 'iterative_multicoil.m' by Tobias Block
if doGridding
    disp('Performing gridding reconstruction...');
    [nx,nc,ntviews,nz]=size(kdata);
    
    % Calculate density compensation (DCF)
    dcf = dcfGridding(ntviews, nx);
    
    % Generate trajectory
    Traj=Trajectory_GoldenAngle(ntviews,nx,1,0);
    % TODO: Check flipping or gradient corrections
    
    % ## Prepare the NUFFT
    fftSize=[nx,nx];
    %         (k,w,phase,shift,imSize, mode)
    FT = NUFFT(Traj(:,1:end), 1, 1, [0,0], fftSize, 2);
    
    % Create image to sum the gridded channels (for sum-of-squares combination)
    if slices == 0
        slices = 1:nz;
    end
    gridding=single(zeros(nx,nx,length(slices)));
    for sl=slices        
        % Loop over channels and calculate image for each coil
        for channel=1:nc
            % Fetch rawdata for slice and channel and multiply with density
            % compensation function
            workKData=dcf.*double(squeeze(kdata(:,channel,:,sl-(slices(1)-1))));
            
            % Run the NUFFT
            workImage=FT'*(workKData);
            
            % Add squared channel image to buffer
             gridding(:,:,sl-(slices(1)-1)) = gridding(:,:,sl-(slices(1)-1)) + abs(workImage.*workImage);
        end
    end
    % Calculate the root (as final part of the sum-of-squares calculation)
    gridding = sqrt(gridding);
    if doFigures
        figure('Name', 'Gridding Solution'), imagescn(imresize(abs(gridding),   [size(gridding,1)*4 size(gridding,2)*4],      'bilinear'),[],[],[],3);
    end
    clear workKData workImage
    save gridding_recon.mat gridding
    disp('...done.');
end
    

%% XD-GRASP Reconstruction
% Adapted from Demo_XDGRASP_NonContrast (NYU Demo provided by Li Feng)

if doXdGrasp
    kdata = permute(kdata, [1,3,4,2]);
    [nx,ntviews,nz,nc]=size(kdata);                            % MCM: seems redundant, but we have a new number of coils since compressing!
%     [nx,nc,ntviews,nz]=size(kdata);                           % Not permuted by coil compression
    
    % Data sorting
    nLinRes=floor(ntviews/ntres);                               % number of lines per respiratory phase
    
    % Sort the trajectory and dcf according to respiratory motion position
    [~,index]=sort(Res_Signal,'descend');
%     index = spokes'; 
    Traj=Traj(:,index);
    dcf=dcf(:,index);
    % Divide k-space trajectory and dcf into motion states
    Traj_u  = zeros(nx,nLinRes,ntres);
    dcf_u   = Traj_u;
    for ii=1:ntres
        Traj_u(:,:,ii)=Traj(:,(ii-1)*nLinRes+1:ii*nLinRes);
        dcf_u(:,:,ii)=dcf(:,(ii-1)*nLinRes+1:ii*nLinRes);        
    end
    
    % ## Prepare the NUFFT
    fftSize=[nx,nx];
    
    % Create empty images to sum the reconstructed channels (for sum-of-squares combination)
    if slices == 0
        slices = 1:nz;
    end
%     gridding_xd = single(zeros(nx,nx,length(slices),ntres));
    gridding_mcxd = single(zeros(nx,nx,length(slices),ntres));
    grasp_xd      = single(zeros(nx,nx,length(slices),ntres));
    
    %%%% Motion resolved gridding reconstruction
%     for sl=slices
%         % Select slice to be reconstructed
%         kdata1=squeeze(kdata(:,:,sl-(slices(1)-1),:));
%         % Sort k-space data according to respiratory motion
%         kdata1=kdata1(:,index,:);
%         % Divide data into motion states
%         kdata_u = zeros(nx,nLinRes,nc,ntres);
%         for ii=1:ntres
%             kdata_u(:,:,:,ii)=kdata1(:,(ii-1)*nLinRes+1:ii*nLinRes,:);
% %             kdata_u(:,:,:,ii)=kdata1(:,(1:nLinRes)+ii-1,:);
%             %         (k,w,phase,shift,imSize, mode)
%             FT_u = NUFFT(squeeze(Traj_u(:,:,ii)), 1, 1, [0,0], fftSize, 2);
%             for channel=1:nc
%                 % Fetch k-data for slice and channel and multiply with density
%                 % compensation function
%                 workKData=dcf_u(:,:,ii).*double(squeeze(kdata_u(:,:,channel,ii)));              
%                 % Run the NUFFT
%                 workImage=FT_u'*(workKData);                
%                 % Add squared channel image to buffer
%                 gridding_xd(:,:,sl-(slices(1)-1),ii) = gridding_xd(:,:,sl-(slices(1)-1),ii) + abs(workImage.*workImage);
%             end
%         end
%     end
%         
%     % Calculate the root (as final part of the sum-of-squares calculation)
%     gridding_xd = sqrt(gridding_xd);
%     gridding_xd = gridding_xd./max(gridding_xd(:));
%     for ii=1:ntres
%         ttl = sprintf('Resp state %d', ii);
%         if doFigures
%           figure('Name', ttl);
%         end
%         imagescn(imresize(gridding_xd(:,:,:,ii), [size(gridding,1)*2 size(gridding,2)*2], 'bilinear'),[],[],[],3);
%     end
%     
    
    %%%% Motion resolved gridding and CS-recon using MCNUFFT 
    b1 = ones(nx, nx, nc);   % temporary solution: implement b1 map estimation!
    for sl=slices
        % Select slice to be reconstructed
        kdata1=squeeze(kdata(:,:,sl-(slices(1)-1),:));
        % Sort k-space data according to respiratory motion
        kdata1=kdata1(:,index,:);
        % Perform density compensation
        kdata1=kdata1.*repmat(sqrt(dcf),[1,1,nc]);
        % Divide data into motion states
        kdata_u = zeros(nx,nLinRes,nc,ntres);
        for ii=1:ntres
            kdata_u(:,:,:,ii)=kdata1(:,(ii-1)*nLinRes+1:ii*nLinRes,:);
        end
        param.E=MCNUFFT(double(Traj_u),double(dcf_u),double(b1));
        param.y=double(kdata_u);
        recon_cs = param.E'*param.y;
        gridding_mcxd(:,:,sl-(slices(1)-1),:)=recon_cs;
        
        param.TV_dim1=TV_Temp;
        param.TVWeight_dim1=max(abs(recon_cs(:)))*0.02;
        param.TVWeight_dim2=0;
        param.nite = 4;
        param.display=1;
        
        tic
        for n=1:2
            recon_cs = CSL1NlCg_XDGRASP(recon_cs,param);
        end
        time=toc;
        time=time/60
        grasp_xd(:,:,sl-(slices(1)-1),:)=recon_cs;
    end
     
    gridding_mcxd = gridding_mcxd./max(abs(gridding_mcxd(:)));
    grasp_xd = grasp_xd/max(abs(grasp_xd(:)));
    save xdgrasp_recon.mat grasp_xd gridding_mcxd
    
    if doFigures
        for ii=1:ntres
            ttl = sprintf('MC Gridding: Resp state %d', ii);
            figure('Name', ttl);
            imagescn(imresize(abs(gridding_mcxd(:,:,:,ii)), [size(gridding_mcxd,1)*2 size(gridding_mcxd,2)*2], 'bilinear'),[],[],[],3);
            ttl = sprintf('GRASP: Resp state %d', ii);
            figure('Name', ttl);
            imagescn(imresize(abs(grasp_xd(:,:,:,ii)), [size(gridding_mcxd,1)*2 size(gridding_mcxd,2)*2], 'bilinear'),[],[],[],3);
        end
    end
end


%% GROG Reconstruction
% Take this from Demo4_Reconstruction.m : TODO!!
% The demo also includes coil sensitivities estimation
if doGrog
    disp('Performing GROG Reconstruction...');
    
    kdata = permute(kdata, [1,3,4,2]);
    if slices == 0
        slices = 1:nz;
    end
%     nz = npar;
    bas=360;% the 2x full FOV is 512 and we only reconstruct 360 here.
    for sl=slices % loop through selected slices
        % Select slice to be reconstructed
        kdata1=squeeze(kdata(:,:,sl-(slices(1)-1),:));
        
%         kdata=kdata1(:,:,sl-(slices(1)-1),:);
        [nx,ntviews,~,nc]=size(kdata1);      % MCM: seems redundant, but we have a new number of coils since compressing!
        % Isn't this wrong? should it be put outside the slice loop?
        
        %Generating trajectory (For GROG only, NOT for general gridding)
        Traj=Trajectory_GoldenAngle_GROG(ntviews,nx);%trajectory
        Traj=Traj*(bas/nx);
        
        %GROG weights calibration
        [Gx,Gy] = GROG.get_Gx_Gy(kdata1,Traj);
        
        %Coil sensitivities estimation
        G = GROG.init(kdata1,Traj,Gx,Gy,0); % GROG initialization
        kref = GROG.interp(kdata1,G,1,0); % GROG interpolation
        kref=CropImg(kref,bas,bas);
        ref=gather(squeeze(ifft2c_mri(kref)));
        b1=adapt_array_2d(ref);clear ref
        b1=single(b1/max(abs(b1(:))));
        
        nt=floor(ntviews/nLinDyn);            % MCM: Number of dynamic time points
        Nqu=nLinDyn;                          % Nyquist spokes, please refer to the GROG-GRASP paper for more details about this.
        % MCM changed from 350
        
        G = GROG.init(kdata1(:,end-Nqu+1:end,:,:),Traj(:,end-Nqu+1:end),Gx,Gy,0);
        D=reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1))]);
        DCF=CropImg(D,bas,bas);
        DCF=repmat(DCF,[1,1,nt,nc]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Data sorting
        Traj=reshape(Traj(:,1:nt*nLinDyn),[nx,nLinDyn,nt]);
        kdata1=reshape(kdata1(:,1:nt*nLinDyn,:,:),[nx,nLinDyn,nt,nc]);
        Res_Signal_tmp=reshape(Res_Signal(1:nt*nline),[nline,nt]);
        [nx,ntviews,nt,nc] = size(kdata1);
        
        % Sorting according to the respiratory motion signal
        for ii=1:nt
            [~,index]=sort(Res_Signal_tmp(:,ii),'descend');
            Traj(:,:,ii)=Traj(:,index,ii);
            kdata1(:,:,ii,:)=kdata1(:,index,ii,:);
        end
        
        [nx,ntviews,nt,nc] = size(kdata1);
        G = GROG.init(kdata1,Traj,Gx,Gy,0);
        DCF_U=reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1)),nt]);
        DCF_U=CropImg(DCF_U,bas,bas);
        
        %     ntres=4;% 4 respiratory phases
        nline2=floor(ntviews/ntres);
        kdata1=reshape(kdata1,[nx,nline2,ntres,nt,nc]);
        kdata1=reshape(kdata1,[nx,nline2,ntres*nt,nc]);
        Traj=reshape(Traj,[nx,nline2,ntres,nt]);
        Traj=reshape(Traj,[nx,nline2,ntres*nt]);
        
        % Do GROG; the preinterpolation of radial data into a Cartesian grid
        G = GROG.init(kdata1,Traj,Gx,Gy,0);
        kdata1 = GROG.interp(kdata1,G,1,0);
        kdata1=CropImg(kdata1,bas,bas);
        
        % calculate GROG weights, as described in the paper
        DCF_U=repmat(DCF_U,[1,1,1,nc]);
        DCF_U(find(DCF_U==0))=1;
        Weighting=DCF./DCF_U;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        kdata1=reshape(kdata1,[bas,bas,ntres,nt,nc]);
        kdata1=permute(kdata1,[1,2,4,5,3]);
        
        SoftWeight=single(zeros(bas,bas,nt,nc,ntres));
        for jj=1:ntres
            SoftWeight(:,:,:,:,jj)=(1/4)^(jj-1);
        end
        
        mask=single(kdata1~=0);
        param.SG=gpuArray(SoftWeight);
        Weighting=repmat(Weighting,[1 1 1 1 ntres]);
        
        param.y=kdata1.*sqrt(Weighting);
        param.E=Emat_GROG2DSG(mask,b1,Weighting);
        recon_cs=param.E'*param.y;
        
        param.TV=TV_Temp;
        Weight1=0.04;
        param.TVWeight=max(abs(recon_cs(:)))*Weight1;
        param.nite = 8;param.display = 1;
        
        clc
        for n=1:3
            recon_cs = CSL1NlCg_4DDCE(recon_cs,param);
        end
        recon_cs=CropImg(recon_cs,nx/2,nx/2);
        data(:,:,:,sl-(slices(1)-1))=gather(abs(single(recon_cs)));
    end
    data=data/max(data(:));
    save result.mat data 

    disp('...done.');
end

%% Plotting data
% imagesc(data(:,:,11,1)),axis image,colormap(gray), axis off

%% Writing data to Dicom slices

%% Next steps?
% - Generate corresponding k-space trajectory
% - Calculating coil sensitivity maps / Coil compression
% 
% - Reconstruct
% - ...
% - Write slices to disk

% Plot data
