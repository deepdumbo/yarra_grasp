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
% ...for global settings & flags
doFigures               = 0;            % Display output in figures
doGpu                   = 0;            % Use GPU where implemented
slices                  = 0;            % Select which slice(s) to reconstruct. 0 = all.

% ...for selecting which parts of the code to run
doLoadData              = 1;            % Load data
doMotionDet             = 1;            % Perform motion detection
doFtz                   = 1;            % Perform FT in z-direction
doGridding              = 0;            % Perform gridding reconstruction
doUnstreaking           = 0;            % Perform coil unstreaking
doCoilCompression       = 1;            % Perform coil compression
doXdGridding            = 0;            % Perform respiratory motion-resolved griddig reconstruction
doXdGrasp               = 1;            % Perform iterative recon using XD-GRASP
doGrog                  = 0;            % Perform iterative recon using GROG
doDicomWrite            = 0;            % Write result as dicom files
 
% ...for reading raw data
loadTwixFile            = true;         % Read file meta data from pre-saved file 'twixData.mat'?
saveTwixFile            = false;        % Save a file twixData.mat after reading in raw data
channels                = 0;            % Select which coil channels to read, for testing purposes/reducing memory load. 0 = all.
spokes                  = 0;        % Select which spokes to read, for testing purposes/reducing memory load. 0 = all.
kzlines                 = 0;            % Select which kz line(s) to read from the data. 0 = all.
removeOS                = false;        % Remove readout oversampling to reduce memory load

% ...for gridding reconstruction
doSaveGridding          = 0;

% ...for coil unstreaking
n1=100;                                 % Number of spokes to generate artifact free images
n2=40;                                  % Number of spokes to generate images with streaking artifacts

% ...for motion detection
nSpokesMotionDet = spokes;         % Number of spokes used for motion detection

% ...for coil compression
ncc=8;                                  % Select to how many coil channels the data should be compressed

% ...for respiratory resolved recon
nresp=4;                                % number of respiratory phases
nLinDyn=0;                              % number of spokes per dynamic time point. 0 = all (ie single time point)
grogType = 'xdgrasp';                   % type of GROG recon: currently 'xdgrasp' or 'dce'
bas=512;                                % Number of columns/rows in cartesian-interpolated (ie GROG) k-space data & output image. 

% ... for writing data as Dicom 
dicomWriteMethod = 0;
out_path = './out';
fileName = 'grogGraspXd';

%% Load data
% Adapted from iterativeradial_multicoil.m (NYU Demo)

if doLoadData
    disp('Loading data...');
    % Thorax in vivo data
    % fileName = 'O:\MarMaa\Data\MR3_Prisma\180329_RaveTest_BioMR0009\meas_MID00200_FID04481_t1_rave_fs_tra_iso1_2.dat';
    % fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180329_RaveTest_BioMR0009/meas_MID00200_FID04481_t1_rave_fs_tra_iso1_2.dat';
    % fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180416_MINAT05/meas_MID00156_FID10853_t1_rave_fs_tra_iso1_2.dat';
    fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180813_MINAT20/meas_MID00256_FID54020_t1_rave_fs_tra_iso1_2.dat';
    
    % Abdomen in vivo data
    % fileName = 'O:\MarMaa\Data\MR3_Prisma\180329_RaveTest_BioMR0009\meas_MID00214_FID04495_t1_rave_fs_tra_iso1_2.dat';
    % fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180329_RaveTest_BioMR0009/meas_MID00214_FID04495_t1_rave_fs_tra_iso1_2.dat';    
    
    % Phantom data
    % fileName='O:\MarMaa\Data\MR3_Prisma\180322_RaveTest_Phantom\meas_MID00228_FID01617_RAVE';
    
    
    % ## Read metadata from the Siemens TWIX file
    if loadTwixFile
        load twixDataMinat20                               % Header info already read in to speed things up
    else
        twix=mapVBVD(fileName);
        
        % For VD software, the TWIX file may contain adjustment data. We only want
        % to look at the image data for now.
        twix = twix{2};
        if saveTwixFile
            save twixDataMinat20.mat twix
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
    
%     if nSpokesMotionDet == 0
        nSpokesMotionDet = length(spokes);
%         nSpokesMotionDet = twix.image.NLin;
%     end
    
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
%     Res_Signal=MC_Step3(Res_Signal,ZIP1);
    
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
%     if doGpu
%         kdata = gpuArray(single(zeros(size(rawdata))));
%     else
        kdata = single(zeros(size(rawdata)));
%     end
    for c = 1:nc
        % Fill up the missing partial Fourier data?
        % See post: https://www.researchgate.net/post/Whats_wrong_with_my_MATLAB_implementation_of_Partial_Fourier_construction_based_on_conjugate_symmetry
        
        % Inverse Fourier transform along z
        kdata(:,c,:,:) = fftshift(ifft(ifftshift(rawdata(:,c,:,:),4),[],4),4);          % Unsure whether order of fftshift and ifftshift is correct here (or whether that matters)
    end
    
    % delete the slices we don't need anymore
    if slices ~= 0
        kdata = kdata(:,:,:,slices);
    end
    
    kdata = gather(kdata);
    
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
    Traj=Trajectory_GoldenAngle(ntviews,nx,1,[0 0]);
    
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
    

%% Perform motion resolved gridding reconstruction
% Based on the demo file 'iterative_multicoil.m' by Tobias Block
if doXdGridding
    % Outline of the XD-Gridding procedure:
    % - Permute k-space to something handy
    % - Get number of lines per resp phase
    % - Loop over slices:
    % -- Select k-space for this slice (FFT has been done in slice dir)
    % -- Get sizes
    % -- Generate full trajectory (can this be place outside slice loop?)
    % -- Estimate coil sensitivities
    % -- Sort data:
    % --- Split trajectory, k-space data and respiratory signal in nt dynamic time points & throw away the last few lines that don't add up to a full dynamic
    % --- Sort trajectory and k-space data according to the respiratory motion signal in each dynamic time point
    % -- Calculate density compensation for each dynamic time point
    % -- Loop over coils
    % --- Perform gridding recon for each coil
    % -- Add up coil elements
    % -- End of slice loop
    % - Normalize image
    % - Permute back k-space data
    
    disp('Performing XD gridding reconstruction...');
    % Permute dimensions of kdata to [nx,ntviews,nz,nc]
    % This allows for easier sorting on ntviews without losing coil info?
    kdata = permute(kdata, [1,3,4,2]);
    [nx,ntviews,nz,nc]=size(kdata);
    
    % Set some initial parameters
    if slices == 0
        slices = 1:nz;
    end
%     nLinDyn = nLinDyn;
%     if nLinDyn == 0
%         nLinDyn = length(spokes);
%     end
%     nt=floor(ntviews/nLinDyn);              % Number of dynamic time points

    % Generate trajectory
    Traj=Trajectory_GoldenAngle(ntviews,nx,1,[0 0]);
    Traj=Traj*(bas/nx);

    % Calculate number of spokes per respiratory phase & crop data
    nSpokesResp = floor(ntviews/nresp);
    Traj=Traj(:,1:nSpokesResp*nresp);
    kdata=kdata(:,1:nSpokesResp*nresp,:,:);
    Res_Signal=Res_Signal(1:nSpokesResp*nresp);
%     [nx,ntviews,nt,nc] = size(kdata1);
            
    % Sort trajectory according to the respiratory motion signal
    [~,index]=sort(Res_Signal,'descend');
    Traj=Traj(:,index);
    
    % Prepare the NUFFT
    fftSize=[nx,nx];
    %         (k,w,phase,shift,imSize, mode)
    for rPhase = 1:nresp
        FT(rPhase) = NUFFT(Traj(:,(1:nSpokesResp)+(rPhase-1)*nSpokesResp), 1, 1, [0,0], fftSize, 2);
    %             FT = nufft(Traj(:,(1:nSpokesResp)+(rPhase-1)*nSpokesResp), 1, 1, [0,0], fftSize, 2);
    end

    % Initialize output image
    out_gridding    = single(zeros(nx,nx,length(slices),nresp));
    
    tic;
    for sl=slices
        fprintf('Reconstructing slice %d of %d\n', (sl-slices(1)+1), length(slices));
        % Select k-space data of slice to be reconstructed
        kdata1 = kdata(:,:,sl-slices(1)+1,:);
        [nx,ntviews,~,nc]=size(kdata1);         % MCM: seems redundant, but we have a new number of coils since compressing!
                                                % Should be inside slice loop, variables are overwritten again below

        % Coil sensitivity estimation: no algorithm yet for gridding data
        b1 = ones(nx, nx, nc);
        
        % Sort k-space data according to the respiratory motion signal
        kdata1=kdata1(:,index,:,:);
        
        % Calculate density compensation for each resp phase 
        dcf = dcfGridding(nSpokesResp, nx);    %% TO DO: Is this correct?
        
        for rPhase = 1:nresp
            % Loop over channels and calculate image for each coil
            workImage = zeros(nx,nx,nc);
            for channel=1:nc
                % Fetch rawdata for slice, channel and resp phase and 
                % multiply with density compensation function
                workKData=dcf.*double(squeeze(kdata1(:, ...
                                                     (1:nSpokesResp)+(rPhase-1)*nSpokesResp, ...
                                                     :, ...
                                                     channel)));
                
                % Run the NUFFT
                workImage(:,:,channel)=FT(rPhase)'*(workKData);
                
                % Add squared channel image to buffer
                
            end
            workImage = abs(workImage.*workImage);
            out_gridding(:,:,sl-(slices(1)-1),rPhase) = sum(workImage,3);
%             out_gridding(:,:,sl-(slices(1)-1),rPhase) = out_gridding(:,:,sl-(slices(1)-1),rPhase) + abs(workImage.*workImage);
        end
    end % slice loop
    t = toc
    % Calculate the root (as final part of the sum-of-squares calculation)
    out_gridding = sqrt(out_gridding);
    % Normalize to 1
    out_gridding = out_gridding/max(out_gridding(:));
%     if doFigures
%         figure('Name', 'Gridding Solution'), imagescn(imresize(abs(out_gridding),   [size(out_gridding,1)*4 size(out_gridding,2)*4],      'bilinear'),[],[],[],3);
%     end
    clear workKData workImage
    save gridding_recon.mat out_gridding
    disp('...done.');
end


%% XD-GRASP Reconstruction
% Adapted from Demo_XDGRASP_NonContrast (NYU Demo provided by Li Feng)

if doXdGrasp
    disp('Performing XD-GRASP reconstruction...');
    
    kdata = permute(kdata, [1,3,4,2]);
    [nx,ntviews,nz,nc]=size(kdata);                            % MCM: seems redundant, but we have a new number of coils since compressing!
%     [nx,nc,ntviews,nz]=size(kdata);                           % Not permuted by coil compression

    % Generate trajectory
    Traj=Trajectory_GoldenAngle(ntviews,nx,1,[0 0]);
    Traj=Traj*(bas/nx);
    
    % Calculate density compensation for each resp phase 
    dcf = dcfGridding(ntviews, nx);
    
    % ## Prepare the NUFFT
    fftSize=[nx,nx];
    
    % Data sorting
    nLinRes=floor(ntviews/nresp);                               % number of lines per respiratory phase
    
    % Sort the trajectory and dcf according to respiratory motion position
    [~,index]=sort(Res_Signal,'descend');
%     index = spokes'; 
    Traj=Traj(:,index);
    dcf=dcf(:,index);
    % Divide k-space trajectory and dcf into motion states
    Traj_u  = zeros(nx,nLinRes,nresp);
    dcf_u   = Traj_u;
    for ii=1:nresp
        Traj_u(:,:,ii)=Traj(:,(ii-1)*nLinRes+1:ii*nLinRes);
        dcf_u(:,:,ii)=dcf(:,(ii-1)*nLinRes+1:ii*nLinRes);   
        FT_u(ii) = NUFFT(squeeze(Traj_u(:,:,ii)), 1, 1, [0,0], fftSize, 2);
    end
    
    
    
    % Create empty images to sum the reconstructed channels (for sum-of-squares combination)
    if slices == 0
        slices = 1:nz;
    end
    gridding_xd   = single(zeros(nx,nx,length(slices),nresp));
    gridding_mcxd = single(zeros(nx,nx,length(slices),nresp));
%     grasp_xd      = single(zeros(nx,nx,length(slices),nresp));
    
    %%%% Motion resolved gridding reconstruction
    tic;
    for sl=slices
        fprintf('Starting XD-Gridding recon slice %d of %d\n', (sl-slices(1)+1), length(slices));
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
                gridding_xd(:,:,sl-(slices(1)-1),ii) = gridding_xd(:,:,sl-(slices(1)-1),ii) + abs(workImage.*workImage);
            end
        end
    end
    
    % Calculate the root (as final part of the sum-of-squares calculation)
    gridding_xd = sqrt(gridding_xd);
    gridding_xd = gridding_xd./max(gridding_xd(:));
    tGridXd = toc;
    
    for ii=1:nresp
        ttl = sprintf('Resp state %d', ii);
        if doFigures
          figure('Name', ttl);
        end
        imagescn(imresize(gridding_xd(:,:,:,ii), [size(gridding_xd,1)*2 size(gridding_xd,2)*2], 'bilinear'),[],[],[],3);
    end
    
    
    %%%% Motion resolved gridding and CS-recon using MCNUFFT 
%     b1 = ones(nx, nx, nc);   % temporary solution: implement b1 map estimation!
%     tXdGrasp = [];
%     tic;
%     for sl=slices
%         fprintf('Starting MC-XD-Gridding recon slice %d of %d\n', (sl-slices(1)+1), length(slices));
%         % Select slice to be reconstructed
%         kdata1=squeeze(kdata(:,:,sl-(slices(1)-1),:));
%         % Sort k-space data according to respiratory motion
%         kdata1=kdata1(:,index,:);
%         % Perform density compensation
%         kdata1=kdata1.*repmat(sqrt(dcf),[1,1,nc]);
%         % Divide data into motion states
%         kdata_u = zeros(nx,nLinRes,nc,nresp);
%         for ii=1:nresp
%             kdata_u(:,:,:,ii)=kdata1(:,(ii-1)*nLinRes+1:ii*nLinRes,:);
%         end
%         param.E=MCNUFFT(double(Traj_u),double(dcf_u),double(b1));
%         param.y=double(kdata_u);
%         recon_cs = param.E'*param.y;
%         gridding_mcxd(:,:,sl-(slices(1)-1),:)=recon_cs;
%         
%         param.TV_dim1=TV_Temp;
%         param.TVWeight_dim1=max(abs(recon_cs(:)))*0.02;
%         param.TVWeight_dim2=0;
%         param.nite = 4;
%         param.display=1;
%         
% %         tic
%         for n=1:2
%             recon_cs = CSL1NlCg_XDGRASP(recon_cs,param);
%         end
% %         time=toc;
% %         time=time/60
%         grasp_xd(:,:,sl-(slices(1)-1),:)=recon_cs;
%         tXdGrasp(sl-(slices(1)-1)) = toc;
%     end
%      
%     gridding_mcxd = gridding_mcxd./max(abs(gridding_mcxd(:)));
%     tGridMcXd = toc;
%     grasp_xd = grasp_xd/max(abs(grasp_xd(:)));
%     save xdgrasp_recon.mat grasp_xd gridding_mcxd
%     
%     if doFigures
%         for ii=1:nresp
%             ttl = sprintf('MC Gridding: Resp state %d', ii);
%             figure('Name', ttl);
%             imagescn(imresize(abs(gridding_mcxd(:,:,:,ii)), [size(gridding_mcxd,1)*2 size(gridding_mcxd,2)*2], 'bilinear'),[],[],[],3);
%             ttl = sprintf('GRASP: Resp state %d', ii);
%             figure('Name', ttl);
%             imagescn(imresize(abs(grasp_xd(:,:,:,ii)), [size(gridding_mcxd,1)*2 size(gridding_mcxd,2)*2], 'bilinear'),[],[],[],3);
%         end
%     end
    
    % Permute dimensions back to default [nx,nc,ntviews,nz]
    kdata = permute(kdata, [1,4,2,3]);
    
    disp('...done.');
end


%% GROG Reconstruction
% Adapted from Demo4_Reconstruction.m from RACER-GRASP_GROG-GRASP demo
% package (NYU Demo provided by Li Feng)
% This section includes coil sensitivities estimation.
if doGrog
    disp('Performing GROG Reconstruction...');
    
    % Permute dimensions of kdata to [nx,ntviews,nz,nc]
%     kdata = permute(kdata, [1,3,4,2]);
    
    % Set some initial parameters
    if slices == 0
        slices = 1:nz;
    end
    if nLinDyn == 0
        nLinDyn = length(spokes);
    end
    nt=floor(ntviews/nLinDyn);              % Number of dynamic time points
    Nqu=min(floor(nx/2 * pi/2), ntviews);   % Nyquist spokes, please refer to the GROG-GRASP paper for more details about this.
                                            % MCM changed from 350
    
    % Initialize output image & recon time measurement
    grogGrasp_xd      = single(zeros(nx,nx,length(slices),nresp));
    tGrog = [];
    for sl=slices % loop through selected slices
        disp(sl)
        tic;
        % Select slice to be reconstructed
        kdata1=kdata(:,:,sl-(slices(1)-1),:);
        [nx,ntviews,~,nc]=size(kdata1);         % MCM: seems redundant, but we have a new number of coils since compressing!
                                                % Should be inside slice loop, variables are overwritten again below
                                               
        % Generate trajectory (For GROG only, NOT for general gridding)
        % This needs to stay within the slice loop; it gets modified during
        % sorting of data in respiratory phases
        Traj=Trajectory_GoldenAngle(ntviews,nx,0,[0 1]);        % GROG expects no normalization and reversed spoke traversal
                                                                % MCM changed from Trajectory_GoldenAngle_GROG()
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
        % motion signal in each dynamic time point
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
        nline2=floor(ntviews/nresp);                            % Number of lines per respiratory bin
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
        
        switch grogType
            case 'xdgrasp'
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
                for n=1:2
                    recon_cs = CSL1NlCg_XDGRASP_Mx(recon_cs,param);
                end
%                 recon_cs=CropImg(recon_cs,nx/2,nx/2);
                grogGrasp_xd(:,:,sl-(slices(1)-1),:)=gather(abs(single(recon_cs)));
            case 'dce'
                % Permute k-space data to correct order (depending on
                % algorithm to be used)
                kdata1=permute(kdata1,[1,2,4,5,3]);
                
                % Calculate weights for respiratory weighted reconstruction
                SoftWeight=single(zeros(bas,bas,nt,nc,nresp));
                for jj=1:nresp
                    SoftWeight(:,:,:,:,jj)=(1/4)^(jj-1);
                end
                if doGpu
                    param.SG=gpuArray(SoftWeight);
                else
                    param.SG=SoftWeight;
                end
                
                % Prepare iterative recon
                param.y=kdata1.*sqrt(Weighting);
                mask=single(kdata1~=0);
                param.E=Emat_GROG2DSG(mask,b1,Weighting);
                recon_cs=param.E'*param.y;                  % Starting point for iterative recon: Fourier transform
                
                % Set parameters for optimizer
                param.TV=TV_Temp;
                Weight1=0.04;
                param.TVWeight=max(abs(recon_cs(:)))*Weight1;
                param.nite = 8;
                param.display = 1;
                
                for n=1:3
                    recon_cs = CSL1NlCg_4DDCE(recon_cs,param);
                end
%                 recon_cs=CropImg(recon_cs,nx/2,nx/2);
                grogGrasp_xd(:,:,:,sl-(slices(1)-1))=gather(abs(single(recon_cs)));
        end
        tGrog(sl-(slices(1)-1)) = toc;
    end %slice loop
    grogGrasp_xd=grogGrasp_xd/max(grogGrasp_xd(:));
%     save result.mat grogGrasp_xd 

    % Crop image in case oversampling was used
    if ~removeOS
        grogGrasp_xd=CropImg(grogGrasp_xd,nx/2,nx/2);
    end
%   Permute dimensions back to default [nx,nc,ntviews,nz]
    kdata = permute(kdata, [1,4,2,3]);

    disp('...done.');
end

%% Plotting data
% imagesc(data(:,:,11,1)),axis image,colormap(gray), axis off

%% Writing data to Dicom slices
if doDicomWrite
    grogGrasp_xd = uint16(grogGrasp_xd./max(grogGrasp_xd(:))*(2^12 - 1));
    
    switch dicomWriteMethod
        case 0 % Even stupider: use NIFTI
            % Create 4D NIFTI structure
            nii = make_nii(grogGrasp_xd, voxSize);          % Uses NIFTI Toolbox by Jimmy Shen
            % Save it
            fName = sprintf('%s_%02d.nii', fileName);
            fName = fullfile(out_path, fName);
            save_nii(nii, fName);                           % Uses NIFTI Toolbox by Jimmy Shen
%             niftiwrite(grogGrasp_xd,fName);
            
        case 1 % The stupid way: just write them without any additional info.

            % Reshape into a 3D volume such that the multiple volumes are stacked on
            % top of eachother
            s = size(grogGrasp_xd);
            grogGrasp_xd = reshape(grogGrasp_xd, [s(1) s(2) prod(s(3:end))]);
            
            infoGrasp = dcmInitInfo();
            infoGrasp.SeriesDescription = 'GrogGrasp'
            infoGrasp = repmat(infoGrasp,size(grogGrasp_xd,3),1);
            imgNumbers = num2cell(1:length(infoGrasp));
            [infoGrasp.InstanceNumber] = imgNumbers{:};

            % Write 2D slices to dcm files
            for i=1:10 %:size(grogGrasp_xd,3)
%                 fName = sprintf('%s_%03d.ima', fileName, i);
                fName = sprintf('%s_%03d.dcm', fileName, i);
                fName = fullfile(out_path, fName);
                dicomwrite(squeeze(grogGrasp_xd(:,:,i)), fName, infoGrasp(i));
            end
        case 2 % Get dicom tags from example series generated by scanner
            % Load 'example' series to copy Dicom tags from
            [~, dcmInfo, ~] = dcmReadSingleVolume(dcmExampleDir, '');
            
            % Create new Dicom info objects
%             infoGrasp        = { dcmInfo{1:length(dcmInfo)/length(bValues_L)} };             % 1 info object for each slice                                        % Necessary for new SeriesInstanceUID
            infoGrasp        = dcmInfo(1:length(slices));             % Use as many info objects as we have reconstructed slices. Making life easy: we'e not selecting _which_ info objects for correct positioning
            infoGrasp        = repmat(infoGrasp,1,nresp);             % 1 info object for each slice at each time point
            infoAddGrasp     = dcmInitInfo(dcmInfo{1});
            
            % Update attributes
            infoAddGrasp.SeriesNumber = 1001;
%             TODO: Add the necessary tags/attributes here. Which ones?
            % - Some tag to make it multi-frame. I can probably choose...
            % Something to do with Measurements? Perhaps I can check in a
            % dynamic series somewhere.
            % - Image comment
            % - Image type...
            % - Series description
            % - ...
%             infoAddGrasp.SeriesDescription = 'ADC2_0000-1200';
            
            % Apply updates to all slices
            infoGrasp  = dcmUpdateInfo(infoGrasp, infoAddGrasp);
            
            % Write slices to file
            directory = out_path;
            if ~exist(directory, 'dir')
                mkdir(directory);
            end
            dcmWriteDwi(uint16(grogGrasp_xd), infoGrasp, out_file, directory);
    end
end


%% Next steps?
% - Generate corresponding k-space trajectory
% - Calculating coil sensitivity maps / Coil compression
% 
% - Reconstruct
% - ...
% - Write slices to disk

% Plot data
