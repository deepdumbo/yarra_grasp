function out_gridding = grog_xdgrasp(in_path, in_file, out_path, temp_path, pars)

% Function for performing a respiratory resolved XD-GRASP reconstruction
% using GROG pre-interpolation, based on NYU demo scripts by Li Feng.
% Follows the convention of Yarra modules for easy integration.
% The input parameter pars can be used to pass optional parameters:
% within Yarra, this appears to be done through a separate config file. We
% will delete it as appropriate.
%
% Written by Marnix Maas (Marnix.Maas@radboudumc.nl), May 2018


%% Set initial parameters
clc

% ## Include packages in subfolders
addpath('../mapVBVD/');                                             % for reading TWIX files
addpath('../NYU Demos/IterativeRadial_v2/');                        % for the NUFFT
addpath('../NYU Demos/IterativeRadial_v2/nufft');                   % for the NUFFT
addpath('../NYU Demos/IterativeRadial_v2/utilities');               % for the NUFFT

addpath('../NYU Demos/XDGRASP_Demo/Code');                          % for iterative recon code
addpath('../NYU Demos/XDGRASP_Demo/Code/nufft_files');
addpath('../NYU Demos/Demo_RACER-GRASP_GROG-GRASP/utils');          % for GROG-specific tools
addpath('../NYU Demos/Demo_RACER-GRASP_GROG-GRASP/nufft_files');    % for GROG-specific tools

if nargin<5
    pars = initReconPars('grog_xdgrasp');
end


%% Load data
% Adapted from iterativeradial_multicoil.m (NYU Demo)

if pars.doLoadData
    disp('Loading data...');
    % ## Read metadata from the Siemens TWIX file
    if pars.loadTwixFile
        load twixData1                               % Header info already read in to speed things up
    else
        fileName = fullfile(in_path,in_file);
        twix=mapVBVD(fileName);
    
        % For VD software, the TWIX file may contain adjustment data. We only want
        % to look at the image data for now.
        twix = twix{2};
        if pars.saveTwixFile
            save twixData1.mat twix
        end
    end
    
    % Set some parameters for reading the raw data
    twix.image.flagRemoveOS = pars.removeOS;         % Remove readout oversampling to reduce memory load
    
    channels = pars.channels;
    if channels == 0
        channels = 1:twix.image.NCha;
    end
    spokes = pars.spokes;
    if spokes == 0
        spokes = 1:twix.image.NLin;
    end
    kzlines = 1:twix.image.NPar;
    
    % Get the k-space data. Data comes in as [samples,channels,spokes,kzlines]
    rawdata = twix.image(:,channels,spokes,kzlines); 
    
    disp('...done.');
end

%% Motion detection
% Adapted from Demo2_MotionDetection.m (NYU Demo provided by Li Feng)
if pars.doMotionDet
    disp('Performing motion detection...');
    if pars.nSpokesMotionDet == 0
        pars.nSpokesMotionDet = length(spokes);
    end						
    [nx nc ntviews nkz]=size(rawdata);
    
    ZIP = squeeze(rawdata(nx/2+2,:,:,:));
    
    % Respiratory motion detection
    ZIP=permute(ZIP,[3,2,1]);                                   % MCM changed from [2,1,3]
    ZIP=abs(fftshift(ifft(ZIP,400,1),1)); % FFT and interpolation along the kz dimension
    
    %Normalization of each projection in each coil element
    ZIP=ProjNorm(ZIP);%Normalization includes temporal smoothing
    
    %There are 3 steps to generate a respiratory motion signal, as shown below    
    % STEP 1: find the coil elements with good representation of respiratory motion
    %         from the late enhancement spokes
    [Coil,Res_Signal_Post]=MC_Step1(ZIP,pars.nSpokesMotionDet);
    
    %STEP 2: Estimate motion signal using PCA from the concatated coil elements
    %Those coil elements were selected in the first step
    [SI,corrm,Res_Signal,ZIP1]=MC_Step2(ZIP,Coil,pars.nSpokesMotionDet,Res_Signal_Post);
    
    %Step 3: You noticed that the signal is not flat, due to the contrast
    %injection. So, now let's estimate the envelop of the signal and substract it
%     Res_Signal=MC_Step3(Res_Signal,ZIP1);
    
    %save the estimated motion signal
%     save Res_Signal.mat Res_Signal
    
    % Free up some space
    clear ZIP ZIP1 
    
    disp('...done.');
end

%% Fourier transform in Z-direction
% So that each 'star' can be reconstructed separately
% This can be done for each coil separately to save memory
% How to treat 'partial fourier' in that direction?
if pars.doFtz
    disp('Performing Fourier Transform in z-direction...');
    [nx, nc, ntviews, nkz] = size(rawdata);
%     if pars.doGpu
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
    if pars.slices ~= 0
        kdata = kdata(:,:,:,pars.slices);
    end
    
    kdata = gather(kdata);
    
    % Free up some memory
    clear rawdata
    
    disp('...done.');
end


%% Coil unstreaking
% Adapted from NYU Demo Demo1_Unstreaking.m
% Should be performed before coil compression
if pars.doUnstreaking
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

if pars.doCoilCompression
    disp('Performing Coil Compression...');
    
    % Permute dimensions:
    % from [nx,nc,ntviews,nz]
    % to   [nx,ntviews,nz,nc]
    kdata = permute(kdata,[1,3,4,2]);
    [nx,ntviews,nz,nc]=size(kdata);

    D=reshape(kdata,nx*nz*ntviews,nc);
    [U,S,V]=svd(D,'econ');
    kdata=single(reshape(D*V(:,1:pars.ncc),nx,ntviews,nz,pars.ncc));
    % save Data/kdata_Unstreaking_CoilCompression.mat kdata
    
    % Free up some memory
    clear D U S V
    
    disp('...done.');
    
    % Permute dimensions back to default [nx,nc,ntviews,nz]
    kdata = permute(kdata, [1,4,2,3]);    
end


%% Perform motion resolved gridding reconstruction
% Based on the demo file 'iterative_multicoil.m' by Tobias Block
if pars.doXdGridding
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
    if pars.slices == 0
        pars.slices = 1:nz;
    end
%     nLinDyn = pars.nLinDyn;
%     if nLinDyn == 0
%         nLinDyn = length(spokes);
%     end
%     nt=floor(ntviews/nLinDyn);              % Number of dynamic time points

    % Generate trajectory
    Traj=Trajectory_GoldenAngle(ntviews,nx,1,[0 0]);
    Traj=Traj*(pars.bas/nx);

    % Calculate number of spokes per respiratory phase & crop data
    nSpokesResp = floor(ntviews/pars.nresp);
    Traj=Traj(:,1:nSpokesResp*pars.nresp);
    kdata=kdata(:,1:nSpokesResp*pars.nresp,:,:);
    Res_Signal=Res_Signal(1:nSpokesResp*pars.nresp);
%     [nx,ntviews,nt,nc] = size(kdata1);
            
    % Sort trajectory according to the respiratory motion signal
    [~,index]=sort(Res_Signal,'descend');
    Traj=Traj(:,index);
    
    % Prepare the NUFFT
    fftSize=[nx,nx];
    %         (k,w,phase,shift,imSize, mode)
    for rPhase = 1:pars.nresp
        FT(rPhase) = NUFFT(Traj(:,(1:nSpokesResp)+(rPhase-1)*nSpokesResp), 1, 1, [0,0], fftSize, 2);
    %             FT = nufft(Traj(:,(1:nSpokesResp)+(rPhase-1)*nSpokesResp), 1, 1, [0,0], fftSize, 2);
    end

    % Initialize output image
    out_gridding    = single(zeros(nx,nx,length(pars.slices),pars.nresp));
    
    tic;
    for sl=pars.slices
        fprintf('Reconstructing slice %d of %d\n', (sl-pars.slices(1)+1), length(pars.slices));
        % Select k-space data of slice to be reconstructed
        kdata1 = kdata(:,:,sl-pars.slices(1)+1,:);
        [nx,ntviews,~,nc]=size(kdata1);         % MCM: seems redundant, but we have a new number of coils since compressing!
                                                % Should be inside slice loop, variables are overwritten again below

        % Coil sensitivity estimation: no algorithm yet for gridding data
        b1 = ones(nx, nx, nc);
        
        % Sort k-space data according to the respiratory motion signal
        kdata1=kdata1(:,index,:,:);
        
        % Calculate density compensation for each resp phase 
        dcf = dcfGridding(nSpokesResp, nx);    %% TO DO: Is this correct?
        
        for rPhase = 1:pars.nresp
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
            out_gridding(:,:,sl-(pars.slices(1)-1),rPhase) = sum(workImage,3);
%             out_gridding(:,:,sl-(pars.slices(1)-1),rPhase) = out_gridding(:,:,sl-(pars.slices(1)-1),rPhase) + abs(workImage.*workImage);
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


%% GROG Reconstruction
% Adapted from Demo4_Reconstruction.m from RACER-GRASP_GROG-GRASP demo
% package (NYU Demo provided by Li Feng)
% This section includes coil sensitivities estimation.

% Outline of the GROG XD-GRASP procedure:
    % - Permute k-space to something handy
    % - Get number of lines per dynamic time point & number of lines per
    % resp phase
    % - Get number of 'Nyquist spokes'. What for again?
    % - Loop over slices:
    % -- Select k-space for this slice (FFT has been done in slice dir)
    % -- Get sizes
    % -- Generate full trajectory (can this be place outside slice loop?)
    % -- Estimate coil sensitivities
    % -- Calculate DCF for 'Nyquist' number of spokes
    % -- Sort data:
    % --- Split trajectory, k-space data and respiratory signal in nt dynamic time points & throw away the last few lines that don't add up to a full dynamic
    % --- Sort trajectory and k-space data according to the respiratory motion signal in each dynamic time point
    % -- Calculate density compensation for each dynamic time point
    % -- Calculate GROG weights
    % -- Prepare k-space data for respiratory and dynamic time point sorted GROG (cartesian interpolation)
    % -- Do GROG; the preinterpolation of radial data into a Cartesian grid
    % -- Perform iterative recon
    % -- End of slice loop
    % - Normalize image
    % - Permute back k-space data
    
if pars.doGrogXdGrasp
    disp('Performing GROG Reconstruction...');
    
    % Permute dimensions of kdata to [nx,ntviews,nz,nc]
    kdata = permute(kdata, [1,3,4,2]);
    
    % Set some initial parameters
    if pars.slices == 0
        pars.slices = 1:nz;
    end
    nLinDyn = pars.nLinDyn;
    if nLinDyn == 0
        nLinDyn = length(spokes);
    end
    nt=floor(ntviews/nLinDyn);              % Number of dynamic time points
    Nqu=min(floor(nx/2 * pi/2), ntviews);   % Nyquist spokes, please refer to the GROG-GRASP paper for more details about this.
                                            % MCM changed from 350
    
    % Initialize output image & recon time measurement
    out_img      = single(zeros(nx,nx,length(pars.slices),pars.nresp));
    tGrog = [];
    for sl=pars.slices % loop through selected slices
        tic;
        fprintf('Reconstructing slice %d of %d\n', (sl-pars.slices(1)+1), length(pars.slices));
        % Select slice to be reconstructed
        kdata1=kdata(:,:,sl-(pars.slices(1)-1),:);
        [nx,ntviews,~,nc]=size(kdata1);         % MCM: seems redundant, but we have a new number of coils since compressing!
                                                % Should be inside slice loop, variables are overwritten again below
                                               
        % Generate trajectory (For GROG only, NOT for general gridding)
        % This needs to stay within the slice loop; it gets modified during
        % sorting of data in respiratory phases
        Traj=Trajectory_GoldenAngle(ntviews,nx,0,[0 1]);        % GROG expects no normalization and reversed spoke traversal
                                                                % MCM changed from Trajectory_GoldenAngle_GROG()
        Traj=Traj*(pars.bas/nx);
        
        %GROG weights calibration
        [Gx,Gy] = GROG.get_Gx_Gy(kdata1,Traj);
        
        %Coil sensitivities estimation
        G = GROG.init(kdata1,Traj,Gx,Gy,0);     % GROG initialization
        kref = GROG.interp(kdata1,G,1,pars.doGpu);   % GROG interpolation; MCM added changed GPU flag (4th param) from 0
        kref=CropImg(kref,pars.bas,pars.bas);
        ref=gather(squeeze(ifft2c_mri(kref)));
        b1=adapt_array_2d(ref);clear ref
        b1=single(b1/max(abs(b1(:))));
        
        % Calculate density compensation for 'Nyquist' number of spokes
        % (see paper, this is the denominator in the Weights calculation)
        % MCM: why does this use the _last_ Nqu spokes of the full
        % trajectory? What's the difference?
        G = GROG.init(kdata1(:,end-Nqu+1:end,:,:),Traj(:,end-Nqu+1:end),Gx,Gy,0);
        D=reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1))]);
        DCF=CropImg(D,pars.bas,pars.bas);
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
        DCF_U=CropImg(DCF_U,pars.bas,pars.bas);
        DCF_U=repmat(DCF_U,[1,1,1,nc]);
        DCF_U(DCF_U==0)=1;
        
        % calculate GROG weights, as described in the paper
        Weighting=DCF./DCF_U;
        Weighting=repmat(Weighting,[1 1 1 1 pars.nresp]);
        
        % Prepare k-space data for respiratory and dynamic time point
        % sorted GROG (cartesian interpolation)
        nline2=floor(ntviews/pars.nresp);                            % Number of lines per respiratory bin
        kdata1=reshape(kdata1,[nx,nline2,pars.nresp,nt,nc]);
        kdata1=reshape(kdata1,[nx,nline2,pars.nresp*nt,nc]);
        Traj=reshape(Traj,[nx,nline2,pars.nresp,nt]);
        Traj=reshape(Traj,[nx,nline2,pars.nresp*nt]);
        
        % Do GROG; the preinterpolation of radial data into a Cartesian grid
        G = GROG.init(kdata1,Traj,Gx,Gy,0);
        kdata1 = GROG.interp(kdata1,G,1,pars.doGpu);                 % MCM added changed GPU flag (4th param) from 0
        kdata1 = CropImg(kdata1,pars.bas,pars.bas);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Iterative reconstruction
        
        kdata1=reshape(kdata1,[pars.bas,pars.bas,pars.nresp,nt,nc]);
        
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
        out_img(:,:,sl-(pars.slices(1)-1),:)=gather(abs(single(recon_cs)));

        tGrog(sl-(pars.slices(1)-1)) = toc;
    end %slice loop
    out_img=out_img/max(out_img(:));
    
%   Permute dimensions back to default [nx,nc,ntviews,nz]
    kdata = permute(kdata, [1,4,2,3]);

    disp('...done.');
end

%% Writing data to Dicom slices
if pars.doDicomWrite
    switch pars.dicomWriteMethod
        case 0 % Even stupider: use NIFTI
            % Create 4D NIFTI structure
            nii = make_nii(out_img, voxSize);          % Uses NIFTI Toolbox by Jimmy Shen
            % Save it
            fName = sprintf('%s_%02d.nii', fileName);
            fName = fullfile(out_path, fName);
            save_nii(nii, fName);                           % Uses NIFTI Toolbox by Jimmy Shen
        case 1 % Write dicom files without any tags: these will be added by postproc module
            if ~exist(out_path, 'dir')
                mkdir(out_path);
            end
            if exist('out_gridding', 'var')
                out1 = uint16(out_gridding*(2^12-1));
                for i=1:size(out1,4) %time points
                    for j=1:size(out1,3) %slices
                        fName = sprintf('grid_slice%d.%d.dcm',j,i);
                        dicomwrite(out1(:,:,j,i), [out_path, '/', fName], 'MultiframeSingleFile', false);
                    end
                end
            end
            
            if exist('out_img', 'var')
                out1 = uint16(out_img*(2^12-1));
                for i=1:size(out1,4) %time points
                    for j=1:size(out1,3) %slices
                        fName = sprintf('slice%d.%d.dcm',j,i);
                        dicomwrite(out1(:,:,j,i), [out_path, '/', fName], 'MultiframeSingleFile', false);
                    end
                end
            end
    end
end



