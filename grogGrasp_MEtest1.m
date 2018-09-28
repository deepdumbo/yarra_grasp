%% Script for performing GROG-GRASP reconstruction of multi-slice, multi-echo single time point RAVE acquisition

clear
clc

%% Set initial parameters
bp 		= './';															% base path for all image recon tools

% ## Include packages in subfolders
addpath(fullfile(bp,'matlabtools/toolboxes/MIRT'));
addpath(fullfile(bp,'matlabtools/toolboxes/MIRT/mri'));                         % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(bp,'matlabtools/toolboxes/MIRT/nufft'));                       % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(bp,'matlabtools/toolboxes/MIRT/utilities'));                   % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(bp,'matlabtools/toolboxes/MIRT/systems'));                     % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(bp,'matlabtools/toolboxes/mapVBVD/'));                         % for reading TWIX files
addpath(fullfile(bp,'matlabtools/toolboxes/NYU/'));                             % for ...stuff
addpath(fullfile(bp,'matlabtools/toolboxes/NYU/imagescn_R2008a/'));             % for plotting images
addpath(fullfile(bp,'matlabtools/toolboxes/NYU/MotionDetection/'));             % for motion detection
addpath(fullfile(bp,'matlabtools/operators'));                                  % for operators like NUFFT, Total Variation, GROG, etc
addpath(fullfile(bp,'matlabtools/tools'));    


% Set parameters:
% ...for global settings & flags
doFigures               = 1;            % Display output in figures
doGpu                   = 0;            % Use GPU where implemented
slices                  = 30:35;            % Select which slice(s) to reconstruct. 0 = all.
doCropImg               = 1;            % Crop images to FOV requested by MR operator

% ...for selecting which parts of the code to run
doLoadData              = 1;            % Load data
doPartialFourier        = 1;            % Perform partial fourier filling of missing k-space data
doSliceInterp           = 1;            % Perform slice interpolation (when slice resolution was <100%)
doMotionDet             = 0;            % Perform motion detection
doFtz                   = 1;            % Perform FT in z-direction
doSliceOversampling     = 1;            % Remove slice-oversampled slices from data
doUnstreaking           = 0;            % Perform coil unstreaking
doCoilCompression       = 1;            % Perform coil compression
doGridding              = 1;            % Perform gridding reconstruction
doGriddingXdGrasp       = 0;            % Perform respiratory motion-resolved griddig reconstruction
doXdGrasp               = 0;            % Perform respiratory motion-resolved iterative recon using XD-GRASP
doGrogGrasp             = 0;            % Perform CS reconstruction using GROG
doGrogXdGrasp           = 0;            % Perform respiratory motion-resolved iterative recon using GROG
doDicomWrite            = 1;            % Write result as dicom files
 
% ...for reading raw data
doLoadTwixFile          = 1;            % Read file meta data from pre-saved file 'twixData.mat'?
doSaveTwixFile          = 0;            % Save a file twixData.mat after reading in raw data
channels                = 0;            % Select which coil channels to read, for testing purposes/reducing memory load. 0 = all.
spokes                  = 1:110;            % Select which spokes to read, for testing purposes/reducing memory load. 0 = all.
partitions              = 0;            % Select which kz line(s) to read from the data. 0 = all.
echoes                  = 0;            % Select which echoes to read, for testing purposes/reducing memory load. 0 = all.
removeOS                = 0;            % Remove readout oversampling to reduce memory load

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


%% Load data
% Adapted from iterativeradial_multicoil.m (NYU Demo)

if doLoadData
    disp('Loading data...');
    % Oesophagus in vivo data
%     fileName = 'O:\MarMaa\Data\MR3_Prisma\180809_RAVEpatients_slokdarm/NanoOes_Ravetestpt01_t2s_RAVE_mGRE3.dat';
%     fileName = 'O:\MarMaa\Data\MR3_Prisma\180809_RAVEpatients_slokdarm/NanoOes_Ravetestpt02_t2s_RAVE_mGRE3.dat';
    fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180809_RAVEpatients_slokdarm/NanoOes_Ravetestpt01_t2s_RAVE_mGRE3.dat';
%     fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180809_RAVEpatients_slokdarm/NanoOes_Ravetestpt02_t2s_RAVE_mGRE3.dat';
    
    % Phantom data
    % fileName='O:\MarMaa\Data\MR3_Prisma\180322_RaveTest_Phantom\meas_MID00228_FID01617_RAVE';
    % fileName='/mrs_data5/MarMaa/Data/MR3_Prisma/180906_RavePhantom/meas_MID00027_FID63140_RAVE_YGXL.dat';
    
    % Name of Twix file to read and/or write
    twixFilePath    = './twixdata/';
    twixFileName    = 'NanoOes_Ravetestpt01.mat';
    [twix, rawdata] = loadDataGoldenAngle(fileName, 'doLoadTwixFile', doLoadTwixFile, ...
                                                    'doSaveTwixFile', doSaveTwixFile, ...
                                                    'twixFilePath', twixFilePath, ...
                                                    'twixFileName', twixFileName, ...
                                                    'removeOS', removeOS,...
                                                    'channels', channels, ...
                                                    'spokes', spokes, ...
                                                    'partitions', partitions, ...
                                                    'echoes', echoes);    
    disp('...done.');
end

%% Partial Fourier
if doPartialFourier
    disp('Doing Partial Fourier...');
    % Ik moet:
    % - Gevraagde slices opzoeken
    % - Opzoeken hoeveel partities daadwerkelijk data moeten bevatten (in
    %   tegenstelling tot 0en). Dit is het aantal lijnen waar naartoe moet
    %   worden ge-PFt.
    %   Te vinden in:
    % --- twix.hdr.Meas.Partitions
    % - Partial Fourier uitvoeren
    % - Partities na oversampling opzoeken. Dit is het aantal partities
    %   waarover straks de FT in z-richting wordt gedaan. Dit kan
    %   bijvoorbeeld via:
    % --- twix.hdr.Config/Meas.NPaftLen, 
    % --- twix.hdr.Config/Meas.RawPar,
    % --- of via twix.hdr.Dicom.flSliceOS)
    % - Zero padden...
    
    
    % Multi-echo partial fourier processing in 2D
    rawdata = pfGoldenAngleMe(rawdata, twix);
    
    
    
    disp('...done.');
end


%% Slice interpolation / zero filling
% Zero filling of k-space in partition direction to accommodate slice
% resolution <100% (ie interpolation)
if doSliceInterp
    disp('Doing Slice interpolation...');
    % This may not use the correct parameters to determine to how many
    % lines to zero-fill...
    nparSliceInterp = twix.hdr.Meas.NPaftLen;    % Number of partitions for Fourier Transform
%     nparSliceInterp = 2*(twix.hdr.Config.LoopStartPar-1) + twix.hdr.Config.LoopLengthPar;    
    rawdata = zeroPadPar_Filt(rawdata,nparSliceInterp);
    
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
    
    [nx, nc, ntviews, npar, ne] = size(rawdata);
    
    % Select k-space data to use for motion detection: only use 1st echo
    ZIP = squeeze(rawdata(nx/2+2,:,:,:,1));
    
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
if doFtz
    disp('Performing Fourier Transform in z-direction...');
    [nx, nc, ntviews, npar, ne] = size(rawdata);
%     if doGpu
%         kdata = gpuArray(single(zeros(size(rawdata))));
%     else
        kdata = single(zeros(size(rawdata)));
%     end
    for e = 1:ne
        for c = 1:nc
            % Inverse Fourier transform along z
            kdata(:,c,:,:,e) = fftshift(ifft(ifftshift(rawdata(:,c,:,:,e),4),[],4),4);          % Unsure whether order of fftshift and ifftshift is correct here (or whether that matters)
        end
    end
    
    % delete the slices we don't need anymore
    if slices ~= 0
        kdata = kdata(:,:,:,slices,:);
    else
        slices = 1:size(kdata,4);
    end
    
    kdata = gather(kdata);
    
    % Free up some memory
    clear rawdata
    
    disp('...done.');
end

%% Slice oversampling
% Remove slices that were inside the slice-oversampled region
if doSliceOversampling
    disp('Performing Slice oversampling...');
    % This may not use the correct parameters to determine to how many
    % lines to remove oversampled lines...
    
    if slices ~=0
        % Slices already removed, do nothing
        slices = 1:size(kdata,4);
    else
        % Find partitions to actually reconstruct (ignoring slice oversampled ones)
        startPartition = twix.hdr.Config.LoopStartPar;
        nparImg = twix.hdr.Config.LoopLengthPar;
        slices = startPartition:startPartition+nparImg-1;
        kdata = kdata(:,:,:,slices,:);
        slices = 1:nparImg;
    end
    
    
%     % Find partitions to actually reconstruct (ignoring slice oversampled ones)
%     startPartition = twix.hdr.Config.LoopStartPar;
%     nparImg = twix.hdr.Config.LoopLengthPar;
% %     nparImg = twix.hdr.Config.NImagePar;
%     slices = startPartition:startPartition+nparImg-1;
%     kdata = kdata(:,:,:,slices,:);
    
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
    % from [nx,nc,ntviews,nz,ne]
    % to   [nx,ntviews,nz,nc,ne]
    kdata = permute(kdata,[1,3,4,2,5]);
    [nx,ntviews,nz,nc,ne]=size(kdata);
    
    % What is best? 
    % - do this for each coil separately
    % - or for all data combined?
    % - coil selection based on 1st echo only; but apply to all echoes

    D=reshape(kdata,nx*ntviews*nz,nc,ne);
    kdata1 = zeros(nx,ntviews,nz,ncc,ne);
    for ec = 1:ne
        [U,S,V]=svd(D(:,:,ec),'econ');
        kdata1(:,:,:,:,ec)=single(reshape(D(:,:,ec)*V(:,1:ncc),nx,ntviews,nz,ncc));
    end
    kdata = kdata1;
    % Free up some memory
    clear D U S V kdata1
    
    disp('...done.');
    
    % Permute dimensions back to default [nx,nc,ntviews,nz,ne]
    kdata = permute(kdata, [1,4,2,3,5]);    
end

%% Perform gridding reconstruction
if doGridding
    disp('Performing gridding reconstruction...');
    out_griddingMe = reconGriddingMeGrasp(kdata);
    
    if doFigures
        for ec = 1:ne
            figure('Name', sprintf('Gridding Solution echo %d', ec)), imagescn(imresize(abs(out_griddingMe(:,:,:,ec)),   [size(out_griddingMe,1)*4 size(out_griddingMe,2)*4],      'bilinear'),[],[],[],3);
        end
    end
    disp('...done.');
end
    

%% Perform motion resolved gridding reconstruction
% Based on the demo file 'iterative_multicoil.m' by Tobias Block
% if doGriddingXdGrasp
%     disp('Performing XD gridding reconstruction...');
%     out_griddingXd = reconGriddingXdGrasp(kdata(:,:,:,slices), Res_Signal, nresp);
%     
%     if doCropImg
%         out_griddingXd = CropImg(out_griddingXd,nImgLin,nImgCol);
%     end
%     
%     if doFigures
%         for ii=1:nresp
%             ttl = sprintf('XD-Gridding: Resp state %d', ii);
%             figure('Name', ttl);
%             imagescn(imresize(abs(out_griddingXd(:,:,:,ii)), [size(out_griddingXd,1)*2 size(out_griddingXd,2)*2], 'bilinear'),[],[],[],3);
%         end
%     end
%     disp('...done.');
% end


%% Perform GRASP Reconstruction



%% Perform GROG-GRASP Reconstruction
% Perform single-volume (ie no dynamic time points, no respiratory
% binning), multi-echo reconstruction using GROG interpolation and
% compressed sensing. This eincludes coil sensitivities estimation.
if doGrogGrasp
    disp('Performing GROG-GRASP Reconstruction...');
    
    [out_grogMeGrasp, tGrogMeGrasp] = reconGrogMeGrasp(kdata(:,:,:,slices,:));
    if doCropImg
        out_grogMeGrasp = CropImg(out_grogMeGrasp,nImgLin,nImgCol);
    end
    
    if doFigures
        for ii=1:ne
            ttl = sprintf('GROG-GRASP: Echo %d', ii);
            figure('Name', ttl);
            imagescn(imresize(abs(out_grogMeGrasp(:,:,:,ii)), [size(out_grogMeGrasp,1)*4 size(out_grogMeGrasp,2)*4], 'bilinear'),[],[],[],3);
        end
    end
    disp('...done.');
end


%% Perform XD-GRASP Reconstruction
% Adapted from Demo_XDGRASP_NonContrast (NYU Demo provided by Li Feng)

% if doXdGrasp
%     disp('Performing XD-GRASP reconstruction...');
%     
%     [out_xdGrasp, out_griddingMcXd, tXdGrasp] = reconXdGrasp(kdata(:,:,:,slices), Res_Signal, nresp);
%     if doCropImg
%         out_xdGrasp      = CropImg(out_xdGrasp,nImgLin,nImgCol);
%         out_griddingMcXd = CropImg(out_griddingMcXd,nImgLin,nImgCol);
%     end
%     
%     if doFigures
%         for ii=1:nresp
%             ttl = sprintf('MC Gridding: Resp state %d', ii);
%             figure('Name', ttl);
%             imagescn(imresize(abs(out_griddingMcXd(:,:,:,ii)), [size(out_griddingMcXd,1)*2 size(out_griddingMcXd,2)*2], 'bilinear'),[],[],[],3);
%             ttl = sprintf('GRASP: Resp state %d', ii);
%             figure('Name', ttl);
%             imagescn(imresize(abs(out_xdGrasp(:,:,:,ii)), [size(out_griddingMcXd,1)*2 size(out_griddingMcXd,2)*2], 'bilinear'),[],[],[],3);
%         end
%     end
%     
%     disp('...done.');
% end


%% Perform GROG XD-GRASP Reconstruction
% Adapted from Demo4_Reconstruction.m from RACER-GRASP_GROG-GRASP demo
% package (NYU Demo provided by Li Feng)
% This section includes coil sensitivities estimation.
% if doGrogXdGrasp
%     disp('Performing GROG XD-GRASP Reconstruction...');
%     
%     [out_grogXdGrasp, tGrogXdGrasp] = reconGrogXdGrasp(kdata(:,:,:,slices), Res_Signal, nresp);
%     if doCropImg
%         out_grogXdGrasp = CropImg(out_grogXdGrasp,nImgLin,nImgCol);
%     end
%     
%     if doFigures
%         for ii=1:nresp
%             ttl = sprintf('GROG-GRASP: Resp state %d', ii);
%             figure('Name', ttl);
%             imagescn(imresize(abs(out_grogXdGrasp(:,:,:,ii)), [size(out_grogXdGrasp,1)*2 size(out_grogXdGrasp,2)*2], 'bilinear'),[],[],[],3);
%         end
%     end
%     disp('...done.');
% end


%% Writing data to Dicom slices
if doDicomWrite
%     grogGrasp_xd = uint16(grogGrasp_xd./max(grogGrasp_xd(:))*(2^12 - 1));
    
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
            if ~exist(out_path, 'dir')
                mkdir(out_path);
            end
            if exist('out_griddingMe', 'var')
                out1 = uint16(abs(out_griddingMe)*(2^12-1));
                for i=1:size(out1,4) %echoes
                    for j=1:size(out1,3) %slices
                        fName = sprintf('grid_slice%d.%d.dcm',j,i);
                        dicomwrite(out1(:,:,j,i), [out_path, '/', fName], 'MultiframeSingleFile', false);
                    end
                end
            end
            
            if exist('out_grogXdGrasp', 'var')
                out1 = uint16(abs(out_grogXdGrasp)*(2^12-1));
                for i=1:size(out1,4) %time points
                    for j=1:size(out1,3) %slices
                        fName = sprintf('slice%d.%d.dcm',j,i);
                        dicomwrite(out1(:,:,j,i), [out_path, '/', fName], 'MultiframeSingleFile', false);
                    end
                end
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


