%% Script for performing GROG-GRASP reconstruction of multi-slice, single time point RAVE acquisition

clear
clc

% ## Include packages in subfolders
addpath('../mapVBVD/');                                             % for reading TWIX files
addpath('../MIRT/irt/mri');                                         % for partial Fourier (from Jeff Fesslers toolbox)
addpath('../MIRT/irt/utilities');                                   % for partial Fourier (from Jeff Fesslers toolbox)
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
doPartialFourier        = 1;            % Perform partial fourier filling of missing k-space data
doMotionDet             = 1;            % Perform motion detection
doFtz                   = 1;            % Perform FT in z-direction
doUnstreaking           = 0;            % Perform coil unstreaking
doCoilCompression       = 1;            % Perform coil compression
doGridding              = 0;            % Perform gridding reconstruction
doGriddingXdGrasp       = 1;            % Perform respiratory motion-resolved griddig reconstruction
doXdGrasp               = 0;            % Perform respiratory motion-resolved iterative recon using XD-GRASP
doGrogXdGrasp           = 0;            % Perform respiratory motion-resolved iterative recon using GROG
doCropImg               = 1;            % Crop images to FOV requested by MR operator
doDicomWrite            = 0;            % Write result as dicom files
 
% ...for reading raw data
loadTwixFile            = true;         % Read file meta data from pre-saved file 'twixData.mat'?
saveTwixFile            = false;        % Save a file twixData.mat after reading in raw data
channels                = 0;            % Select which coil channels to read, for testing purposes/reducing memory load. 0 = all.
spokes                  = 1:300;        % Select which spokes to read, for testing purposes/reducing memory load. 0 = all.
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
doErr = 1;                              % Select if calculation of compression error should be performed

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
%     fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180813_MINAT20/meas_MID00256_FID54020_t1_rave_fs_tra_iso1_2.dat';
    fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180827_MINAT22/meas_MID00166_FID59004_t1_rave_fs_tra_iso1_2_fov260.dat';
    
    
    % Abdomen in vivo data
    % fileName = 'O:\MarMaa\Data\MR3_Prisma\180329_RaveTest_BioMR0009\meas_MID00214_FID04495_t1_rave_fs_tra_iso1_2.dat';
    % fileName = '/mrs_data5/MarMaa/Data/MR3_Prisma/180329_RaveTest_BioMR0009/meas_MID00214_FID04495_t1_rave_fs_tra_iso1_2.dat';    
    
    % Phantom data
    % fileName='O:\MarMaa\Data\MR3_Prisma\180322_RaveTest_Phantom\meas_MID00228_FID01617_RAVE';
    
    
    % ## Read metadata from the Siemens TWIX file
    if loadTwixFile
        load twixDataMinat22_short.mat                               % Header info already read in to speed things up
    else
        twix=mapVBVD(fileName);
        
        % For VD software, the TWIX file may contain adjustment data. We only want
        % to look at the image data for now.
        twix = twix{2};
        if saveTwixFile
            save twixDataMinat22_short.mat twix
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
    
    % Some general parameters for image recon (shouldn't be here...)
    nImgLin = twix.hdr.Config.ImageLines;
    nImgCol = twix.hdr.Config.ImageColumns;    
    
    % Get the k-space data. Data comes in as [samples,channels,spokes,kzlines]
    rawdata = twix.image(:,channels,spokes,kzlines);
    
    disp('...done.');
end

%% Partial Fourier
if doPartialFourier
    disp('Doing Partial Fourier...');
    % Get size of raw data
    [nx, nc, ntviews, npar] = size(rawdata);
    
    % Get full number of k-space lines in partition direction
    nparFull = twix.hdr.Config.NPaftLen;
    
    if npar ~= nparFull
        % Fill missing k-space lines using POCS algorithm:
        % Permute raw data to [nx, npar, ntviews, nc]
        rawdata = permute(rawdata, [1 4 3 2]);
        % Reshape raw data to [nx, npar, ntvies*nc]
        rawdata = reshape(rawdata, [nx npar ntviews*nc]);
        % Loop over all nx-npar 'blades' and perform partial Fourier filling
        kRec = single(zeros(nx, nparFull, ntviews*nc));
        for i=1:ntviews*nc
            [~, kRec(:,:,i)] = ir_mri_partial_fourier_3d(rawdata(:,:,i), [nx nparFull], 'pf_location', [0 0 0], 'niter', 5);
        end
        % Reshape back to [nx, npar, ntviews, nc]
        kRec = reshape(kRec, [nx nparFull ntviews nc]);
        % Permute back to original [nx, nc, ntviews, npar]
        kRec = permute(kRec, [1 4 3 2]);
        
    end
    
    rawdata = kRec;
    clear kRec
    
    % Find partitions to actually reconstruct (ignoring slice oversampled ones)
    startPartition = twix.hdr.Config.LoopStartPar;
    nPartitions = twix.hdr.Config.LoopLengthPar;
    
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
        % Inverse Fourier transform along z
        kdata(:,c,:,:) = fftshift(ifft(ifftshift(rawdata(:,c,:,:),4),[],4),4);          % Unsure whether order of fftshift and ifftshift is correct here (or whether that matters)
    end
    
    % delete the slices we don't need anymore
    if slices ~= 0
        kdata = kdata(:,:,:,slices);
    else
        slices = 1:nparFull;
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
%     Hoe dit werkt:
% - Beschrijf elk datapunt als een vector in coil-ruimte (elke coil
% vertegenwoordigt een as in deze ruimte)
% - Schrijf de data daartoe om naar een Nd (datapunten) x Nc (coils)
% matrix, die de datapunten dus als rijvectoren bevat.
% - Bereken de SVD van deze ruimte.
% - De hieruit volgende matrix V bevat de z.g. 'right singular vectors',
% waarvan de kolommen de eenheidsvectoren bevatten van een nieuwe basis
% voor de datavectoren. Deze basisvectoren zijn lineaire combinaties van de
% oorspronkelijke eenheidsvectoren, die zo zijn gekozen dat de variatie in
% de data langs de betreffende as steeds zo groot mogelijk is, afgezien van
% die langs eerdere assen. (Dus variatie langs as 1 is grootst, langs as 2
% kleiner, etc, maar altijd maximum bereikbare van de overgebleven
% dimensies). In dit geval kunnen deze basisvectoren worden gezien als
% 'gecombineerde spoelelementen'.
% Het product van een data- (=rij-)vector met een kolom van V geeft de
% coordinaat van die vector langs de betreffende V-as. Product van de
% rijvector met matrix V geeft dus dezelfde vector uitgedrukt in de basis
% V, en product D*V geeft alle datavectoren uitgedrukt in de basis V.
% - Waar V de nieuwe basisvectoren bevat op volgorde van data-variatie,
% bevat het SVD-resultaat S de grootte van die variaties. Latere elementen
% van S bevatten lagere waarden, die dus tot uitdrukking brengen dat de
% data langs de bijbehorende V-vector weinig variatie (=informatie) bevat.
% Dit betekent weer dat weglaten van die componenten weinig gevolgen heeft
% voor de daadwerkelijke inhoud van de data
% - Bekijk dus de 'singular values' (uit S), en besluit welke wel en niet
% zullen worden meegenomen
% - Vermenigvuldig D met de kolommen van V die de meeste informatie
% bevatten: data is gecomprimeerd tot een aantal 'gecombineerde
% spoelelementen', dat dus lager is dan het totaal.
    
    
    [kdata_out, S, err] = coilCompress(kdata, ncc, doErr);
    disp('...done.');
end

%% Perform gridding reconstruction
% Based on the demo file 'iterative_multicoil.m' by Tobias Block
if doGridding
    disp('Performing gridding reconstruction...');
    out_gridding = reconGriddingGrasp(kdata);
    
    if doFigures
        figure('Name', 'Gridding Solution'), imagescn(imresize(abs(gridding),   [size(gridding,1)*4 size(gridding,2)*4],      'bilinear'),[],[],[],3);
    end
    disp('...done.');
end
    

%% Perform motion resolved gridding reconstruction
% Based on the demo file 'iterative_multicoil.m' by Tobias Block
if doGriddingXdGrasp
    disp('Performing XD gridding reconstruction...');
    out_griddingXd = reconGriddingXdGrasp(kdata(:,:,:,slices), Res_Signal, nresp);
    
    if doCropImg
        out_griddingXd = CropImg(out_griddingXd,nImgLin,nImgCol);
    end
    
    if doFigures
        for ii=1:nresp
            ttl = sprintf('XD-Gridding: Resp state %d', ii);
            figure('Name', ttl);
            imagescn(imresize(abs(out_griddingXd(:,:,:,ii)), [size(out_griddingXd,1)*2 size(out_griddingXd,2)*2], 'bilinear'),[],[],[],3);
        end
    end
    disp('...done.');
end


%% Perform XD-GRASP Reconstruction
% Adapted from Demo_XDGRASP_NonContrast (NYU Demo provided by Li Feng)

if doXdGrasp
    disp('Performing XD-GRASP reconstruction...');
    
    [out_xdGrasp, out_griddingMcXd, tXdGrasp] = reconXdGrasp(kdata(:,:,:,slices), Res_Signal, nresp);
    if doCropImg
        out_xdGrasp      = CropImg(out_xdGrasp,nImgLin,nImgCol);
        out_griddingMcXd = CropImg(out_griddingMcXd,nImgLin,nImgCol);
    end
    
    if doFigures
        for ii=1:nresp
            ttl = sprintf('MC Gridding: Resp state %d', ii);
            figure('Name', ttl);
            imagescn(imresize(abs(out_griddingMcXd(:,:,:,ii)), [size(out_griddingMcXd,1)*2 size(out_griddingMcXd,2)*2], 'bilinear'),[],[],[],3);
            ttl = sprintf('GRASP: Resp state %d', ii);
            figure('Name', ttl);
            imagescn(imresize(abs(out_xdGrasp(:,:,:,ii)), [size(out_griddingMcXd,1)*2 size(out_griddingMcXd,2)*2], 'bilinear'),[],[],[],3);
        end
    end
    
    disp('...done.');
end


%% Perform GROG XD-GRASP Reconstruction
% Adapted from Demo4_Reconstruction.m from RACER-GRASP_GROG-GRASP demo
% package (NYU Demo provided by Li Feng)
% This section includes coil sensitivities estimation.
if doGrogXdGrasp
    disp('Performing GROG XD-GRASP Reconstruction...');
    
    [out_grogXdGrasp, tGrogXdGrasp] = reconGrogXdGrasp(kdata(:,:,:,slices), Res_Signal, nresp);
    if doCropImg
        out_grogXdGrasp = CropImg(out_grogXdGrasp,nImgLin,nImgCol);
    end
    
    if doFigures
        for ii=1:nresp
            ttl = sprintf('GROG-GRASP: Resp state %d', ii);
            figure('Name', ttl);
            imagescn(imresize(abs(out_grogXdGrasp(:,:,:,ii)), [size(out_grogXdGrasp,1)*2 size(out_grogXdGrasp,2)*2], 'bilinear'),[],[],[],3);
        end
    end
    disp('...done.');
end


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


