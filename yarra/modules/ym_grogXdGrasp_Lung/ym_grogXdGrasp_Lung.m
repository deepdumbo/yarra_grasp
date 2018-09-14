function [out_grogXdGrasp, out_griddingXd] = ym_grogXdGrasp_Lung(in_path, in_file, out_path, temp_path, pars)

% Function for performing a respiratory resolved XD-GRASP reconstruction
% using GROG pre-interpolation, based on NYU demo scripts by Li Feng.
% Follows the convention of Yarra modules for easy integration.
% The input parameter pars can be used to pass optional parameters:
% within Yarra, this appears to be done through a separate config file. We
% will delete it as appropriate.
%
% Written by Marnix Maas (Marnix.Maas@radboudumc.nl), May 2018


%% Parse inputs
clc
if nargin<5
    pars = initReconPars('grogXdGrasp_lung');
end


%% Set initial parameters
bp 		= '../../../';															% base path for all image recon tools

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
addpath(fullfile(bp,'matlabtools/demos/NYU/RACER-GRASP_GROG-GRASP/'));          % for plotting images
addpath(fullfile(bp,'matlabtools/demos/NYU/RACER-GRASP_GROG-GRASP/utils'));     % for GROG-specific tools
% addpath('../NYU Demos/Demo_RACER-GRASP_GROG-GRASP/nufft_files');            	% for GROG-specific tools
addpath(fullfile(bp,'matlabtools/operators'));                                  % for operators like NUFFT, Total Variation, GROG, etc
addpath(fullfile(bp,'matlabtools/tools'));                                      % for home-written helper and recon functions



%% Load data
% Adapted from iterativeradial_multicoil.m (NYU Demo)

if pars.doLoadData
    disp('Loading data...');
    
    fileName = fullfile(in_path,in_file);
    [twix, rawdata] = loadDataGoldenAngle(fileName, 'doLoadTwixFile', pars.doLoadTwixFile, ...
                                                    'doSaveTwixFile', pars.doSaveTwixFile, ...
                                                    'twixFileName', pars.twixFileName, ...
                                                    'removeOS', pars.removeOS,...
                                                    'channels', pars.channels, ...
                                                    'spokes', pars.spokes, ...
                                                    'partitions', pars.partitions, ...
                                                    'echoes', pars.echoes);    
                                                
    % Some general parameters for image recon (shouldn't be here...)
    nImgLin = twix.hdr.Config.ImageLines;
    nImgCol = twix.hdr.Config.ImageColumns;
    
    disp('...done.');
end


%% Partial Fourier
if pars.doPartialFourier
    disp('Doing Partial Fourier...');
    
    % Multi-echo partial fourier processing in 2D
    rawdata = pfGoldenAngleMe(rawdata, twix);
    
    disp('...done.');
end


%% Slice interpolation / zero filling
% Zero filling of k-space in partition direction to accommodate slice
% resolution <100% (ie interpolation)
if pars.doSliceInterp
    disp('Doing Slice interpolation...');
    nparSliceInterp = twix.hdr.Meas.NPaftLen;    % Number of partitions for Fourier Transform
    rawdata = zeroPadPar_Filt(rawdata,nparSliceInterp);
    
    disp('...done.');
end


%% Motion detection
if pars.doMotionDet
    disp('Performing motion detection...');
    
    Res_Signal = motionDetGrasp(rawdata, doContrastCorr, pars.doFigures);
    
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


%% Slice oversampling
% Remove slices that were inside the slice-oversampled region
if pars.doSliceOversampling
    disp('Performing Slice oversampling...');
    
    if pars.slices ~=0
        % Slices already removed, do nothing
        slices = 1:size(kdata,4);
    else
        % Find partitions to actually reconstruct (ignoring slice oversampled ones)
        startPartition = twix.hdr.Config.LoopStartPar;
        nparImg = twix.hdr.Config.LoopLengthPar;
        slices = startPartition:startPartition+nparImg-1;
        kdata = kdata(:,:,:,slices);
        slices = 1:nparImg;
    end
    
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
    disp('Performing XD gridding reconstruction...');
    out_griddingXd = reconGriddingXdGrasp(kdata(:,:,:,slices), Res_Signal, pars.nresp);
    
    if pars.doCropImg
        out_griddingXd = CropImg(out_griddingXd,nImgLin,nImgCol);
    end
    
    if pars.doFigures
        for ii=1:pars.nresp
            ttl = sprintf('XD-Gridding: Resp state %d', ii);
            figure('Name', ttl);
            imagescn(imresize(abs(out_griddingXd(:,:,:,ii)), [size(out_griddingXd,1)*2 size(out_griddingXd,2)*2], 'bilinear'),[],[],[],3);
        end
    end
    disp('...done.');
end


%% Perform GROG XD-GRASP Reconstruction
% Adapted from Demo4_Reconstruction.m from RACER-GRASP_GROG-GRASP demo
% package (NYU Demo provided by Li Feng)
% This section includes coil sensitivities estimation.
if pars.doGrogXdGrasp
    disp('Performing GROG XD-GRASP Reconstruction...');
    
    [out_grogXdGrasp, tGrogXdGrasp] = reconGrogXdGrasp(kdata(:,:,:,slices), Res_Signal, pars.nresp);
    if pars.doCropImg
        out_grogXdGrasp = CropImg(out_grogXdGrasp,nImgLin,nImgCol);
    end
    
    if pars.doFigures
        for ii=1:pars.nresp
            ttl = sprintf('GROG-GRASP: Resp state %d', ii);
            figure('Name', ttl);
            imagescn(imresize(abs(out_grogXdGrasp(:,:,:,ii)), [size(out_grogXdGrasp,1)*2 size(out_grogXdGrasp,2)*2], 'bilinear'),[],[],[],3);
        end
    end
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
                out1 = uint16(abs(out_griddingXd)*(2^12-1));
                for i=1:size(out1,4) %time points
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
    end
end



