function [out_grogMwGrasp, kdata, Res_Signal, tGrogMwGrasp] = ym_grogMwGrasp_PancDCE(in_path, in_file, out_path, temp_path, pars)

% Function for performing a respiratory motion-weighted dynamic GRASP 
% reconstruction using GROG pre-interpolation, based on NYU demo scripts by 
% Li Feng.
% Follows the convention of Yarra modules for easy integration.
% The input parameter pars can be used to pass optional parameters. When
% called from Yarra (i.e. without the 'pars' argument), a set of default
% parameters is obtained by calling the accompanying 'initReconPars.m'
% file.
%
% Written by Marnix Maas (Marnix.Maas@radboudumc.nl), May 2018


%% Parse inputs
clc
if nargin<5
    hostname = getHostName();
    switch hostname
        case 'rdbiomr'
            disp('Running on rdbiomr, using full data set');
            pars = initReconPars('full');
        case 'rdcuda'
            disp('Running on rdcuda, using reduced data set to save memory');
            pars = initReconPars('low_mem');
        otherwise
            disp('Running on an unknown machine, using reduced data set to save memory');
            pars = initReconPars('low_mem');
    end
    
end


%% Set initial parameters

% ## Include packages in subfolders
addpath(fullfile(pars.bp,'matlabtools/toolboxes/MIRT'));
addpath(fullfile(pars.bp,'matlabtools/toolboxes/MIRT/mri'));                         % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(pars.bp,'matlabtools/toolboxes/MIRT/nufft'));                       % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(pars.bp,'matlabtools/toolboxes/MIRT/utilities'));                   % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(pars.bp,'matlabtools/toolboxes/MIRT/systems'));                     % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(pars.bp,'matlabtools/toolboxes/mapVBVD/'));                         % for reading TWIX files
addpath(fullfile(pars.bp,'matlabtools/toolboxes/NYU/'));                             % for ...stuff
addpath(fullfile(pars.bp,'matlabtools/toolboxes/NYU/imagescn_R2008a/'));             % for plotting images
addpath(fullfile(pars.bp,'matlabtools/toolboxes/NYU/MotionDetection/'));             % for motion detection
addpath(fullfile(pars.bp,'matlabtools/toolboxes/NYU/Optimizers/'));                  % for (compressed sensing) optimizers
addpath(fullfile(pars.bp,'matlabtools/toolboxes/Nifti toolbox'));                    % for creating & storing output in NIFTI-format 
addpath(fullfile(pars.bp,'matlabtools/operators'));                                  % for operators like NUFFT, Total Variation, GROG, etc
addpath(fullfile(pars.bp,'matlabtools/tools'));                                      % for home-written helper and recon functions



%% Load data
% Adapted from iterativeradial_multicoil.m (NYU Demo)

if pars.doLoadData
    disp('Loading data...');
    
    fileName = fullfile(in_path,in_file);
    [twix, rawdata] = loadDataGoldenAngle(fileName, 'doLoadTwixFile', pars.doLoadTwixFile, ...
                                                    'doSaveTwixFile', pars.doSaveTwixFile, ...
                                                    'twixFilePath', pars.twixFilePath, ...
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
    nparPF = twix.hdr.Meas.Partitions;
    rawdata = pfGoldenAngle(rawdata, nparPF);
    
    disp('...done.');
end


%% Slice interpolation / zero filling
% Zero filling of k-space in partition direction to accommodate slice
% resolution <100% (ie interpolation)
if pars.doSliceInterp
    disp('Doing Slice interpolation...');
    nparSliceInterp = twix.hdr.Meas.NPaftLen;    % Number of partitions for Fourier Transform
    rawdata = zeroPadPar_Filt(rawdata, nparSliceInterp);
    
    disp('...done.');
end


%% Motion detection
if pars.doMotionDet
    disp('Performing motion detection...');
    Res_Signal = motionDetGrasp(rawdata, pars.doContrastCorr, pars.doFigures);
    
    disp('...done.');
end

%% Fourier transform in Z-direction
% So that each 'star' can be reconstructed separately
% This can be done for each coil separately to save memory
if pars.doFtz
    disp('Performing Fourier Transform in z-direction...');
    [nx, nc, ntviews, nkz] = size(rawdata);
    
    % Pre-allocate memory
	kdata = single(zeros(size(rawdata)));

    for c = 1:nc        
        % Inverse Fourier transform along z
        kdata(:,c,:,:) = fftshift(ifft(ifftshift(rawdata(:,c,:,:),4),[],4),4);          % Unsure whether order of fftshift and ifftshift is correct here (or whether that matters)
    end
    
    % Keep only the slices that were requested to be reconstructed
    if pars.slices ~= 0
        kdata = kdata(:,:,:,pars.slices);
    end
    
    % Free up some memory
    clear rawdata
    
    disp('...done.');
end


%% Slice oversampling
% Remove slices that were inside the slice-oversampled region
if pars.doSliceOversampling
    disp('Performing Slice oversampling...');
    
    if pars.slices ~=0
        % Slices have already been removed, do nothing
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
    
    kdata = coilUnstreak(kdata, pars.n1, pars.n2, pars.doFigs);
        
    disp('...done.');
end

%% Coil compression
% Adapted from Demo1_Unstreaking.m (NYU Demo provided by Li Feng)

if pars.doCoilCompression
    disp('Performing Coil Compression...');
    
    kdata = coilCompress(kdata, pars.ncc);
    
    disp('...done.');
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


%% Perform Motion-weighted GROG-GRASP Reconstruction
% Adapted from Demo4_Reconstruction.m from RACER-GRASP_GROG-GRASP demo
% package (NYU Demo provided by Li Feng)
% This section includes coil sensitivities estimation.
if pars.doGrogMwGrasp
    disp('Performing GROG MW-GRASP Reconstruction...');
    
    [out_grogMwGrasp, tGrogMwGrasp] = reconGrogMwDceGrasp(kdata(:,:,:,slices), Res_Signal, pars);
        
    if pars.doCropImg
        if nImgLin<size(out_grogMwGrasp,1) && nImgCol<size(out_grogMwGrasp,2)
            out_grogMwGrasp = CropImg(out_grogMwGrasp,nImgLin,nImgCol);
        end
    end
    
    if pars.doFigures
        nt = floor(ntviews/pars.nLinDyn);
        for ii=1:nt
            ttl = sprintf('GROG-GRASP: dyn time point %d', ii);
            figure('Name', ttl);
            imagescn(imresize(abs(out_grogMwGrasp(:,:,:,ii)), [size(out_grogMwGrasp,1)*2 size(out_grogMwGrasp,2)*2], 'bilinear'),[],[],[],3);
        end
        
%         for ii=nIterOut
%             ttl = sprintf('GROG-GRASP: %d iterations', ii);
%             figure('Name', ttl);
%             imagescn(imresize(abs(out_grogMwGrasp(:,:,:,1,ii)), [size(out_grogMwGrasp,1)*2 size(out_grogMwGrasp,2)*2], 'bilinear'),[],[],[],3);
%         end
    end
    disp('...done.');
    fName = 'out_grogMwGrasp.mat';
    fName = fullfile(out_path, fName);
    save(fName, 'out_grogMwGrasp', 'tGrogMwGrasp');
end



%% Writing data to Dicom slices
if pars.doDicomWrite
    switch pars.dicomWriteMethod
        case 0 % Use NIFTI: can be used as intermediate step before 
               % conversion to Dicom (e.g. using Mevislab)
            % Calculate voxel size
            voxSizeX = twix.hdr.Config.ReadFoV / twix.hdr.Config.NImageCols;
            voxSizeY = twix.hdr.Config.PhaseFoV / twix.hdr.Config.NImageLins;
            slabThickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1, 1}.dThickness;
            voxSizeZ = slabThickness / twix.hdr.Config.NImagePar;
            voxSize = [voxSizeX voxSizeY voxSizeZ]
            % Create 4D NIFTI structure
            out1 = uint16(abs(out_grogMwGrasp)*(2^12-1));
            nii = make_nii(out1, voxSize);          % Uses NIFTI Toolbox by Jimmy Shen
            % Save it
            fName = 'out_img.nii';
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
            
            if exist('out_grogMwGrasp', 'var')
                out1 = uint16(abs(out_grogMwGrasp)*(2^12-1));
                for i=1:size(out1,4) %time points
                    for j=1:size(out1,3) %slices
                        fName = sprintf('slice%d.%d.dcm',j,i);
                        dicomwrite(out1(:,:,j,i), [out_path, '/', fName], 'MultiframeSingleFile', false);
                    end
                end
            end
    end
end



