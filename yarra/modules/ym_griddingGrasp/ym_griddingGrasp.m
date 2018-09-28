function out_gridding = ym_griddingGrasp(in_path, in_file, out_path, temp_path, pars)

% Function for performing a gridding reconstruction of golden-angle radial
% sparse data, both single-echo and multi-echo. No motion compensation is
% performed. Based on NYU demo scripts.
% Follows the convention of Yarra modules for easy integration.
% The input parameter pars can be used to pass optional parameters:
% within Yarra, this appears to be done through a separate config file. We
% will delete it as appropriate.
%
% Written by Marnix Maas (Marnix.Maas@radboudumc.nl), September 2018


%% Parse inputs
clc
if nargin<5
    % TODO: introduce a way of finding out how much memory is available,
    % and choose the recon mode as appropriate
    pars = initReconPars('low_mem');
end


%% Set initial parameters
bp 		= '../../';															% base path for all image recon tools

% ## Include packages in subfolders
addpath(fullfile(bp,'matlabtools/toolboxes/MIRT'));
addpath(fullfile(bp,'matlabtools/toolboxes/MIRT/mri'));                         % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(bp,'matlabtools/toolboxes/MIRT/nufft'));                       % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(bp,'matlabtools/toolboxes/MIRT/utilities'));                   % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(bp,'matlabtools/toolboxes/MIRT/systems'));                     % for partial Fourier (from Jeff Fesslers toolbox)
addpath(fullfile(bp,'matlabtools/toolboxes/mapVBVD/'));                         % for reading TWIX files
addpath(fullfile(bp,'matlabtools/toolboxes/NYU/'));                             % for ...stuff
addpath(fullfile(bp,'matlabtools/toolboxes/NYU/imagescn_R2008a/'));             % for plotting images
% addpath(fullfile(bp,'matlabtools/toolboxes/NYU/MotionDetection/'));             % for motion detection
addpath(fullfile(bp,'matlabtools/operators'));                                  % for operators like NUFFT, Total Variation, GROG, etc
addpath(fullfile(bp,'matlabtools/tools'));                                      % for home-written helper and recon functions



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
    
    % Partial fourier processing in 2D
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
    rawdata = zeroPadPar_Filt(rawdata,nparSliceInterp);
    
    disp('...done.');
end


%% Motion detection
% if pars.doMotionDet
%     disp('Performing motion detection...');
%     doContrastCorr = 0;
%     Res_Signal = motionDetGrasp(rawdata, doContrastCorr, pars.doFigures);
%     
%     disp('...done.');
% end

%% Fourier transform in Z-direction
% So that each 'star' can be reconstructed separately
% This can be done for each coil separately to save memory
if pars.doFtz
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
    
    % Get rid of redundant dimensions (e.g. in case of single-echo)
    kdata = squeeze(kdata);
    
    % delete the slices we don't need anymore
    if pars.slices ~= 0
        kdata = kdata(:,:,:,pars.slices,:);
    else
        pars.slices = 1:size(kdata,4);
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


%% Perform gridding reconstruction
if pars.doGridding
    disp('Performing gridding reconstruction...');
    if ne>1
        out_gridding = reconGriddingMeGrasp(kdata);
    else
        out_gridding = reconGriddingGrasp(kdata);
    end
    
    if pars.doFigures
        for ec = 1:ne
            figure('Name', sprintf('Gridding Solution echo %d', ec)), imagescn(imresize(abs(out_gridding(:,:,:,ec)),   [size(out_gridding,1)*4 size(out_gridding,2)*4],      'bilinear'),[],[],[],3);
        end
    end
    disp('...done.');
end


%% Perform motion resolved gridding reconstruction
% % Based on the demo file 'iterative_multicoil.m' by Tobias Block
% if pars.doXdGridding
%     disp('Performing XD gridding reconstruction...');
%     out_griddingXd = reconGriddingXdGrasp(kdata(:,:,:,slices), Res_Signal, pars.nresp);
%     
%     if pars.doCropImg
%         out_griddingXd = CropImg(out_griddingXd,nImgLin,nImgCol);
%     end
%     
%     if pars.doFigures
%         for ii=1:pars.nresp
%             ttl = sprintf('XD-Gridding: Resp state %d', ii);
%             figure('Name', ttl);
%             imagescn(imresize(abs(out_griddingXd(:,:,:,ii)), [size(out_griddingXd,1)*2 size(out_griddingXd,2)*2], 'bilinear'),[],[],[],3);
%         end
%     end
%     disp('...done.');
% end




%% Writing data to Dicom slices
if pars.doDicomWrite
    switch pars.dicomWriteMethod
        case 1 % Write dicom files without any tags: these will be added by postproc module
            if ~exist(out_path, 'dir')
                mkdir(out_path);
            end
            if exist('out_gridding', 'var')
                out_gridding = uint16(abs(out_gridding)*(2^12-1));
                if ne>1
                for i=1:size(out_gridding,4) %echoes
                    for j=1:size(out_gridding,3) %slices
                        fName = sprintf('grid_slice%d.%d.dcm',j,i);
                        dicomwrite(out_gridding(:,:,j,i), [out_path, '/', fName], 'MultiframeSingleFile', false);
                    end
                end
                else
                    for i=1:size(out_gridding,3) %slices
                        fName = sprintf('grid_slice%d.dcm',i);
                        dicomwrite(out_gridding(:,:,i), [out_path, '/', fName], 'MultiframeSingleFile', false);
                    end
                end
            end
    end %switch
end



