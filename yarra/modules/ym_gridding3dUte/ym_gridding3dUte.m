function [out_img, out_img_resp] = ym_gridding3dUte(in_path, in_file, out_path, temp_path, pars)

% Function for performing a gridding reconstruction from 3D koosh-ball UTE
% data.
%
% Written by Marnix Maas (Marnix.Maas@radboudumc.nl), January 2019


%% Parse inputs
clc
if nargin<5
    hostname = getHostName();
    switch hostname
        case 'rdbiomr'
            disp('Running on rdbiomr, using full data set');
            pars = reconPars('full');
        case 'rdcuda'
            disp('Running on rdcuda, using reduced data set to save memory');
            pars = reconPars('low_mem');
        otherwise
            disp('Running on an unknown machine, using reduced data set to save memory');
            pars = reconPars('low_mem');
    end
    
end


%% Set initial parameters

% ## Include packages in subfolders
addpath(fullfile(pars.bp,'private/UKW/demos/gridding/mapVBVD - new'));         % for reading TWIX files.
addpath(fullfile(pars.bp,'private/UKW/tools'));                             % for UKW-specific tools

% addpath(fullfile(pars.bp,'matlabtools/toolboxes/NYU/imagescn_R2008a/'));             % for plotting images
addpath(fullfile(pars.bp,'ImageRecon/matlabtools/toolboxes/Nifti_toolbox'));                    % for creating & storing output in NIFTI-format 
addpath(fullfile(pars.bp,'ImageRecon/matlabtools/operators'));                                  % for operators like NUFFT, Total Variation, GROG, etc
addpath(fullfile(pars.bp,'ImageRecon/matlabtools/tools'));                                      % for home-written helper and recon functions

run (fullfile(pars.bp,'ImageRecon/matlabtools/toolboxes/MIRT/setup.m'));

%% Load data

if pars.doLoadData
    disp('Loading data...');
    
    fileName = fullfile(in_path,in_file);
    [twix, rawdata] = loadData(fileName, 'doLoadTwixFile', pars.doLoadTwixFile, ...
                                         'doSaveTwixFile', pars.doSaveTwixFile, ...
                                         'doUnsorted', pars.doUnsorted, ...
                                         'twixFilePath', pars.twixFilePath, ...
                                         'twixFileName', pars.twixFileName, ...
                                         'removeOS', pars.removeOS,...
                                         'channels', pars.channels, ...
                                         'spokes', pars.spokes, ...
                                         'partitions', pars.partitions, ...
                                         'echoes', pars.echoes, ...
                                         'repeats', pars.repeats);    
    % WARNING: The routine for reading 3D Koosh ball data (ie 'unsorted')
    % depends on a different version of mapVBVD. Perhaps the golden angle
    % data can be loaded with this same version, but this should be
    % checked.
          
    disp('...done.');
end

%% Trajectory
if pars.doCalcTrajectory
    disp('Calculating k-space trajectory...');
    % Get data dimensions
    [nx,~,nspokes] = size(rawdata);
    order       = load(fullfile(pars.orderFilePath, pars.orderFileName));
    fldNames    = fieldnames(order);
    fldName     = fldNames{1};
    order       = order.(fldName);
    order(1)    = 1;                        % MCM Only works for the random order provided by UKW!
    traj        = traj_3dUteRandom(nx, nspokes, fullfile(pars.protFilePath, pars.protFileName), order);
    disp('...done.');
end
    

%% Perform motion detection
if pars.doMotionDet
    disp('Performing motion detection...');
    
    breathStates = motionDet3dUte_dc(rawdata, pars);
    
    disp('...done.');
end

%% Coil compression
% To be implemented for 3D Koosh ball?

try
    if pars.doCoilCompression
        disp('Performing Coil Compression...');
        
        rawdata = coilCompress(rawdata, pars.ncc);
        
        disp('...done.');
    end
catch eCoilCompression
    fprintf(2,'Error during coil compression. Message:\n%s',eCoilCompression.message);
end


%% Perform gridding reconstruction
try
    if pars.doGridding
        disp('Performing gridding reconstruction...');
        
        % Get image properties (FOV and matrix size)
        imgProperties.FOVx      = twix.hdr.Config.ReadFoV;
        imgProperties.FOVz      = imgProperties.FOVx;
        imgProperties.size_x    = 2 * twix.hdr.Config.BaseResolution;
        imgProperties.size_z    = imgProperties.size_x;
        out_img = reconGridding3dUte(rawdata, imgProperties, pars);
        
        disp('...done.');
    else
        out_img = [];
    end
catch eGridding
    fprintf(2,'Error during gridding reconstruction. Message:\n%s',eGridding.message);
end

%% Perform resp motion resolved gridding reconstruction
try
    if pars.doRespResolvedRecon
        disp('Performing respiratory motion resolved reconstruction...');
        
        % Get image properties (FOV and matrix size)
        imgProperties.FOVx      = twix.hdr.Config.ReadFoV;
        imgProperties.FOVz      = imgProperties.FOVx;
        imgProperties.size_x    = 2 * twix.hdr.Config.BaseResolution;
        imgProperties.size_z    = imgProperties.size_x;
        
        % Sort raw data & trajectory into resp motion states & save to
        % temp-files
        
        % Reconstruct each motion state
        out_img_resp = zeros(imgProperties.size_x, imgProperties.size_x, imgProperties.size_z, pars.nresp);
        for rs = 1:pars.nresp
            disp(['Reconstructing motion state ', num2str(rs), '/', num2str(pars.nresp)]);
            spokes = find(breathStates(rs,:)>0);
            out_img_resp(:,:,:,rs) = reconItSense3dUte(rawdata(:,:,spokes), traj(spokes,:,:), imgProperties, pars);
        end
        
        disp('...done.');
    end
catch eRespGridding
    fprintf(2,'Error during respiratory resolved gridding reconstruction. Message:\n%s',eRespGridding.message);
end



%% Writing data to Dicom slices
try
if pars.doImageFileWrite
    switch pars.imageFileWriteMethod
        case 0 % Use NIFTI: can be used as intermediate step before 
               % conversion to Dicom (e.g. using Mevislab)
            % Calculate voxel size
            voxSizeX = twix.hdr.Config.ReadFoV / twix.hdr.Config.NImageCols;
            voxSizeY = twix.hdr.Config.PhaseFoV / twix.hdr.Config.NImageLins;
            slabThickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1, 1}.dThickness;
            voxSizeZ = slabThickness / twix.hdr.Config.NPar;
            voxSize = [voxSizeX voxSizeY voxSizeZ];
            if exist('out_img', 'var')
                if ~isempty(out_img)
                    % Create 4D NIFTI structure
                    out1 = uint16(abs(out_img)*(2^12-1));
                    nii = make_nii(out1, voxSize);          % Uses NIFTI Toolbox by Jimmy Shen
                    % Save it
                    fName = 'out_img.nii';
                    fName = fullfile(out_path, fName);
                    save_nii(nii, fName);                           % Uses NIFTI Toolbox by Jimmy Shen
                end
            end
            if exist('out_img_resp', 'var')
                if ~isempty(out_img_resp)
                    % Create 4D NIFTI structure
                    out1 = uint16(abs(out_img_resp)*(2^12-1));
                    nii = make_nii(out1, voxSize);          % Uses NIFTI Toolbox by Jimmy Shen
                    % Save it
                    fName = 'out_img_resp.nii';
                    fName = fullfile(out_path, fName);
                    save_nii(nii, fName);                           % Uses NIFTI Toolbox by Jimmy Shen
                end
            end
        case 1 % Write dicom files without any tags: these will be added by postproc module
            if ~exist(out_path, 'dir')
                mkdir(out_path);
            end
            if exist('out_img', 'var')
                if ~isempty(out_img)
                    out1 = uint16(abs(out_img)*(2^12-1));
                    for i=1:size(out1,4) %time points
                        for j=1:size(out1,3) %slices
                            fName = sprintf('gridding_slice%d.%d.dcm',j,i);
                            dicomwrite(out1(:,:,j,i), [out_path, '/', fName], 'MultiframeSingleFile', false);
                        end
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
catch eImageFileWrite
    fprintf(2,'Error during image file writing. Message:\n%s',eImageFileWrite.message);
end



