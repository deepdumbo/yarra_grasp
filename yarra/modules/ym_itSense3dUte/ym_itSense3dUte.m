function [out_img_respSense, out_img_respGridding] = ym_itSense3dUte(in_path, in_file, out_path, temp_path, pars)

% Function for performing an iterative Sense reconstruction from 3D 
% koosh-ball UTE data.
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
addpath(fullfile(pars.bp,'private/UKW/demos/gridding/mapVBVD - new'));          % for reading TWIX files.
addpath(fullfile(pars.bp,'private/UKW/tools'));                                 % for UKW-specific tools
addpath(fullfile(pars.bp,'private/UKW/operators'));                             % for UKW recon operators

% addpath(fullfile(pars.bp,'toolboxes/NYU/imagescn_R2008a/'));             % for plotting images
addpath(fullfile(pars.bp,'ImageRecon/matlabtools/toolboxes/Nifti_toolbox'));                    % for creating & storing output in NIFTI-format 
% addpath(fullfile(pars.bp,'ImageRecon/matlabtools/operators'));                                  % for operators like NUFFT, Total Variation, GROG, etc
addpath(fullfile(pars.bp,'ImageRecon/matlabtools/tools'));                      % for home-written helper and recon functions

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
    
    % The matrix breathStates contains info on which spokes belong to which
    % motion states. Row indices represent motion states, column indices
    % represent spokes: if breathStates(i,j)>0, then spoke j has been
    % assigned to motion state i. Each spoke can be assigned to multiple
    % motion states, i.e. overlap can exist.
    breathStates = motionDet3dUte_dc(rawdata, pars);
    fNameBS   = 'UTE_breathingstates';
    save(fullfile(temp_path, fNameBS),'breathStates');
    
    disp('...done.');
end


%% Perform data sorting
% Sort rawdata and k-space trajectory into respiratory motion states
% determined in previous step
if pars.doDataSort
    disp('Performing data sorting...');
    
    channels = pars.channels;
    if channels == 0
        channels = 1:twix.image.NCha;
    end
    
    for rs =1:pars.nresp
        % Select spokes
        spokes      = breathStates(rs,:)>0;
        
        % Sort & save trajectory
        trajTmp     = traj(spokes,pars.initialPoint:end,:);                 % Select spokes to use & cut off first few points of each spoke
        trajTmp     = permute(trajTmp,[3 2 1]);                             % Permute dimensions to [3(xyz), nx, spokes]
%         trajTmp     = trajTmp(:,pars.initialPoint:end,:);                   % Cut off first few points of each spoke
        trajTmp     = trajTmp(:,:);                                         % Concatenate all spokes along 2nd dimension
        trajResp    = permute(trajTmp,[2 1]);                               % Permute to [points*spokes, 3]
        
        fNameTraj   = ['traj_UTE_breathingstate_',num2str(rs)];
        save(fullfile(temp_path, fNameTraj),'trajResp','-v7.3');
        
        % Sort & save raw data
        rawResp     = rawdata(pars.initialPoint:end,channels,spokes);       % Select spokes & channels to use, and cut first few points off of each spoke
        rawResp     = permute(rawResp,[2 1 3]);                             % Permute dimensions to [channels, nx, spokes]
        rawResp     = rawResp (:,:);                                        % Concatenate all spokes along 2nd dimension
        
        fNameRaw    = ['rawdata_UTE_breathingstate_',num2str(rs)];
        save(fullfile(temp_path, fNameRaw),'rawResp','-v7.3');
        
    end
    
%     clear trajTmp trajectory  
    
    if pars.doClearRaw
        clear rawResp rawdata
    end
    
    
    disp('...done.');
end


%% Calculate density compensation
if pars.doDcf
    disp('Calculating density compensation...');
    
    % Get image properties (FOV and matrix size)
    FOVx        = twix.hdr.Config.ReadFoV;
    FOVz        = FOVx;
    size_x      = 2 * twix.hdr.Config.BaseResolution;
    size_z      = size_x;
%     msize       = size (traj,2); % Should traj be permuted first? This  number should amount to about 249
        
    for rs = 1:pars.nresp        
        fprintf('... for resp state %d/%d\n', rs, pars.nresp);        
        mask = [];  %MCM: does not seem to be  used in densitycompensation() anyway...
        
        % Select which spokes to use and reshape trajectory as needed
        spokes  = breathStates(rs,:)>0;
        trajTmp = traj(spokes,:,:);
        trajTmp = permute(trajTmp,[3 2 1]);
        trajTmp = trajTmp(:,pars.initialPoint:end,:);
        msize   = size (trajTmp,2);
        trajTmp = trajTmp(:,:);
        trajTmp = permute(trajTmp,[2 1]);
        
        % Calculate density compensation
        % Note: densitycompensation() expects the k-space trajectory in the
        % form Np x 3, where Np is the number of k-space *points* sampled,
        % i.e. all points in the trajectory are passed as one long list,
        % instead of separated out in different spokes.
        w = densitycompensation(trajTmp, mask, FOVx, FOVz, size_x, size_z, msize, pars.nIterPcg);
        
        fNameW = ['w_breathingstate_',num2str(rs)];
        save(fullfile(temp_path, fNameW),'w','-v7.3');
        
        clear kspace4Recon;
        clear w;
        
        
    end
    disp('...done.');
end


%% Calculate coil sensitivity maps
if pars.doCoilSensitivities
    disp('Calculating coil sensitivity maps...');
    
%     respStates = pars.respState4CoilSens;
    respState4CoilSens = ceil(pars.nresp/2);
    for rs = respState4CoilSens
        %load rawdata
        load(fullfile(temp_path, ['rawdata_UTE_breathingstate_' num2str(rs)]), 'rawResp');
        
        %load trajectory
        load(fullfile(temp_path, ['traj_UTE_breathingstate_' num2str(rs)]), 'trajResp');
        
        %load density compensation
        load(fullfile(temp_path, ['w_breathingstate_' num2str(rs)]), 'w');
        
        %Calculate low resolution reconstruction
        disp('...calculating low res recon...');
        channels = 1:size(rawResp,1);
        reconlowres = lowresrecon(rawResp, trajResp, w, mask, FOVx, FOVz, size_x, size_z, channels , msize);

        % Save low resolution reconstruction
        disp('...saving low res recon...');
        fNameLowRes = ['reconLowRes_',num2str(rs)];                         % Why do we need to save this? We're not loading it anywhere... can be removed
        save(fullfile(temp_path, fNameLowRes),'reconlowres','-v7.3');
        
        %Determine coil sensitivities map
        disp('...calculating coil sensitivities...');
        [~, cmps] = CoilSensitivities(reconlowres);
        fNameCmps = ['cmps_breathingstate_',num2str(rs)];
        save(fullfile(temp_path, fNameCmps),'cmps','-v7.3');
        clear reconlowres;    
    end
    
    disp('...done.');
end


%% Perform gridding reconstruction
% try
%     if pars.doGridding
%         disp('Performing gridding reconstruction...');
%         
%         % Get image properties (FOV and matrix size)
%         imgProperties.FOVx      = twix.hdr.Config.ReadFoV;
%         imgProperties.FOVz      = imgProperties.FOVx;
%         imgProperties.size_x    = 2 * twix.hdr.Config.BaseResolution;
%         imgProperties.size_z    = imgProperties.size_x;
%         out_img = reconGridding3dUte(rawdata, imgProperties, pars);
%         
%         disp('...done.');
%     else
%         out_img = [];
%     end
% catch eGridding
%     fprintf(2,'Error during gridding reconstruction. Message:\n%s\n',eGridding.message);
% end

%% Perform resp motion resolved gridding reconstruction
try
    if pars.doRespResolvedGridding
        disp('Performing respiratory motion resolved reconstruction...');
        
        % Get image properties (FOV and matrix size)
        imgProperties.FOVx      = twix.hdr.Config.ReadFoV;
        imgProperties.FOVz      = imgProperties.FOVx;
        imgProperties.size_x    = 2 * twix.hdr.Config.BaseResolution;
        imgProperties.size_z    = imgProperties.size_x;
        
        % Reconstruct each motion state
        out_img_respGridding = zeros(imgProperties.size_x, imgProperties.size_x, imgProperties.size_z, pars.nresp);
        for rs = 1:pars.nresp
            disp(['Reconstructing motion state ', num2str(rs), '/', num2str(pars.nresp)]);
            spokes = find(breathStates(rs,:)>0);
            out_img_respGridding(:,:,:,rs) = reconGridding3dUte(rawdata(:,:,spokes), traj(spokes,:,:), imgProperties, pars);
        end
        
        disp('...done.');
    end
catch eRespGridding
    fprintf(2,'Error during respiratory resolved gridding reconstruction. Message:\n%s\n',eRespGridding.message);
end


%% Perform resp motion resolved Iterative Sense reconstruction
try
    if pars.doRespResolvedItSense
        disp('Performing respiratory motion resolved Iterative SENSE reconstruction...');
        
        % Get image properties (FOV and matrix size)
        imgProperties.FOVx      = twix.hdr.Config.ReadFoV;
        imgProperties.FOVz      = imgProperties.FOVx;
        imgProperties.size_x    = 2 * twix.hdr.Config.BaseResolution;
        imgProperties.size_z    = imgProperties.size_x;
        
        
        %load coil sensitivity map
%         disp ('loading cmps: ATTENTION HARD CODED TO 4!!!!!!!!!')
        fNameCmps = ['cmps_breathingstate_',num2str(respState4CoilSens),'.mat'];
        fprintf('...loading coil sensitivity maps from %s\n', fNameCmps);
        load(fullfile(temp_path, fNameCmps), 'cmps');
        
        for rs = 1:pars.nresp
            fprintf('... for resp state %d/%d\n', rs, pars.nresp);  
            
            %load rawdata
            fName = ['rawdata_UTE_breathingstate_' num2str(rs)];
            fprintf('...loading raw data from %s\n', fName);
            load(fullfile(temp_path, fName), 'rawResp');
            
            %load trajectory
            fName = ['traj_UTE_breathingstate_' num2str(rs)];
            fprintf('...loading trajectory from %s\n', fName);
            load(fullfile(temp_path, fName), 'trajResp');
            
            %load density compensation
            fName = ['w_breathingstate_' num2str(rs)];
            fprintf('...loading density compensation from %s\n', fName);
            load(fullfile(temp_path, fName), 'w');
            
            %run reconstruction
            nc = size(rawResp,1);
            fprintf('...starting SENSE recon for motion state %d/%d using %d channels\n', rs, pars.nresp, nc);
            b = sense(rawResp, trajResp, w, cmps, imgProperties);  %%GM
            fprintf('...saving SENSE recon of motion state %d\n', rs);
            filename = ['recon_sense_breathingstate_',num2str(rs)];
            save(fullfile(temp_path, filename),'b','-v7.3');            
        end
        
        % Combine motion states
        fprintf('...combining motion states\n');
        size_x = imgProperties.size_x;
        out_img_respSense=zeros(size_x,size_x,size_x,pars.nresp);       % Was 8,256,256,256
        for rs=1:pars.nresp
            load(fullfile(temp_path, ['recon_sense_breathingstate_',num2str(rs),'.mat']), 'b');
            out_img_respSense(:,:,:,rs)=b;
        end
        
        % Crop image (remove oversampling)
        
        % Normalize image
        out_img_respSense = out_img_respSense./max(abs(out_img_respSense(:)));
        
        % Save image as .mat-file        
        filename = 'recon_sense_all';
        save(fullfile(out_path, filename),'out_img_respSense','-v7.3');
        
        disp('...done.');
    end
catch eRespItSense
    fprintf(2,'Error during respiratory resolved iterative sense reconstruction. Message:\n%s\n',eRespItSense.message);
end


%% Writing data to Dicom slices
try
if pars.doImageFileWrite
    disp('Writing images to file(s)...');
    switch pars.imageFileWriteMethod
        case 0 % Use NIFTI: can be used as intermediate step before 
               % conversion to Dicom (e.g. using Mevislab)
            % Calculate voxel size
            voxSizeX = twix.hdr.Config.ReadFoV / twix.hdr.Config.NImageCols;
            voxSizeY = twix.hdr.Config.PhaseFoV / twix.hdr.Config.NImageLins;
            slabThickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1, 1}.dThickness;
            voxSizeZ = slabThickness / twix.hdr.Config.NPar;
            voxSize = [voxSizeX voxSizeY voxSizeZ];
            if exist('out_img_respSense', 'var')
                if ~isempty(out_img_respSense)
                    % Create 4D NIFTI structure
                    out1 = uint16(abs(out_img_respSense)*(2^12-1));
                    nii = make_nii(out1, voxSize);          % Uses NIFTI Toolbox by Jimmy Shen
                    % Save it
                    fName = 'out_img_respSense.nii';
                    fName = fullfile(out_path, fName);
                    save_nii(nii, fName);                           % Uses NIFTI Toolbox by Jimmy Shen
                end
            end
            if exist('out_img_respGridding', 'var')
                if ~isempty(out_img_respGridding)
                    % Create 4D NIFTI structure
                    out1 = uint16(abs(out_img_respGridding)*(2^12-1));
                    nii = make_nii(out1, voxSize);          % Uses NIFTI Toolbox by Jimmy Shen
                    % Save it
                    fName = 'out_img_respGridding.nii';
                    fName = fullfile(out_path, fName);
                    save_nii(nii, fName);                           % Uses NIFTI Toolbox by Jimmy Shen
                end
            end
        case 1 % Write dicom files without any tags: these will be added by postproc module
            if ~exist(out_path, 'dir')
                mkdir(out_path);
            end
            if exist('out_img_respSense', 'var')
                if ~isempty(out_img_respSense)
                    out1 = uint16(abs(out_img_respSense)*(2^12-1));
                    for i=1:size(out1,4) %time points
                        for j=1:size(out1,3) %slices
                            fName = sprintf('sense_slice%d.%d.dcm',j,i);
                            dicomwrite(out1(:,:,j,i), [out_path, '/', fName], 'MultiframeSingleFile', false);
                        end
                    end
                end
            end
    end
    disp('...done.');
end
catch eImageFileWrite
    fprintf(2,'Error during image file writing. Message:\n%s',eImageFileWrite.message);
end



