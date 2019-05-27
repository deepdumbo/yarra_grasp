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
            msg = ('Running on rdbiomr, using full data set');
            logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
            pars = reconPars('full');
        case 'rdcuda'
            msg = ('Running on rdcuda, using reduced data set to save memory');
            logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
            pars = reconPars('low_mem');
        otherwise
            msg = ('Running on an unknown machine, using reduced data set to save memory');
            logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
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
    msg = ('Loading data...');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
    
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
          
    msg = ('...done.');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
end

%% Trajectory
if pars.doCalcTrajectory
    msg = ('Calculating k-space trajectory...');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
    
    % Get data dimensions
    [nx,~,nspokes] = size(rawdata);
    order       = load(fullfile(pars.orderFilePath, pars.orderFileName));
    fldNames    = fieldnames(order);
    fldName     = fldNames{1};
    order       = order.(fldName);
    order(1)    = 1;                        % MCM Only works for the random order provided by UKW!
    traj        = traj_3dUteRandom(nx, nspokes, fullfile(pars.protFilePath, pars.protFileName), order);
    disp('...done.');
    
    % Get file name for SimulationProtocol (for ramp sampling)
    % If not explicitly provided, construct it first
    if isempty(pars.protFileName)
        
        % scanner type
        sysType = strrep(twix.hdr.Dicom.ManufacturersModelName, '_', '');
        
        % gradient mode
        gm = twix.hdr.MeasYaps.sGRADSPEC.ucMode;
        switch gm
            case 8
                gradMode = 'Perf';
            case 4
                gradMode = 'Fast';
            case 2
                gradMode = 'Normal';
        end
        fov = twix.hdr.Config.ReadFoV;
        baseRes = twix.hdr.Config.BaseResolution;
        dt = twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1};
        bw = round(1E9/(dt*baseRes*2));
        
        simProtFileName = sprintf('SimulationProtocol_%s_gm%s_fov%d_r%d_bw%d.txt',sysType, gradMode, fov, baseRes, bw);
        msg = sprintf('No SimulationProtocol file specified, using %s', simProtFileName);
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
    else
        simProtFileName = pars.protFileName;
    end
        
    % Calculate trajectory
%     traj        = traj_3dUteRandom(nx, nspokes, fullfile(pars.protFilePath, pars.protFileName), order);
    traj        = traj_3dUteRandom(nx, nspokes, fullfile(pars.protFilePath, simProtFileName), order);
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
    
    msg = ('...done.');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
end


%% Perform data sorting
% Sort rawdata and k-space trajectory into respiratory motion states
% determined in previous step
if pars.doDataSort
    msg = ('Performing data sorting...');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
    
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
    
    
    msg = ('...done.');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
end


%% Calculate density compensation
if pars.doDcf
    msg = ('Calculating density compensation...');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
    
    % Get image properties (FOV and matrix size)
    FOVx        = twix.hdr.Config.ReadFoV;
    FOVz        = FOVx;
    size_x      = 2 * twix.hdr.Config.BaseResolution;
    size_z      = size_x;
%     msize       = size (traj,2); % Should traj be permuted first? This  number should amount to about 249
        
    for rs = 1:pars.nresp        
        msg = sprintf('... for resp state %d/%d\n', rs, pars.nresp);
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
        
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
    msg = ('...done.');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
end


%% Calculate coil sensitivity maps
if pars.doCoilSensitivities
    msg = ('Calculating coil sensitivity maps...');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
    
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
        msg = ('...calculating low res recon...');
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
        channels = 1:size(rawResp,1);
        reconlowres = lowresrecon(rawResp, trajResp, w, mask, FOVx, FOVz, size_x, size_z, channels , msize);

        % Save low resolution reconstruction
        msg = ('...saving low res recon...');
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
        fNameLowRes = ['reconLowRes_',num2str(rs)];                         % Why do we need to save this? We're not loading it anywhere... can be removed
        save(fullfile(temp_path, fNameLowRes),'reconlowres','-v7.3');
        
        %Determine coil sensitivities map
        msg = ('...calculating coil sensitivities...');
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
        [~, cmps] = CoilSensitivities(reconlowres);
        fNameCmps = ['cmps_breathingstate_',num2str(rs)];
        save(fullfile(temp_path, fNameCmps),'cmps','-v7.3');
        clear reconlowres;    
    end
    
    msg = ('...done.');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
end


%% Perform gridding reconstruction
% try
%     if pars.doGridding
%         msg = ('Performing gridding reconstruction...');
%         logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
%         
%         % Get image properties (FOV and matrix size)
%         imgProperties.FOVx      = twix.hdr.Config.ReadFoV;
%         imgProperties.FOVz      = imgProperties.FOVx;
%         imgProperties.size_x    = 2 * twix.hdr.Config.BaseResolution;
%         imgProperties.size_z    = imgProperties.size_x;
%         out_img = reconGridding3dUte(rawdata, imgProperties, pars);
%         
%         msg = ('...done.');
%         logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
%     else
%         out_img = [];
%     end
% catch eGridding
%     msg = sprintf(2,'Error during gridding reconstruction. Message:\n%s\n',eGridding.message);
%     logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
% end

%% Perform resp motion resolved gridding reconstruction
try
    if pars.doRespResolvedGridding
        msg = ('Performing respiratory motion resolved reconstruction...');
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
        
        % Get image properties (FOV and matrix size)
        imgProperties.FOVx      = twix.hdr.Config.ReadFoV;
        imgProperties.FOVz      = imgProperties.FOVx;
        imgProperties.size_x    = 2 * twix.hdr.Config.BaseResolution;
        imgProperties.size_z    = imgProperties.size_x;
        
        % Reconstruct each motion state
        out_img_respGridding = zeros(imgProperties.size_x, imgProperties.size_x, imgProperties.size_z, pars.nresp);
        for rs = 1:pars.nresp
            msg = (['Reconstructing motion state ', num2str(rs), '/', num2str(pars.nresp)]);
            logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
            spokes = find(breathStates(rs,:)>0);
            out_img_respGridding(:,:,:,rs) = reconGridding3dUte(rawdata(:,:,spokes), traj(spokes,:,:), imgProperties, pars);
        end
        
        msg = ('...done.');
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
    end
catch eRespGridding
    msg = sprintf(2,'Error during respiratory resolved gridding reconstruction. Message:\n%s\n',eRespGridding.message);
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
end


%% Perform resp motion resolved Iterative Sense reconstruction
try
    if pars.doRespResolvedItSense
        msg = ('Performing respiratory motion resolved Iterative SENSE reconstruction...');
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
        
        % Get image properties (FOV and matrix size)
        imgProperties.FOVx      = twix.hdr.Config.ReadFoV;
        imgProperties.FOVz      = imgProperties.FOVx;
        imgProperties.size_x    = 2 * twix.hdr.Config.BaseResolution;
        imgProperties.size_z    = imgProperties.size_x;
        
        
        %load coil sensitivity map
%         msg =  ('loading cmps: ATTENTION HARD CODED TO 4!!!!!!!!!')
%         logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
        fNameCmps = ['cmps_breathingstate_',num2str(respState4CoilSens),'.mat'];
        msg = sprintf('...loading coil sensitivity maps from %s\n', fNameCmps);
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
        load(fullfile(temp_path, fNameCmps), 'cmps');
        
        for rs = 1:pars.nresp
            msg = sprintf('... for resp state %d/%d\n', rs, pars.nresp);
            logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
            
            %load rawdata
            fName = ['rawdata_UTE_breathingstate_' num2str(rs)];
            msg = sprintf('...loading raw data from %s\n', fName);
            logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
            load(fullfile(temp_path, fName), 'rawResp');
            
            %load trajectory
            fName = ['traj_UTE_breathingstate_' num2str(rs)];
            msg = sprintf('...loading trajectory from %s\n', fName);
            logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
            load(fullfile(temp_path, fName), 'trajResp');
            
            %load density compensation
            fName = ['w_breathingstate_' num2str(rs)];
            msg = sprintf('...loading density compensation from %s\n', fName);
            logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
            load(fullfile(temp_path, fName), 'w');
            
            %run reconstruction
            nc = size(rawResp,1);
            msg = sprintf('...starting SENSE recon for motion state %d/%d using %d channels\n', rs, pars.nresp, nc);
            logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
            b = sense(rawResp, trajResp, w, cmps, imgProperties);  %%GM
            msg = sprintf('...saving SENSE recon of motion state %d\n', rs);
            logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
            filename = ['recon_sense_breathingstate_',num2str(rs)];
            save(fullfile(temp_path, filename),'b','-v7.3');            
        end
        
        % Combine motion states
        msg = sprintf('...combining motion states\n');
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
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
        
        msg = ('...done.');
        logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
    end
catch eRespItSense
    msg = sprintf(2,'Error during respiratory resolved iterative sense reconstruction. Message:\n%s\n',eRespItSense.message);
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
end


%% Writing data to Dicom slices
try
if pars.doImageFileWrite
    msg = ('Writing images to file(s)...');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
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
    msg = ('...done.');
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
end
catch eImageFileWrite
    msg = sprintf(2,'Error during image file writing. Message:\n%s',eImageFileWrite.message);
    logRecon(msg, fullfile(temp_path,pars.logFileName), pars.doShowLogMsg);
end



