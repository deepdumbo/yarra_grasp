classdef reconPars
    %RECONPARS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Parameters:
        % ...for global settings & flags
        bp                      = '~/';         % base path for all image recon tools
        modulePath                              % Path of this module. This parameter is useful, because a call from a YarraServer will result in a different working directory (ie. /yarra/work) compared to a call directly from matlab
        doGpu                   = 0;            % Use GPU where implemented
        slices                  = 0;            % Select which slice(s) to reconstruct. 0 = all.
        doCropImg               = 1;            % Crop images to FOV requested by MR operator
        doFigures               = 0;            % Display output in figures
        logFileName             = 'recon.log';  % Name of log file
        doShowLogMsg            = 1;            % Show log messages on screen
        
        % ...for selecting which parts of the code to run
        doLoadData              = 1;            % Load data
        doCalcTrajectory        = 1;            % Calculate k-space trajectory
        doMotionDet             = 1;            % Perform motion detection
        doDataSort              = 1;            % Perform data sorting
        doDcf                   = 1;            % Calculate density compensation
        doCoilSensitivities     = 1;            % Calculate coil sensitivity maps
        doGridding              = 0;            % Perform gridding recon
        doRespResolvedGridding  = 0;            % Perform respiratory motion resolved Gridding recon
        doRespResolvedItSense   = 1;            % Perform respiratory motion resolved Iterative Sense recon
        doImageFileWrite        = 1;            % Write reconstructed image data to disk
        
        % ...for reading raw data
        doLoadTwixFile          = 0;            % Read file meta data from pre-saved file 'twixData.mat'?
        doSaveTwixFile          = 0;            % Save a file twixData.mat after reading in raw data
        doUnsorted              = 1;            % Read data in 'unsorted' way: needed for handling 3D Koosh ball UKW-style
        twixFilePath                            % Path to Twix data files
        twixFileName            = '';
        channels                = 0;            % Select which coil channels to read, for testing purposes/reducing memory load. 0 = all.
        spokes                  = 0;            % Select which spokes to read, for testing purposes/reducing memory load. 0 = all.
        partitions              = 0;            % Select which partitions to read, for testing purposes/reducing memory load. 0 = all.
        echoes                  = 0;            % Select which echoes to read, for testing purposes/reducing memory load. 0 = all.
        repeats                 = 0;            % Select which repeats to read, for testing purposes/reducing memory load. 0 = all.
        removeOS                = 0;            % Remove readout oversampling to reduce memory load
        ignoreSeg               = 1;            % Merge 'segments' together to reduce memory load
        
        % ...for calculating the trajectory
        protFilePath
        protFileName
        orderFilePath
        orderFileName           = 'new_random_order.mat';
        
        % ...for motion detection
%         filesource              = '/mnt/resfilsp03/data/nis_home/mrs/marnix/projects/software/ImageReconstruction/yarra/modules/ym_gridding3dUte';  %% CHECK
        transientphase          = 2000; % CHECK
        dcPoints                = 2:6;          % Select which points in each spoke should be used to determine the DC signal
        sgfFrameLength          = 801;         % Frame length of the DC signal smoothing filter
        sgfOrder                = 3;           % Polynomial order of the DC signal smoothing filter
        blcOrder                = 8;            % Polynomial order of DC signal baseline correction
        ncx                     = 3;            % Number of coil elements to combine for resp motion detection
        
        % ...for respiratory binning
        nresp                   = 4;  %%%%%%%%%%%%%%%% CHECK
        MinPe4Reco              = 0; %%%%%%%%%%%%%%%% CHECK: was 45000
        
        % ...for data selection and sorting
        initialPoint            = 0;            % First point to be used within each spoke: used to eliminate artifacts due to digital filter. 0=autoCalculate based on digFilterDelay
        digFilterDelay          = 13;           % Duration in us of the initial ADC period during which the digital filter causes instable output. Determines the number of initial points in each FID to be discarded. 
                                                % Empirically found to be ~13us for PrismaFit.
        doClearRaw              = 1;            % Clear rawdata after data sorting (to save memory)
        
        % ...for density compensation calculation
        nIterPcg                = 10;           % number of iterations for PCG algorithm
        
        % ...for coil sensitivity estimation
        
        % ...for static gridding recon
        doFermiFilter          = 0;            % Apply a Fermi filter to the reconstruction operator
        
        % ...for iterative SENSE recon
        doFlip                  = 1;            % Flip output image in 1st dimension
        
        % ... for writing output data to file(s)
        imageFileWriteMethod   = 1;             % 0 = NIFTI, 1 = DICOM
%         out_path = '';
    end
    
    methods
        function pars = reconPars(reconType)
            % Constructor
            pars.modulePath              = fullfile(pars.bp, 'ReconTools/yarra/modules/ym_itSense3dUte/');
            pars.twixFilePath            = fullfile(pars.bp, 'twixdata');
            pars.protFilePath            = fullfile(pars.modulePath, 'SimulationProtocols');
            pars.protFileName            = '';
            pars.orderFilePath           = fullfile(pars.modulePath, './');
            pars.orderFileName           = 'new_random_order.mat';
            
            if nargin == 1
                switch reconType
                    case 'full'
                    case 'low_mem'
                        pars.spokes                  = 0;           % Select which spokes to read, for testing purposes/reducing memory load. 0 = all.
                        pars.echoes                  = 1;           % Select which echoes to read, for testing purposes/reducing memory load. 0 = all.
                        pars.repeats                 = 1:30;       % Select which repeats to read, for testing purposes/reducing memory load. 0 = all.
                        pars.removeOS                = 0;           % Remove readout oversampling to reduce memory load
                        
                    case 'test'
                        % Set parameters:
                        % ...for global settings & flags
%                         pars.bp                      = '~/';        % base path for all image recon tools
                        pars.doGpu                   = 0;           % Use GPU where implemented
                        pars.slices                  = 0;           % Select which slice(s) to reconstruct. 0 = all.
                        pars.doCropImg               = 1;           % Crop images to FOV requested by MR operator
                        pars.doFigures               = 0;           % Display output in figures
                        
                        % ...for selecting which parts of the code to run
                        pars.doLoadData              = 1;           % Load data
                        pars.doCalcTrajectory        = 1;           % Calculate k-space trajectory
                        pars.doMotionDet             = 1;           % Perform motion detection
                        pars.doDataSort              = 1;           % Perform data sorting
                        pars.doDcf                   = 1;           % Calculate density compensation
                        pars.doCoilSensitivities     = 1;           % Calculate coil sensitivity maps
                        pars.doGridding              = 0;           % Perform gridding recon
                        pars.doRespResolvedGridding  = 0;           % Perform respiratory motion resolved Gridding recon
                        pars.doRespResolvedItSense   = 1;           % Perform respiratory motion resolved Iterative Sense recon
                        pars.doImageFileWrite        = 1;           % Write reconstructed image data to disk
                        
                        % ...for reading raw data
                        pars.doLoadTwixFile          = 0;           % Read file meta data from pre-saved file 'twixData.mat'?
                        pars.doSaveTwixFile          = 0;           % Save a file twixData.mat after reading in raw data
%                         pars.twixFilePath            = fullfile(pars.bp,'twixdata');			% Path to Twix data files
                        pars.twixFilePath            = '/nis_home/mrs/marnix/projects/software/ImageReconstruction/twixdata';			% Path to Twix data files
                        pars.twixFileName            = 'MINAT46_3dUte.mat';
                        pars.channels                = 1:2;           % Select which coil channels to read, for testing purposes/reducing memory load. 0 = all.
                        pars.spokes                  = 0;           % Select which spokes to read, for testing purposes/reducing memory load. 0 = all.
                        pars.partitions              = 0;           % Select which partitions to read, for testing purposes/reducing memory load. 0 = all.
                        pars.echoes                  = 0;           % Select which echoes to read, for testing purposes/reducing memory load. 0 = all.
                        pars.repeats                 = 1:30;       % Select which repeats to read, for testing purposes/reducing memory load. 0 = all.
                        pars.removeOS                = 0;           % Remove readout oversampling to reduce memory load
                        
                        % ...for calculating the trajectory
                        pars.protFilePath            = '~/ReconTools/yarra/modules/ym_itSense3dUte/SimulationProtocols';
                        pars.protFileName            = '';          % Leave empty for automatic generation based on acquisition parameters
                        pars.orderFilePath           = '~/ReconTools/yarra/modules/ym_itSense3dUte';
                        
                        % ...for motion detection
%                         pars.filesource              = '/mnt/resfilsp03/data/nis_home/mrs/marnix/projects/software/ImageReconstruction/yarra/modules/ym_gridding3dUte';  %% CHECK
                        pars.transientphase          = 200;        % Number of initial spokes to discard for motion detection, was 2000
                        pars.dcPoints                = 2:6;         % Select which points in each spoke should be used to determine the DC signal
                        pars.sgfFrameLength          = 801;         % Frame length of the DC signal smoothing filter
                        pars.sgfOrder                = 3;           % Polynomial order of the DC signal smoothing filter
                        pars.blcOrder                = 8;           % Polynomial order of DC signal baseline correction
                        pars.ncx                     = 3;           % Number of coil elements to combine for resp motion detection
                        
                        % ...for respiratory binning
                        pars.nresp                   = 2;           % Number of respiratory states
                        pars.MinPe4Reco              = 0;       % Minimum number of spokes per respiratory state, was 12500
%                         pars.MinPe4Reco              = 50000;       % Minimum number of spokes per respiratory state, was 12500
%                         pars.MinPe4Reco              = 2000*max(pars.repeats)/pars.nresp;
                        
                        % ...for data selection and sorting
                        pars.initialPoint           = 0;            % First point to be used within each spoke: used to eliminate artifacts due to digital filter. 0=autoCalculate based on digFilterDelay
                        pars.digFilterDelay          = 13;           % Duration in us of the initial ADC period during which the digital filter causes instable output. Determines the number of initial points in each FID to be discarded. 
                                                % Empirically found to be ~13us for PrismaFit.
%                         pars.initialPoint           = 8;            % First point to be used within each spoke: used to eliminate artifacts due to digital filter
%                         pars.doClearRaw              = 1;            % Clear rawdata after data sorting (to save memory)
%                         
                        % ...for density compensation calculation
                        pars.nIterPcg               = 2;           % number of iterations for PCG algorithm, default 10
                        
                        % ...for coil sensitivity estimation
                        
                        % ...for static gridding recon
                        pars.doFermiFilter          = 0;            % Apply a Fermi filter to the reconstruction operator

                        % ...for iterative SENSE recon
                        
                        % ... for writing output data to file(s)
                        pars.imageFileWriteMethod   = 1;
%                         pars.out_path = '';
                        
                end % switch
            end %if nargin==1
        end % constructor
        
    
    end % methods
end % classdef

