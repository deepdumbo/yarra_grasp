function pars = initReconPars(reconType)
% Function for initializing reconstruction parameters for advanced
% iterative reconstruction algorithms. Input parameter recontype specifies
% the type of recon the parameter struct pars should be initialized for
% 
% Written by Marnix Maas (Marnix.Maas@radboudumc.nl), May 2018



switch reconType
    case 'full'
        % Set parameters:
        % ...for global settings & flags
        pars.doGpu                   = 0;           % Use GPU where implemented
        pars.slices                  = 0;           % Select which slice(s) to reconstruct. 0 = all.
        pars.doCropImg               = 1;           % Crop images to FOV requested by MR operator
        pars.doFigures               = 0;           % Display output in figures
        
        % ...for selecting which parts of the code to run
        pars.doLoadData              = 1;           % Load data
        pars.doPartialFourier        = 1;			% Perform partial fourier data filling
        pars.doSliceInterp           = 1;			% Perform interpolation in slice direction by filtered zero-filling 
        pars.doMotionDet             = 0;           % Perform motion detection
        pars.doFtz                   = 1;           % Perform FT in z-direction
        pars.doSliceOversampling     = 1;			% Cut out / do not reconstruct slices in the oversampled region
        pars.doUnstreaking           = 0;           % Perform coil unstreaking
        pars.doCoilCompression       = 1;           % Perform coil compression
        pars.doGridding              = 1;           % Perform motion-averaged gridding recon
        pars.doXdGridding            = 0;           % Perform respiratory resolved gridding recon (eg for debugging resp sorting issues)
        pars.doDicomWrite            = 1;           % Write reconstructed data to disk as dicom files
        
        % ...for reading raw data
        pars.doLoadTwixFile          = 0;           % Read file meta data from pre-saved file 'twixData.mat'?
        pars.doSaveTwixFile          = 0;           % Save a file twixData.mat after reading in raw data
		pars.twixFilePath            = '';			% Path to Twix data files
        pars.twixFileName            = '';
        pars.channels                = 0;           % Select which coil channels to read, for testing purposes/reducing memory load. 0 = all.
        pars.spokes                  = 0;           % Select which spokes to read, for testing purposes/reducing memory load. 0 = all.
        pars.echoes                  = 0;           % Select which echoes to read, for testing purposes/reducing memory load. 0 = all.
        pars.partitions              = 0;           % Select which partitions to read, for testing purposes/reducing memory load. 0 = all.
        pars.removeOS                = 0;           % Remove readout oversampling to reduce memory load
        
        % ...for coil unstreaking
        pars.n1=100;                                % Number of spokes to generate artifact free images: 0 = all.
        pars.n2=40;                                 % Number of spokes to generate images with streaking artifacts
        
%         % ...for motion detection
%         if pars.spokes == 0
%             pars.nSpokesMotionDet = pars.spokes;
%         else
%             pars.nSpokesMotionDet = length(pars.spokes);% Number of spokes used for motion detection
%         end
        
        % ...for coil compression
        pars.ncc=8;                                 % Select to how many coil channels the data should be compressed
        
%         % ...for respiratory resolved recon
%         pars.nresp=4;                              % number of respiratory phases
%         if pars.spokes == 0
%             pars.nLinDyn = 0;
%         else
%             pars.nLinDyn=length(pars.spokes);           % number of spokes per dynamic time point. 0 = all (ie 1 time point)
%         end
        pars.bas=512;                               % Number of columns/rows in cartesian-interpolated (ie GROG) k-space data & output image.
        
        % ... for writing data as Dicom 
        pars.dicomWriteMethod = 1;
        pars.out_path = '';
        
    
    case 'low_mem'
        % Version for testing under low-memory conditions (e.g. less
        % spokes)
        % Set parameters:
        % ...for global settings & flags
        pars.doGpu                   = 0;           % Use GPU where implemented
        pars.slices                  = 0;           % Select which slice(s) to reconstruct. 0 = all.
        pars.doCropImg               = 1;           % Crop images to FOV requested by MR operator
        pars.doFigures               = 0;           % Display output in figures
        
        % ...for selecting which parts of the code to run
        pars.doLoadData              = 1;           % Load data
        pars.doPartialFourier        = 1;			% Perform partial fourier data filling
        pars.doSliceInterp           = 1;			% Perform interpolation in slice direction by filtered zero-filling 
        pars.doMotionDet             = 0;           % Perform motion detection
        pars.doFtz                   = 1;           % Perform FT in z-direction
        pars.doSliceOversampling     = 1;			% Cut out / do not reconstruct slices in the oversampled region
        pars.doUnstreaking           = 0;           % Perform coil unstreaking
        pars.doCoilCompression       = 1;           % Perform coil compression
        pars.doGridding              = 1;           % Perform motion-averaged gridding recon
        pars.doXdGridding            = 0;           % Perform respiratory resolved gridding recon (eg for debugging resp sorting issues)
        pars.doDicomWrite            = 1;           % Write reconstructed data to disk as dicom files
        
        % ...for reading raw data
        pars.doLoadTwixFile          = 0;           % Read file meta data from pre-saved file 'twixData.mat'?
        pars.doSaveTwixFile          = 0;           % Save a file twixData.mat after reading in raw data
		pars.twixFilePath            = '';			% Path to Twix data files
        pars.twixFileName            = '';
        pars.channels                = 0;           % Select which coil channels to read, for testing purposes/reducing memory load. 0 = all.
        pars.spokes                  = 1:28;           % Select which spokes to read, for testing purposes/reducing memory load. 0 = all.
        pars.echoes                  = 0;           % Select which echoes to read, for testing purposes/reducing memory load. 0 = all.
        pars.partitions              = 0;           % Select which partitions to read, for testing purposes/reducing memory load. 0 = all.
        pars.removeOS                = 0;           % Remove readout oversampling to reduce memory load
        
        % ...for coil unstreaking
        pars.n1=100;                                % Number of spokes to generate artifact free images: 0 = all.
        pars.n2=40;                                 % Number of spokes to generate images with streaking artifacts
        
%         % ...for motion detection
%         if pars.spokes == 0
%             pars.nSpokesMotionDet = pars.spokes;
%         else
%             pars.nSpokesMotionDet = length(pars.spokes);% Number of spokes used for motion detection
%         end
        
        % ...for coil compression
        pars.ncc=8;                                 % Select to how many coil channels the data should be compressed
        
%         % ...for respiratory resolved recon
%         pars.nresp=4;                              % number of respiratory phases
%         if pars.spokes == 0
%             pars.nLinDyn = 0;
%         else
%             pars.nLinDyn=length(pars.spokes);           % number of spokes per dynamic time point. 0 = all (ie 1 time point)
%         end
        pars.bas=512;                               % Number of columns/rows in cartesian-interpolated (ie GROG) k-space data & output image.
        
        % ... for writing data as Dicom 
        pars.dicomWriteMethod = 1;
        pars.out_path = '';
        
		
	case 'test'
        % Set parameters:
        % ...for global settings & flags
        pars.doGpu                   = 0;           % Use GPU where implemented
        pars.slices                  = 30:31;           % Select which slice(s) to reconstruct. 0 = all.
        pars.doCropImg               = 1;           % Crop images to FOV requested by MR operator
        pars.doFigures               = 1;           % Display output in figures
        
        % ...for selecting which parts of the code to run
        pars.doLoadData              = 1;           % Load data
        pars.doPartialFourier        = 1;			% Perform partial fourier data filling
        pars.doSliceInterp           = 1;			% Perform interpolation in slice direction by filtered zero-filling 
        pars.doMotionDet             = 0;           % Perform motion detection
        pars.doFtz                   = 1;           % Perform FT in z-direction
        pars.doSliceOversampling     = 1;			% Cut out / do not reconstruct slices in the oversampled region
        pars.doUnstreaking           = 0;           % Perform coil unstreaking
        pars.doCoilCompression       = 1;           % Perform coil compression
        pars.doGridding              = 1;           % Perform motion-averaged gridding recon
        pars.doXdGridding            = 0;           % Perform respiratory resolved gridding recon (for debugging resp sorting issues)
        pars.doDicomWrite            = 1;           % Write reconstructed data to disk as dicom files
        
        % ...for reading raw data
        pars.doLoadTwixFile          = 1;           % Read file meta data from pre-saved file 'twixData.mat'?
        pars.doSaveTwixFile          = 0;           % Save a file twixData.mat after reading in raw data
        pars.twixFilePath            = './../twixdata';			% Path to Twix data files
        pars.twixFileName            = 'NanoOes_Ravetestpt01.mat';
        pars.channels                = 0;           % Select which coil channels to read, for testing purposes/reducing memory load. 0 = all.
        pars.spokes                  = 1:110;           % Select which spokes to read, for testing purposes/reducing memory load. 0 = all.
        pars.echoes                  = 1;           % Select which echoes to read, for testing purposes/reducing memory load. 0 = all.
        pars.partitions              = 0;           % Select which partitions to read, for testing purposes/reducing memory load. 0 = all.
        pars.removeOS                = 0;           % Remove readout oversampling to reduce memory load
        
        % ...for coil unstreaking
        pars.n1=100;                                % Number of spokes to generate artifact free images: 0 = all.
        pars.n2=40;                                 % Number of spokes to generate images with streaking artifacts
        
%         % ...for motion detection
%         if pars.spokes == 0
%             pars.nSpokesMotionDet = pars.spokes;
%         else
%             pars.nSpokesMotionDet = length(pars.spokes);% Number of spokes used for motion detection
%         end
        
        % ...for coil compression
        pars.ncc=8;                                 % Select to how many coil channels the data should be compressed
        
%         % ...for respiratory resolved recon
%         pars.nresp=4;                              % number of respiratory phases
%         if pars.spokes == 0
%             pars.nLinDyn = 0;
%         else
%             pars.nLinDyn=length(pars.spokes);           % number of spokes per dynamic time point. 0 = all (ie 1 time point)
%         end
        pars.bas=512;                               % Number of columns/rows in cartesian-interpolated (ie GROG) k-space data & output image.
        
        % ... for writing data as Dicom 
        pars.dicomWriteMethod = 1;
        pars.out_path = '';

end