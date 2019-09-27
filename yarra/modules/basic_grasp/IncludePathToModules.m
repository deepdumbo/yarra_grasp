function IncludePathToModules()
    % Set paths that contain the necessary modules and data
    pars.bp = '~/project/ReconTools/';
    addpath(fullfile(pars.bp,'matlabtools/toolboxes/MIRT'));
    addpath(fullfile(pars.bp,'matlabtools/toolboxes/MIRT/mri'));                         % for partial Fourier (from Jeff Fesslers toolbox)
    addpath(fullfile(pars.bp,'matlabtools/toolboxes/MIRT/nufft'));                       % for partial Fourier (from Jeff Fesslers toolbox)
    addpath(fullfile(pars.bp,'matlabtools/toolboxes/MIRT/utilities'));                   % for partial Fourier (from Jeff Fesslers toolbox)
    addpath(fullfile(pars.bp,'matlabtools/toolboxes/MIRT/systems'));                     % for partial Fourier (from Jeff Fesslers toolbox)
    addpath(fullfile(pars.bp,'matlabtools/toolboxes/mapVBVD/'));                         % for reading TWIX files
    addpath(fullfile(pars.bp,'matlabtools/toolboxes/NYU/'));                             % for ...stuff
    addpath(fullfile(pars.bp,'matlabtools/toolboxes/NYU/imagescn_R2008a/'));             % for plotting images
    addpath(fullfile(pars.bp,'matlabtools/toolboxes/NYU/MotionDetection/'));             % for motion detection
    addpath(fullfile(pars.bp,'matlabtools/toolboxes/Nifti_toolbox'));                    % for creating & storing output in NIFTI-format 
    addpath(fullfile(pars.bp,'matlabtools/operators'));                                  % for operators like NUFFT, Total Variation, GROG, etc
    addpath(fullfile(pars.bp,'matlabtools/tools'));    
end