pars.bp = '~/project/ReconTools/'

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
% addpath(fullfile(pars.bp,'matlabtools/toolboxes/NYU/Optimizers/'));                  % for (compressed sensing) optimizers
% addpath(fullfile(pars.bp,'matlabtools/demos/NYU/RACER-GRASP_GROG-GRASP/'));          % for plotting images
% addpath(fullfile(pars.bp,'matlabtools/demos/NYU/RACER-GRASP_GROG-GRASP/utils'));     % for GROG-specific tools
% addpath('../NYU Demos/Demo_RACER-GRASP_GROG-GRASP/nufft_files');            	% for GROG-specific tools
addpath(fullfile(pars.bp,'matlabtools/toolboxes/Nifti_toolbox'));               % for creating & storing output in NIFTI-format 
addpath(fullfile(pars.bp,'matlabtools/operators'));                                  % for operators like NUFFT, Total Variation, GROG, etc
addpath(fullfile(pars.bp,'matlabtools/tools'));                                      % for home-written helper and recon functions

work_path = '~/project/temp/basic_grasp';


meas_file = 'meas_MID00588_FID193705_t2s_RAVE_mGRE3_head.dat';
meas_file2 = 'meas_MID00075_FID192066_t2s_RAVE_mGRE3.dat';


temp_path = '~/project/temp/trash';
output_path = '~/project/temp/basic_grasp'; 

mode_file = 'modes/GRASP_basic.mode';

data = yarra_GRASP0611(work_path, meas_file, output_path,...
    temp_path, mode_file)