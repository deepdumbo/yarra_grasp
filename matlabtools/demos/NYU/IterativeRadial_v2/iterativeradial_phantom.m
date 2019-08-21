%##########################################################################
%
%  Matlab source code for MRI reconstruction method described in:
%  
%  Block KT, Uecker M, Frahm J. 
%  Undersampled radial MRI with multiple coils. Iterative image 
%  reconstruction using a total variation constraint. 
%  Magn Reson Med. 2007 Jun;57(6):1086-98
% 
%  The algorithm in this file corresponds to Figure 3.
%
%  The source code uses the following external packages:
%    - NUFFT toolkit by Jeffrey Fessler 
%      (http://www.eecs.umich.edu/~fessler/)
%    - NUFFT operator by Miki Lustig 
%      (http://www.eecs.berkeley.edu/~mlustig/Software.html)
%    - Siemens TWIX file reader by Philipp Ehses
%    - Poblano Toolbox by Sandia National Laboratories 
%      (http://software.sandia.gov/trac/poblano)
%    - MRI Phantom by Ronald Ouwekerk
%      (http://www.mathworks.com/matlabcentral/fileexchange/1759-mriphantom)
%
%  If you are using this code, please cite the publication listed above.
%
%  Version 25.03.15.
%
%##########################################################################

% ## Clear previous data and variables
clc; 
close all;
clear all;

% ## Set reconstruction parameters

% Note: The penalty weight lambdaPhantom and the stop criterion for the 
% iterations stopTolPhantom needs to be adjusted depending on the value 
% range of the data, and on the undersampling level
param.lambdaPhantom=0.0001;
param.stopTolPhantom=1e-7;
param.iterationsPhantom=150;

% Retrospective undersampling of golden-angle data (0=full data set)
param.undersampleSpokesTo=55;
% Ordering scheme (1=Golden Angle, 0=linear)
param.orderingGoldenAngle=1;
% Downsampling of readout oversampling to reduce matrix size (1=on, 0=off)
param.readoutDownsampling=1;
% Numerical simulation of k-space data (will overwrite scan data)
param.simulateData=1;

param.filename='data/phantom.dat';

% ## Include packages in subfolders
addpath('nufft');                % for the NUFFT
addpath('utilities');            % for the NUFFT
addpath('graph');                % for the NUFFT
addpath('imagescn_R2008a/');     % for plotting images
addpath('read_meas_VBVD/');      % for reading TWIX files
addpath('poblano_toolbox_1.1');  % for numerical optimization
addpath('mriphantom');           % for phantom simulation


% ## Read the data from the Siemens TWIX file
data=mapVBVD(param.filename);

% For VD software, the TWIX file may contain adjustmentdata. We only want
% to look at the image data for now.
data=data.image;

% Read the sequence params from the file
baseresolution=double(data.NCol);
spokes=double(data.NLin);
channels=double(data.NCha);
slices=1;

fprintf('\n');
fprintf('Radial Iterative Reconstruction\n');
fprintf('-------------------------------\n\n');
fprintf('Spokes = %d\n',         spokes        );
fprintf('Channels = %d\n',       channels      );
fprintf('Baseresolution = %d\n', baseresolution);
fprintf('\n');

% Read the k-space data. Data comes in as [samples,channels,spokes,slice]
rawdata=data().image;


% ## Downsampling of readout data (reduces matrix sizes)
if (param.readoutDownsampling==1)
    tempraw=fftshift(fft(ifftshift(rawdata(:,:,:,1),1),[],1),1);
    tempraw=tempraw(end/4+1:3*end/4,:,:,:);
    rawdata=ifftshift(ifft(fftshift(tempraw,1),[],1),1);
    baseresolution=baseresolution/2;
end


% ## Overwrite sequence params for retrospective undersampling
if (param.undersampleSpokesTo~=0)
    spokes=double(param.undersampleSpokesTo);
end


% ## Prepare the density compensation function (DCF)
dcfRow=ones(baseresolution,1);
for i=1:baseresolution
   dcfRow(i)=abs(baseresolution/2-(i-0.5));
end
dcfRow=pi/spokes*dcfRow;
dcf=repmat(dcfRow,1,spokes);


% ## Calculate angles for the radial trajectory
if (param.orderingGoldenAngle==1)
    % For the Golden-Angle scheme
    GA = 111.246117975/180*pi; 
    phi = [pi/2:GA:GA*spokes];
else
    % For the linear ordering mode
    phi = zeros(1,spokes);
    for i=1:spokes
        phi(i)=pi/spokes*(i-1);
        % Flip every second spoke, as done in the sequence
        if (mod(i,2) == 0)
            phi(i)=phi(i)+pi;
        end
    end
end
phi = mod(phi,2*pi);


% ## Calculate the k-space trajectory points 
if (param.readoutDownsampling==0)
    % Without readout downsampling
    rho = linspace(0,baseresolution-1,baseresolution)' - (baseresolution-1)/2;
else
    % With readout downsampling
    rho = linspace(0,baseresolution-1,baseresolution)' - (baseresolution-0.5)/2;
end

% Norm to range -0.5 to 0.5 for the NUFFT
rho = rho/baseresolution;

% Generate vector with k-space coordinates (as complex values kx + i ky)
k = double(rho*exp(-1j*phi));


% ## Simulate k-space data with numerical phantom (will overwrite
% measurement data with numerical Shepp-Logan data)
if (param.simulateData==1)
    fprintf('Simulating data...\n');
    phantomdata=mriphantom(imag(-1*k*baseresolution/(0.7*pi)),real(-1*k*baseresolution/(0.7*pi)));
    
    % Reduce number of channels to 1 for the simulation
    channels=double(1);

    for ic=1:channels
        rawdata(:,ic,1:spokes,1)=phantomdata;
    end
    fprintf('Done.\n');
end


% ## Prepare the NUFFT
fftSize=[baseresolution,baseresolution];
FT = NUFFT(k, 1, 1, [0,0], fftSize, 2);


% ## Perform gridding reconstruction for comparison
fprintf('Calculating gridding reconstruction...');

% Create image to sum the gridded channels (for sum-of-squares combination)
gridding=zeros(baseresolution,baseresolution);
    
% Loop over channels and calculate image for each coil
for channel=1:channels
    % Fetch rawdata for slice and channel and multiply with density
    % compensation function
    workRawData=dcf.*double(squeeze(rawdata(:,channel,1:spokes,1)));

    % Run the NUFFT
    workImage=FT'*(workRawData);

    % Add squared channel image to buffer
    gridding = gridding + abs(workImage.*workImage);
end

% Calculate the root (as final part of the sum-of-squares calculation)
gridding = sqrt(gridding);
fprintf('done.\n');


% ## Iterative reconstruction

% Pass information to the optimizer / cost function
param.nvar=2*baseresolution^2;
param.FT=FT;
param.dcf=dcf;
param.rawdata=rawdata;
param.br=baseresolution;
param.spokes=spokes;
param.channels=channels;

% Global counter for displaying the number of evaluations
global iterationCounter
iterationCounter=0;

% Create figure for displaying the iterations
figure('Name', 'Iterations')

% ## Initialize the optimizer

% Get the rawdata for the cost function
param.y=double(squeeze(rawdata(:,1,1:spokes,1)));
% Start optimizer with empty image
x0=zeros(param.nvar,1);
% Set stopping tolerance
options.StopTol=param.stopTolPhantom;

% First run the optimizer for 10 iterations without penalty terms
param.enablePenalty=0;
options.MaxIters=10;
out = lbfgs(@(x) costfunction_phantom(x,param), x0, options);

% Now enable the penalty terms and run the optimizer for the remaining
% iterations.
param.enablePenalty=1;
options.MaxIters=param.iterationsPhantom-10;
out = lbfgs(@(x) costfunction_phantom(x,param), out.X, options);

% Reshape the result vector into 2D image format
iterative=vec_to_img(out.X,param.br);


% ## Now show the results
fprintf('Finished.');
close;

% Get magnitude and crop center part of images. Show the reconstructed 
% images with 4x bilinear interpolation (looks better)
figure('Name', 'Gridding FFT'), imagescn(imresize(fliplr(abs(fftshift(fft2(fftshift(gridding))))),   [size(gridding,1)*4 size(gridding,2)*4],      'bilinear'),[],[],[],3);
figure('Name', 'Iterative FFT'),imagescn(imresize(fliplr(abs(fftshift(fft2(fftshift(iterative))))),   [size(gridding,1)*4 size(gridding,2)*4],      'bilinear'),[],[],[],3);

figure('Name', 'Gridding Solution'), imagescn(imresize(fliplr(abs(gridding)),   [size(gridding,1)*4 size(gridding,2)*4],      'bilinear'),[],[],[],3);
figure('Name', 'Iterative Solution'),imagescn(imresize(fliplr(abs(iterative)),  [size(iterative,1)*4 size(iterative,2)*4],    'bilinear'),[],[],[],3);
