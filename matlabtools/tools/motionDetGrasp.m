%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to perform motion detection on golden-angle radial sparse
% k-space data.
% Selects profiles along kz of the central k-space points (ie kx=ky=0) of
% each spoke and creates z-projection images by iFFT of each profile. Then
% exctracts the respiratory motion signal using PCA methods.
% The option doContrastCorr determines whether a correction for the
% profiles for contrast injection should be performed (1=yes,0=no). Default
% is 1.
% The option doFig (0/1) controls whether profiles and respiratory signal 
% should be plotted in a figure. Default = 0.
%
% USE: Res_Signal = motionDetGrasp(rawdata, doContrastCorr)
% 
% Written by Marnix Maas (Marnix.Maas@radboudumc.nl), 14-sep-2018
% Adapted from Li Feng's NYU demo Demo2_MotionDetection.m


function Res_Signal = motionDetGrasp(rawdata, doContrastCorr, doFig)

if nargin<3
    doFig = 0;
end
if nargin<2
    doContrastCorr = 0;
end

[nx, ntviews, ~, ~, ~]=size(rawdata);       %MCM TODO: handle multi-echo data. 
nSpokesMotionDet = ntviews;                 %MCM TODO: determine this according to useful properties & add input parameters

% Get central k-space profiles along kz
ZIP = squeeze(rawdata(nx/2+2,:,:,:));

% Respiratory motion detection
ZIP=permute(ZIP,[3,2,1]);                                   % MCM changed from [2,1,3]
ZIP=abs(fftshift(ifft(ZIP,400,1),1)); % FFT and interpolation along the kz dimension

%Normalization of each projection in each coil element
ZIP=ProjNorm(ZIP);%Normalization includes temporal smoothing
if doFig
    figure,imagesc(abs(ZIP(:,:,15))),axis image, axis off, colormap(gray),title('Respiratory Motion')
end

%There are 3 steps to generate a respiratory motion signal, as shown below
% STEP 1: find the coil elements with good representation of respiratory motion
%         from the late enhancement spokes
[Coil,Res_Signal_Post]=MC_Step1(ZIP,nSpokesMotionDet,doFig);

%STEP 2: Estimate motion signal using PCA from the concatated coil elements
%Those coil elements were selected in the first step
[SI,corrm,Res_Signal,ZIP1]=MC_Step2(ZIP,Coil,nSpokesMotionDet,Res_Signal_Post,doFig);

%Step 3: You noticed that the signal is not flat, due to the contrast
%injection. So, now let's estimate the envelop of the signal and substract it
if doContrastCorr
    Res_Signal=MC_Step3(Res_Signal,ZIP1,doFig);
end