%%%%%%%%v%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This demo demonstrates the automatic detection of a respiratory motion signal 
%from a stack-of-stars DCE dataset


%HUUH WAAROM DOEN WE DIT DAN NOG KLOJOS
%Note that the ZIP.mat file was generated using the centers of k-space
%after performing a FFT along the slice dimension.

%Due to space limits, the full 3D k-space were not uploaded.

% Li Feng, NYU, 12/18/2017
% Li.Feng@nyumc.org

clear
clc
addpath Data
addpath nufft_files
addpath utils
close all
load Data/ZIP.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Suppose the size of your stack-of-stars k-space is
% [nx nviews nz nc]=size(kdata);

% HERE
% nx = readout points
% nviews = the total number of spokes
% nz = the number of slices (NOTE, this dimension is kz, not z)
% nc = the number of coil elements

% then, you can get your ZIP with
%ivom: nx/2+2 is a single point! Selects the kspace center.
%ivom: hoogste signaal nemen ipv ~midden?
% ZIP = kdata(nx/2+2,:,:,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Respiratory motion detection
%ivom: after below permutation, kdata/ZIP has format 
%ivom: kdata(nz, spokes, ncoils), ie, in the case of ZIP:
%ivom: 1144 spokes (spokes) with stack-of-stars 'depth' of 38 (nz) probed
%ivom: with 20 separate coils. 
%ivom: HERE ZIP [nx=1, nviews, nz, nc]
ZIP=permute(ZIP,[2,1,3]);
%ivom: HERE ZIP [nviews, nx, 

%ivom: Zero-pad up to 400 elements (!). FFT along k_z dimension w/ only the
%ivom: the centers of k-space included in ZIP (see above and Feng, 2018).

%ivom: From help fft ...
%ivom:  "fft(X,N) is the N-point fft, padded with zeros if X has less
%ivom:  than N points and truncated if it has more."
%ivom: Big zero-pad! From 38 to 400 elements. 
ZIP1=abs(fftshift(ifft(ZIP,400,1),1)); % FFT and interpolation along the kz dimension

%Normalization of each projection in each coil element
ZIP=ProjNorm(ZIP1);%Normalization includes temporal smoothing
%figure,imagesc(abs(ZIP(:,:,15))),axis image, axis off, colormap(gray),title('Respiratory Motion')

%There are size3 steps to generate a respiratory motion signal, as shown below

n1=800; 
% the last 800 spokes were used as the late enhancement phase for motion detection
% as described in the paper

% STEP 1: find the coil elements with good representation of respiratory motion
%         from the late enhancement spokes
%ivom: late enhancement phase: contrast is constant, little variation due
%ivom: to injection.
disp(size(ZIP));
[Coil,Res_Signal_Post]=MC_Step1(ZIP,n1);
fprintf('Coil...');
disp(Coil);
%STEP 2: Estimate motion signal using PCA from the concatated coil elements
%Those coil elements were selected in the first step
[SI,corrm,Res_Signal,ZIP1]=MC_Step2(ZIP,Coil,n1,Res_Signal_Post);

%Step 3: You noticed that the signal is not flat, due to the contrast
%injection. So, now let's estimate the envelop of the signal and substract it
Res_Signal=MC_Step3(Res_Signal,ZIP1);

%save the estimated motion signal
%save Data/Res_Signal.mat Res_Signal


