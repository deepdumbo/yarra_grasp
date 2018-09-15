%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to perform Coil Unstreaking as described in Li Feng's paper
% [...]
%
% Marnix Maas, RadboudUMC, 15-9-2018
% Adapted from Li Feng's Demo1_Unstreaking.m
% 
% Parameters:
% - kdata:        k-space data *after Fourier Transform in z-direction!*
% - n1:
% - n2:
% - doFigs:       flag to control whether result should be displayed.
% Default = 0
%
% USE: kdata = coilUnstreak(kdata, pars)

function kdata = coilUnstreak(kdata, n1, n2, doFigs)

if nargin<4
    doFigs = 0;
end

% Permute dimensions to prepare for Coil Unstreaking
kdata = permute(kdata,[1,3,4,2]);

[nx,ntviews,nz,nc]=size(kdata);
% Note that the third dimension is z, NOT kz. (A FFT was performed already)

% Generate trajectory (For GROG only, NOT for general gridding)
Traj=Trajectory_GoldenAngle_ME(ntviews, nx, ne, [], 'normalize',1,'flip',[0 0]);

% Calculate density compensation
DensityComp = dcfGridding(ntviews, nx);

Ref=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n1,nx));                 % artifact-free image
Img=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n2,nx)*n1/n2);           % image with streaks

%As described in the paper, the Diff image is calculated as the 2x FOV
Diff=abs(Ref-Img);

%The Ref and Img are then cropped to the 1x FOV
Ref=Ref(nx/4+1:end-nx/4,nx/4+1:end-nx/4,:,:);
Img=Img(nx/4+1:end-nx/4,nx/4+1:end-nx/4,:,:);

data=sos(Ref,4);
data(:,:,:,2)=sos(Img,4);

%calculating the streak ratio
clear StreakRatio
for ii=1:nc
    StreakRatio(ii,1)=norm(col(Diff(:,:,:,ii)))/norm(col(Ref(:,:,:,ii)));
end
StreakRatio=StreakRatio/min(StreakRatio);
if doFigs
    figure,plot(StreakRatio) % plot streak ratio for each coil element
end

%find the coil elements whose streak ratio greater than 1.3
%unstreaking is performed only for these coils, as described in the paper
Coil=find(StreakRatio>1.3)';

% Do unstreaking
StreakRatio=repmat(StreakRatio(Coil),[1,nx,nz,ntviews]);
StreakRatio=permute(StreakRatio,[2,4,3,1]);
kdata(:,:,:,Coil)=kdata(:,:,:,Coil)./StreakRatio;

% reconstruct images again for comparison
Ref=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n1,nx/2));               % artifact-free image
Img=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n2,nx/2)*n1/n2);         % image with streaks

data(:,:,:,3)=sos(Ref,4);
data(:,:,:,4)=sos(Img,4);
data=data./max(data(:));

% Note that the images displayed below are before iterative reconstruction
% Thus, they both have strearking artifacts.
% However, note that the lower image has significantly less streaks
if doFigs
    subplot(2,1,1)
    imagesc(abs(data(:,:,8,2))),axis off, axis square, colormap(gray),title('before unstreaking')
    subplot(2,1,2)
    imagesc(abs(data(:,:,8,4))),axis off, axis square, colormap(gray),title('after unstreaking')
end
