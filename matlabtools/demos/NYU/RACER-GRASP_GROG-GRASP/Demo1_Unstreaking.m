%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This demo demonstrates the calculation of streak ratio of each coil
%element and the unstreaking process for this sample data

%Note that only 12 slices of the data were uploaded
%for demonstration purpose due the space limit.

% Li Feng, NYU, 12/18/2017
% Li.Feng@nyumc.org

clear
clc
addpath Data
addpath nufft_files
addpath utils
load Data/kdata.mat
load Data/Trajectory.mat

n1=1000;% using 1000 spokes to generate artifact free images
n2=40; % using 40 spokes to generate images with streaking artifacts

[nx,ntviews,nz,nc]=size(kdata);
% Note that the third dimension is z, NOT kz. (A FFT was performed already)

Ref=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n1,nx)); % artifact-free image
Img=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n2,nx)*n1/n2); % image with streaks

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
figure,plot(StreakRatio) % plot streak ratio for each coil element

%find the coil elements whose streak ratio greater than 1.3
%unstreaking is performed only for these coils, as described in the paper
Coil=find(StreakRatio>1.3)';

% Do unstreaking
StreakRatio=repmat(StreakRatio(Coil),[1,nx,nz,ntviews]);
StreakRatio=permute(StreakRatio,[2,4,3,1]);
kdata(:,:,:,Coil)=kdata(:,:,:,Coil)./StreakRatio;

% reconstruct images again for comparison
Ref=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n1,nx/2)); % artifact-free image
Img=abs(ReconForUnstreaking(kdata,Traj,DensityComp,n2,nx/2)*n1/n2); % image with streaks

data(:,:,:,3)=sos(Ref,4);
data(:,:,:,4)=sos(Img,4);
data=data./max(data(:));

% Note that the images displayed below are before iterative reconstruction
% Thus, they both have strearking artifacts.
% However, note that the lower image has significantly less streaks
subplot(2,1,1)
imagesc(abs(data(:,:,8,2))),axis off, axis square, colormap(gray),title('before unstreaking')
subplot(2,1,2)
imagesc(abs(data(:,:,8,4))),axis off, axis square, colormap(gray),title('after unstreaking')


%Now compress the unstreaked kpsace kdata
ncc=8;% compress to 8 elements
D=reshape(kdata,nx*nz*ntviews,nc);
[U,S,V]=svd(D,'econ');
kdata=single(reshape(D*V(:,1:ncc),nx,ntviews,nz,ncc));

save Data/kdata_Unstreaking_CoilCompression.mat kdata
