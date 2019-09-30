%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This demo demonstrates the detection of the bolus signal as described in
%the paper

% Li Feng, NYU, 12/18/2017
% Li.Feng@nyumc.org

clear
clc
addpath Data
addpath nufft_files
addpath utils
close all
load Data/kdata_Unstreaking_CoilCompression.mat
load Data/Trajectory.mat
load Data/b1.mat % coil sensitivite maps

%%%%%%%%%%%%%%%%%%%%%%
b1=double(b1);
nline=55; % 55 spokes per frame
step=9;   % step size 9 spokes
ndiff=nline-step;
close all

[nx,ntviews,nz,nc]=size(kdata);
nt=40;% reconstruct only 40 frames for saving time
Cut=5;% discard the first 5 frames
kdata_Under=zeros(nx,nline,nz,nc,nt);
Traj_Under=zeros(nx,nline,nt);
DensityComp_Under=zeros(nx,nline,nt);

filter=kaiser(nx,2);% do some filtering, for bolus tracking only.
kdata1=kdata.*repmat(filter,[1,ntviews,nz,nc]);

% Sort the k-space data
for ii=Cut+1:Cut+nt
    kdata_Under(:,:,:,:,ii-Cut)=kdata1(:,(ii-1)*step+1:ii*step+ndiff,:,:);
    Traj_Under(:,:,ii-Cut)=Traj(:,(ii-1)*step+1:ii*step+ndiff);
    DensityComp_Under(:,:,ii-Cut)=DensityComp(:,(ii-1)*step+1:ii*step+ndiff);
end

E=MCNUFFT3D(Traj_Under,DensityComp_Under,b1);
y=kdata_Under.*repmat(sqrt(permute(DensityComp_Under,[1 2,4,5,3])),[1,1,nz,nc,1]);
recon=abs(E'*y);
recon=recon(65:end-64,65:end-64,:,:);
recon=recon./max(recon(:));

% Do some spatial smooth filtering on the reconstructed images
% This is used for the bolus detection ONLY
for ii=1:size(recon,4)
    for jj=1:size(recon,3)
        recon(:,:,jj,ii)=medfilt2(recon(:,:,jj,ii),[3,3]);
    end
end

%%%%%%%%%%%%%%%
%find the best mask for each selected slice
clear maskAO
for tt=1:size(recon,3);
    clear circularity
    
    %Height filter
    data=squeeze(recon(:,:,tt,:));
    [nx,ny,nt]=size(data);
    [Y,I] = max(data,[],3);
    height=Y-mean(data(:,:,1:3),3);
    
    %Slope filter
    [Y,I] = max(data,[],3);
    Rt=Y./mean(data(:,:,1:3),3);
    slope=(Rt-1)./I;
    slope=slope./max(slope(:));
    
    mask=slope.*height;
    mask(find(isnan(mask)==1))=0;
    mask=mask>max(mask(:))*0.15;
    
    % labelling process as described in the paper
    % to find the mask that matches the aorta.
    Label=bwlabel(mask);
    for NumObj=1:max(Label(:))
        tmp=(Label==NumObj);
        tmp=imfill(tmp,'holes');
        if sum(tmp(:))<20 | sum(tmp(:))>300
            mask(Label==NumObj)=0;
        end
    end
    Label=bwlabel(mask);
    for NumObj=1:max(Label(:))
        tmp=(Label==NumObj);
        tmp=imfill(tmp,'holes');
        Label(find(tmp~=0))=NumObj;
    end
    if sum(Label(:))~=0
        [B,L] = bwboundaries(mask,'noholes');
        stats = regionprops(L,'Area','Centroid','MajorAxisLength','MinorAxisLength');
        for NumObj = 1:length(B)
            boundary = B{NumObj};
            delta_sq = diff(boundary).^2;
            perimeter = sum(sqrt(sum(delta_sq,2)));
            area = stats(NumObj).Area;
            circularity(NumObj) = 4*pi*area/perimeter^2*(stats(NumObj).MinorAxisLength/stats(NumObj).MajorAxisLength);
        end
        mask=(Label==find(circularity==max(circularity)));
    else
        mask=zeros(nx,ny);
    end
    maskAO(:,:,tt)=mask;
end
maskAO=repmat(maskAO,[1,1,1,nt]);

Signal_Cluster=squeeze(sum(sum(abs(recon.*maskAO),1),2));
Signal_Cluster=Signal_Cluster';

%Clusting to find the bolus signal (the dominate signal pattern)
thresh = 0.97;
[Signal, cluster] = CoilClustering(Signal_Cluster, thresh);
tmp=find(cluster~=0);
while length(tmp)<=2
    thresh=thresh-0.01;
    [Signal, cluster] = CoilClustering(Signal_Cluster, thresh);
    tmp=find(cluster~=0);
end

% Find the peak of the bolus signal, which will be the point of peak
% contrast enhancement
tmp=diff(Signal);tmp=tmp/max(tmp(:));
BolusPeak_Index=find(tmp==1);
while tmp(BolusPeak_Index+1)>0.2
    BolusPeak_Index=BolusPeak_Index+1;
end
BolusPeak_Index=BolusPeak_Index+1;

figure,plot(Signal);
hold on; plot(BolusPeak_Index,Signal(BolusPeak_Index),'ro');

%calculate the number of spoke that reaches the peak enhancement
AOSpoke=(BolusPeak_Index+Cut-1)*step+ceil(step/2);
save Data/Bolus.mat AOSpoke Signal BolusPeak_Index recon Signal_Cluster
close all