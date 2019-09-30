%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Matlab source code for XD-GRASP liver MRI 
%  as described in:
%  
%  Feng L, Axel L, Chandarana H, Block KT, Sodickson DK, Otazo R. 
%  XD-GRASP: Golden-angle radial MRI with reconstruction of 
%  extra motion-state dimensions using compressed sensing
%  Magn Reson Med. 2016 Feb;75(2):775-88
% 
%
%  The source code uses the following external packages:
%    - NUFFT toolkit by Jeffrey Fessler 
%      (http://www.eecs.umich.edu/~fessler/)
%    - Non-linear conjugate gradient algorithm by Miki Lustig
%      (http://www.eecs.berkeley.edu/~mlustig/Software.html)
%    - Coil clustering code by Tao Zhang 
%      (http://web.stanford.edu/~tzhang08/software.html)

%  If you use this code, please cite the above publication.
%
%  (c) Li Feng, 2016, New York University
%  Li.Feng@nyumc.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('Dataset');  
addpath('Code');   
addpath('Code/nufft_files');

clear
clc

%load data
load data_NonContrast.mat
% kc is the central profiles of the stack-of-stars k-space (kx=ky=0).

% kdata is the k-space data. It is only one slice selected from the 3D stack-of-stars
% datasets after a FFT along the kz. In order to save recon time and
% memory, the full stack-of-stars data is not included. Email me if you
% want a fully 3D dataset.

% b1 is the coil sensitivity maps of the selected slice

% k is the radial k-space trajectory and w is the corresponding density compensation

[nz,ntviews,nc]=size(kc);
nx=size(kdata,1);
%nz: number of slice
%ntviews: number of acquired spokes
%nc: number of coil elements
%nx: readout point of each spoke (2x oversampling included)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Respiratory motion detection

% Generate the z-projection profiles
%ZIP: Projection profiles along the Z dimension with interpolation.
ZIP=abs(fftshift(ifft(kc,400,1),1));

%Normalization of each projection in each coil element
for ii=1:nc
    for jj=1:ntviews
        maxprof=max(ZIP(:,jj,ii));
        minprof=min(ZIP(:,jj,ii));
        ZIP(:,jj,ii)=(ZIP(:,jj,ii)-minprof)./(maxprof-minprof);
    end
end
figure,imagesc(abs(ZIP(:,:,5))),axis image,colormap(gray), axis off

% Perform PCA on each coil element
close all
kk=1;clear PCs
for ii=1:nc
    tmp=permute(ZIP(:,:,ii),[1,3,2]);
    tmp=abs(reshape(tmp,[size(tmp,1)*size(tmp,2),ntviews])');
    
    covariance=cov(tmp);
    [tmp2, V] = eig(covariance);
    V = diag(V);
    [junk, rindices] = sort(-1*V);
    V = V(rindices);
    tmp2 = tmp2(:,rindices);
    PC = (tmp2' * tmp')';
    
    % Take the first two principal components from each coil element.
    for jj=1:2
        tmp3=smooth(PC(:,jj),6,'lowess'); % do some moving average smoothing
        
        %Normalize the signal for display
        tmp3=tmp3-min(tmp3(:));
        tmp3=tmp3./max(tmp3(:));
        PCs(:,kk)=tmp3;kk=kk+1;
%         %plot the estimated signal
%         imagesc(abs(ZIP(:,:,ii))),axis image,colormap(gray), axis off
%         hold on
%         plot(-tmp3(:)*100+220,'r')
%         hold off
%         pause
    end
end    
               
close all
% Do coil clusting to find the respiratory motion signal
% Function obtained from Tao Zhang (http://web.stanford.edu/~tzhang08/software.html)
thresh = 0.95;
[Res_Signal, cluster] = CoilClustering(PCs, thresh);

%Normalize the signal for display
Res_Signal=Res_Signal-min(Res_Signal(:));
Res_Signal=Res_Signal./max(Res_Signal(:));

% Plot the respiratory motion signal on the projection profiles
imagesc(abs(ZIP(:,:,5))),axis image,colormap(gray), axis off,title('Respiratory Motion')
hold on
plot(-Res_Signal(:)*100+220,'r')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data sorting

ntres=4;% number of respiratory phases
nline=floor(ntviews/ntres);

% Sort the k-space data and trajectory according to respiratory motion 
[~,index]=sort(Res_Signal,'descend');
kdata=kdata(:,index,:);
k=k(:,index,:);
w=w(:,index);
kdata=kdata.*repmat(sqrt(w),[1,1,nc]);

clear kdata_u k_u w_u Res_Signal_u 
for ii=1:ntres
    kdata_u(:,:,:,ii)=kdata(:,(ii-1)*nline+1:ii*nline,:);
    k_u(:,:,ii)=k(:,(ii-1)*nline+1:ii*nline);
    w_u(:,:,ii)=w(:,(ii-1)*nline+1:ii*nline);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recon

param.E=MCNUFFT(double(k_u),double(w_u),double(b1));
param.y=double(kdata_u);
recon_cs=param.E'*param.y;
data_gridding=recon_cs/max(abs(recon_cs(:)));

param.TV_dim1=TV_Temp;
param.TVWeight_dim1=max(abs(recon_cs(:)))*0.02;
param.TVWeight_dim2=0;
param.nite = 4;param.display=1;

clc
tic
for n=1:2
    recon_cs = CSL1NlCg_XDGRASP(recon_cs,param);
end
time=toc;
time=time/60
data_xdgrasp=recon_cs/max(abs(recon_cs(:)));

figure, imshow3(abs(data_gridding),[0 .8],[1,4]);title('Grdding Motion State 1-4')
figure, imshow3(abs(data_xdgrasp),[0 .8],[1,4]);title('XD-GRASP Motion State 1-4')


