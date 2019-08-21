%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Matlab source code for GRASP DCE-liver MRI as described in:
%  
%  Feng L, Grimm R Block KT, Chandarana H, Kim S, Xu J, Axel L, Sodickson DK, Otazo R. 
%  Golden-angle radial sparse parallel MRI: Combination of compressed sensing, 
%  parallel imaging, and golden-angle radial sampling for fast and flexible dynamic 
%  volumetric MRI
%  Magn Reson Med. 2014 Sep;72(3):707-17
% 
%
%  The source code uses the following external packages:
%    - NUFFT toolkit by Jeffrey Fessler 
%      (http://www.eecs.umich.edu/~fessler/)
%    - Non-linear conjugate gradient algorithm by Miki Lustig
%      (http://www.eecs.berkeley.edu/~mlustig/Software.html)

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
load data_DCE.mat
% kc is not used in this demo as respiratory motion signal is not needed.

% kdata is the k-space data. It is only one slice selected from the 3D stack-of-stars
% datasets after a FFT along the kz. In order to save recon time and
% memory, the full stack-of-stars data is not included. Email me if you
% want a fully 3D dataset.

% b1 is the coil sensitivity maps of the selected slice

% k is the radial k-space trajectory and w is the corresponding density compensation

[nx,ntviews,nc]=size(kdata);
%ntviews: number of acquired spokes
%nc: number of coil elements
%nx: readout point of each spoke (2x oversampling included)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data sorting

nline=34; % number of spokes in each ontrast-enhanced phase
nt=floor(ntviews/nline);% number of contrast-enhanced phases

kdata=kdata.*repmat(sqrt(w),[1,1,nc]);

clear kdata_u k_u w_u Res_Signal_u 
for ii=1:nt
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
param.TVWeight_dim1=max(abs(recon_cs(:)))*0.03;
param.TVWeight_dim2=0;
param.nite = 5;param.display=1;

clc
tic
for n=1:3
    recon_cs = CSL1NlCg_XDGRASP(recon_cs,param);
end
time=toc;
time=time/60
data_grasp=recon_cs/max(abs(recon_cs(:)));

figure, imshow3(abs(data_gridding(:,:,1:4:end)),[0 .5],[1,8]);title('Grdding')
figure, imshow3(abs(data_grasp(:,:,1:4:end)),[0 .5],[1,8]);title('GRASP')


