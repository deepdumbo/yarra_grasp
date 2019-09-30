%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Matlab source code for XD-GRASP DCE-MRI of the liver
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

clear
clc

%load data
load data_DCE.mat
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: find the coil elements with good respiratory motion display
%         from the late enhancement spokes
ZIP1=ZIP(:,501:end,:);
[nz,ntviews,nc]=size(ZIP1);

% Perform PCA on each coil element
close all
kk=1;clear PCs
for ii=1:nc
    ii
    tmp=permute(ZIP1(:,:,ii),[1,3,2]);
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
        tmp3=smooth(PC(:,jj),6,'lowess');
        
        %Normalize the signal for display
        tmp3=tmp3-min(tmp3(:));
        tmp3=tmp3./max(tmp3(:));
        PCs(:,kk)=tmp3;kk=kk+1;
%         %plot the estimated signal
%         imagesc(abs(ZIP1(:,:,ii))),axis image,colormap(gray), axis off
%         hold on
%         plot(-tmp3(:)*100+220,'r')
%         hold off
%         pause
    end
end    
               
close all
% Coil clusting to find the respiratory motion signal
thresh = 0.97;
[Res_Signal, cluster] = CoilClustering(PCs, thresh);

%Normalize the signal for display
Res_Signal=Res_Signal-min(Res_Signal(:));
Res_Signal=Res_Signal./max(Res_Signal(:));
imagesc(abs(ZIP1(:,:,5))),axis image,colormap(gray), axis off
hold on
plot(-Res_Signal(:)*100+220,'r')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
hold off

cluster=abs(cluster(1:2:end))+abs(cluster(2:2:end));
% cluster gives a coil index with ~=0 indicating "good" coils with good
% respiratory motion display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Estimating respiratory motion from the "good" coil elements

%Perform PCA on the stack of "good" coil elements
[nz,ntviews,nc]=size(ZIP);
SI=permute(ZIP(:,:,[find(cluster~=0)]),[1,3,2]);
SI=abs(reshape(SI,[size(SI,1)*size(SI,2),ntviews])');

covariance=cov(SI);
[PC, V] = eig(covariance);
V = diag(V);
[junk, rindices] = sort(-1*V);
V = V(rindices);
PC = PC(:,rindices);
SI = (PC' * SI')';

% %Check the first 3 components
% for ii=1:3
%     figure,plot(SI(:,ii))
%     pause
% end
% close all

%Do some smoothing
for ii=1:3
    SI(:,ii)=smooth(SI(:,ii),6,'lowess');
end

% Pick the best principal component
%ideally, this could be picked by evaluating the frequency peak of each PC
Res_Signal=SI(:,3);

%Normalization for display
Res_Signal=-Res_Signal; % For this dataset, the respiratory motion is up-side down
                        % so put it back. This can be automatically
                        % detected since the end-expiratory phases has
                        % stronger signal intensity in the entire FOV due
                        % to the decreased size of lung volume

Res_Signal=Res_Signal-min(Res_Signal(:));
Res_Signal=Res_Signal./max(Res_Signal(:));
imagesc(abs(ZIP(:,:,6))),axis image,axis image,colormap(gray), axis off
hold on
plot(-Res_Signal(:)*100+220,'r')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

% Estimate the envelope of the signal (contrast enhancement + respiration)
% This part needs to be further improved
ft = fittype( 'smoothingspline' );
opts = fitoptions( ft );
opts.SmoothingParam = 0.015;

t=16;
[Peak,Peak_Index]=findpeaks(double(Res_Signal),'MINPEAKDISTANCE',t);
figure,plot(Res_Signal); 
hold on; plot(Peak_Index,Res_Signal(Peak_Index),'ro');

% substract the estimated envelope
[xData, yData] = prepareCurveData(Peak_Index,Res_Signal(Peak_Index));
[fitresult, gof] = fit( xData, yData, ft, opts );
cfval = coeffvalues(fitresult);
ftmax = feval(ft,cfval(1),(1:ntviews)');
Res_Signal=Res_Signal-ftmax;

[Peak,Peak_Index]=findpeaks(double(Res_Signal),'MINPEAKDISTANCE',t);
figure,plot(Res_Signal); 
hold on; plot(Peak_Index,Res_Signal(Peak_Index),'ro');

% display the final signal on the projections
imagesc(abs(ZIP(:,:,5))),axis image,axis image,colormap(gray), axis off
hold on
plot(-Res_Signal*150+100,'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data sorting

nline=100 % number of spoeks for each contrast-enhanced phase
nt=floor(ntviews/nline); % number of contrast-enhanced phases
ntres=4; % number of respiratory phases
nline2=floor(nline/ntres); % number of spokes in each phases after respiratory sorting

kdata=kdata.*repmat(sqrt(w),[1,1,nc]);

clear kdata_u k_u w_u Res_Signal_u 
for ii=1:nt
    kdata_u(:,:,:,ii)=kdata(:,(ii-1)*nline+1:ii*nline,:);
    k_u(:,:,ii)=k(:,(ii-1)*nline+1:ii*nline);
    w_u(:,:,ii)=w(:,(ii-1)*nline+1:ii*nline);
    Res_Signal_u(:,ii)=Res_Signal((ii-1)*nline+1:ii*nline);
end

clear kdata_u1 k_u1 w_u1  
for ii=1:nt
    tmp1=kdata_u(:,:,:,ii);
    tmp2=k_u(:,:,ii);
    tmp3=w_u(:,:,ii);
    [~,index]=sort(Res_Signal_u(:,ii),'descend');
    tmp1=tmp1(:,index,:);
    tmp2=tmp2(:,index);
    tmp3=tmp3(:,index);
    for jj=1:ntres
        kdata_u1(:,:,:,jj,ii)=tmp1(:,(jj-1)*nline2+1:jj*nline2,:);
        k_u1(:,:,jj,ii)=tmp2(:,(jj-1)*nline2+1:jj*nline2,:);
        w_u1(:,:,jj,ii)=tmp3(:,(jj-1)*nline2+1:jj*nline2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recon
param.E=MCNUFFT_DCE(double(k_u1),double(w_u1),double(b1));
param.y=double(kdata_u1);
recon_cs=param.E'*param.y;
data_gridding=recon_cs/max(abs(recon_cs(:)));

param.TV_dim1=TV_Temp_Dim1;
param.TVWeight_dim1=max(abs(recon_cs(:)))*0.03;
param.TV_dim2=TV_Temp_Dim2;
param.TVWeight_dim2=max(abs(recon_cs(:)))*0.015;
param.nite = 5;param.display=1;

clc
tic
for n=1:3
    recon_cs = CSL1NlCg_XDGRASP(recon_cs,param);
end
time=toc;
time=time/60
data_xdgrasp=recon_cs/max(abs(recon_cs(:)));

% Display the end-expiratory phase
figure, imshow3(abs(squeeze(data_gridding(:,:,1,1:6))),[0 .8],[1,6]);title('Grdding')
figure, imshow3(abs(squeeze(data_xdgrasp(:,:,1,1:6))),[0 .8],[1,6]);title('XD-GRASP')

% Display different motion states
figure, imshow3(abs(squeeze(data_xdgrasp(:,:,:,5))),[0 .8],[1,4]);title('XD-GRASP')
