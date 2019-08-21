%Motion detection: step 1

%Li Feng, NYU, 12/18/2017

function [Coil,Res_Signal_Post]=MC_Step1(ZIP,n1);

ntviews=size(ZIP,2);
close all
k=1;N=n1;
clear Nav

% Do PCA along each coil element and to get the first two PCs from each
% coil, as described in the paper
%ivom: Loop 2 data sequences projected onto principal components (k)
%ivom: for each of 20 (jj) coils.  
for ii=1:size(ZIP,3)
    SI=ZIP(:,end-N+1:end,ii)';
    covariance=cov(SI);
    [PC, V] = eig(covariance);
    V = diag(V);
    [junk, rindices] = sort(-1*V); %ivom: or sort(V, 'descend')
    V = V(rindices);
    %ivom: PC are the new basis vectors with most dynamics.
    %ivom: PCs represent linear combinations of the time-series at z-points
    %ivom: which are 'relevant'. 
    PC = PC(:,rindices);
    %ivom: size(SI) = [800 spokes, 400 samples/spoke]
    SI = (PC' * SI')';
    for jj=1:2
        %ivom: smooth data and save to array Nav
        tmp=smooth(SI(:,jj),6,'lowess');
        tmp=tmp-min(tmp(:));
        tmp=tmp./max(tmp(:));
        %ivom: this works! creates new Nav array. Nav contains
        %ivom: data along first 2 PC's for 20 coils. 
        %ivom: why do some coils show such different behaviour?
        %ivom: (Even anti-correlation??) Are those coils in opposite
        %ivom: positions?
        Nav(:,k)=tmp;k=k+1;
    end
end
figure,imagesc(abs(ZIP(:,end-n1+1:end,15))),axis image, axis off, colormap(gray),title('Respiratory Motion')
hold on
plot(Nav(:,1)*100+220, 'b')
close all

% coil clustering
%ivom: coil clustering is an SVD from a cross-correlation at zero
%ivom: time lag of all 40 time series (20 coils, 2 principal components).
%ivom: The covariance matrix is almost fully diagonal, with ~5
%ivom: off-diagonal elements. 
%ivom: Based of SOM of Feng paper, I think coil clustering selects the time curves
%ivom: That are very highly correlated and averages these as the
%ivom: respiratory signal. 
% code obtained from Tao Zhang, stanford
thresh = 0.97;
[Res_Signal, cluster] = CoilClustering(Nav, thresh);
tmp=find(cluster~=0);
%ivom: if the threshold is too high and no cluster is found, decrease it. 
while length(tmp)<=2
    thresh=thresh-0.01;
    [Res_Signal, cluster] = CoilClustering(Nav, thresh);
    tmp=find(cluster~=0);
end
thresh
%ivom: set minimum to zero for Res_Signal array
Res_Signal=Res_Signal-min(Res_Signal(:));
%ivom: set maximum to one through normalization
Res_Signal_Post=Res_Signal./max(Res_Signal(:));

% find the "good" coil elements used for estimation of all motion later
cluster=abs(reshape(cluster,[2,size(Nav,2)/2]));
cluster=sum(cluster,1);
Coil=find(cluster>0);

% if mean(Res_Signal_Post)<0.5
%     Res_Signal_Post=-Res_Signal_Post;
%     Res_Signal_Post=Res_Signal_Post-min(Res_Signal_Post(:));
%     Res_Signal_Post=Res_Signal_Post./max(Res_Signal_Post(:));
% end

figure,imagesc(abs(ZIP(:,end-n1+1:end,15))),axis image, axis off, colormap(gray),title('Respiratory Motion')
hold on
plot(-Res_Signal_Post(:)*100+220,'r')
