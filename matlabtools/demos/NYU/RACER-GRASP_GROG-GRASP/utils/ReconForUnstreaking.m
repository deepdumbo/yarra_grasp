%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to perform simple gridding reconstruction
% Used for streak ratio calculation

%Li Feng, NYU, 12/18/2017

function ref=ReconForUnstreaking(kdata,Traj,DensityComp,N,Bas);

Traj=Traj(:,end-N+1:end);
DensityComp=DensityComp(:,end-N+1:end);;
kdata=kdata(:,end-N+1:end,:,:);
[nx,ntviews,nz,nc]=size(kdata);
kdata=double(kdata.*repmat(sqrt(DensityComp),[1,1,nz,nc]));

Img_Dim=[Bas,Bas,nz];
param.E = MCNUFFT3D(double(Traj),double(DensityComp),ones(Img_Dim));clear ref
% Pre-allocate memory for ref here. How much is needed?
for ch=1:nc
    ch
    ref(:,:,:,ch)=param.E'*kdata(:,:,:,ch);
end
