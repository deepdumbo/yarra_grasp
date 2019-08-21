function  res = Emat_GROG2DSG(mask,b1,W1)

% 2D GROG operator
% Input
% mask: k-space mask
% b1: coil sensitivity maps
% W1: GROG weighting matrix

%Li Feng, NYU, 12/18/2017

res.adjoint = 0;
res.mask = mask;
res.b1=repmat(permute(b1,[1,2,4,3]),1,1,size(mask,3),1);
res.W1=sqrt(W1);
res = class(res,'Emat_GROG2DSG');

