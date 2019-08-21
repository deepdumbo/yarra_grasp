function  res = TV_Temp_Dim2()

% inite-differencing operator along respiratory dimension.
%
% Li Feng 2016

res.adjoint = 0;
res = class(res,'TV_Temp_Dim2');

