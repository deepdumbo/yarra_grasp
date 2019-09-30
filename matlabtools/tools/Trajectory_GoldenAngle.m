%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to generate 2D golden-angle radial trajectory
%
% Marnix Maas, RadboudUMC, 5-4-2018
% Adapted from Li Feng's Trajectory_GoldenAngle_GROG.m
% 
% Differences with GROG trajectory:
% - Optional normalization of k-space coordinates by flag normalize.
% Default is normalize=1.
% - Optional flipping of rotation angle and direction of k-space traversal
% along each spoke, by 2-element flag vector flip = [flipRotationAngle
% flipSpokeTraversal]. Default is no flips (ie flip = [0 0]).
% - Optional gradient correction factor gradcorr, according to Tobias
% Block's iterativeradial_multicoil.m. Default is no gradient correction
% (ie gradcorr = 1)
%
% USE: Traj=Trajectory_GoldenAngle(ntviews, nx, normalize=1, flip=[0 0], gradcorr=1)
% ivom: I think normalization is incorrect for GROG! Not done in demo, and
% wrong results in my cases. 

function Traj=Trajectory_GoldenAngle(ntviews, nx, normalize, flip, gradcorr)
if nargin<5
    gradcorr = 1;
end
if nargin<4
    flip = [0 0];
end
if nargin<3
    normalize = 1;
end

Gn = (1 + sqrt(5))/2;

% Calculate angles
% radian = mod((0:(ntviews-1))*a*pi/180 + pi/2,                                                                                                                                                                                                                                                                                                                                      2*pi);     % Added pi/2
phi = mod( pi/2 + (0:(ntviews-1))*pi/Gn, 2*pi);         % Angles in radians

% Calculate coordinates along each spoke
rho = linspace(0,nx-1,nx) - (nx-1)/2;
% rho = [-floor(nx/2):floor(nx/2)];
% rho=rho(1:nx)+0.5;
if gradcorr
    rho = rho - (1-gradcorr)/2;                         % This reproduces behavior of Tobias Block's iterativeradial_multicoil.m; 
                                                        % not necessary to make it this complex, but keeps things numerically the same
end

if normalize
    rho = rho/nx;                                       % Normalization
end

if flip(2)
    rho = -rho;
end

if flip(1)
    Traj = double(rho'*exp(-1j*phi));
else
    Traj = double(rho'*exp(1j*phi));
end

return