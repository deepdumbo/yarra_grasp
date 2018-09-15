%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to generate 2D golden-angle radial trajectory, optionally with
% multiple echoes
%
% Marnix Maas, RadboudUMC, 30-8-2018
% Adapted from Li Feng's Trajectory_GoldenAngle_GROG.m
% 
% Parameters:
% - ntviews:      number of radial vies
% - nx:           number of points along each spoke
% - ne:           number of echoes (default=1)
% - blipAngleDeg: blip angle between spokes in degrees
% 
% Optional parameters can be set using varargin:
% - 'normalize': flag for normalization of k-space coordinates (default 1)
% - 'flip':      flipping of rotation angle and direction of k-space traversal
% along each spoke, by 2-element flag vector flip = [flipRotationAngle
% flipSpokeTraversal]. Default is no flips (ie flip = [0 0]).
% - 'gradcorr':  gradient correction factor according to Tobias Block's 
% iterativeradial_multicoil.m. Default is no gradient correction
% (ie gradcorr = 1)
% 
%
% USE: Traj=Trajectory_GoldenAngle_ME(ntviews, nx, ne, blipAngleDeg,
% varargin)

function Traj=Trajectory_GoldenAngle_ME(ntviews, nx, ne, blipAngleDeg, varargin)

% Parse input arguments
p = parseInputs(ntviews, nx, ne, blipAngleDeg, varargin);
parse(p, ntviews, nx, ne, blipAngleDeg, varargin{:});

gradcorr    = p.Results.gradcorr;
flip        = p.Results.flip;
normalize   = p.Results.normalize;

% if nargin<3
%     ne = 1;
% end

Gn = (1 + sqrt(5))/2;

% Calculate angles
phi = mod( pi/2 + (0:(ntviews-1))*pi/Gn, 2*pi);         % Angles in radians

% Calculate coordinates along each spoke 
rho = linspace(0,nx-1,nx) - (nx-1)/2;

% Optional: apply gradient correction
if gradcorr
    rho = rho - (1-gradcorr)/2;                         % This reproduces behavior of Tobias Block's iterativeradial_multicoil.m; 
                                                        % not necessary to make it this complex, but keeps things numerically the same
end
% Optional: Normalize spoke length
if normalize
    rho = rho/nx;
end
% Optional: Swap traversal direction of spokes
if flip(2)
    rho = -rho;
end

% Create 2-dimensional trajectory (1st echo)
if flip(1)
    Traj = double(rho'*exp(-1j*phi));
else
    Traj = double(rho'*exp(1j*phi));
end

% Make it 3D, with echoes in the 3rd dimension
Traj = repmat(Traj, [1 1 ne]);

% Create vector of blip angles to apply to different echoes
dPhi = pi*blipAngleDeg/180 * ((1:ne)-1);       % deg to rad

for i=1:ne
    Traj(:,:,i) = Traj(:,:,i) * exp(1i*dPhi(i)) * exp(1i*pi*(i-1));
    % Last factor is necessary because even echoes are traversed in
    % opposite direction compared to odd ones.
end

% Multiply trajectory with blip angles
% Traj = bsxfun(@times,Traj,dPhi);

% Traj = zeros(nx, ntviews, ne);
% for i = 1:ne
%     if pars.flip(1)
%         Traj(:,:,i) = double(rho(i,:)'*exp(-1j*phi(i,:)));
%     else
%         Traj(:,:,i) = double(rho(i,:)'*exp(1j*phi(i,:)));
%     end
% end

end % Trajectory_GoldenAngle_ME()

function p = parseInputs(ntviews, nx, ne, blipAngleDeg, varargin)
p = inputParser;
% default values
defNe               = 1;
defBlipAngleDeg     = 0;
defNormalize        = 1;
defFlip             = [0 0];
defGradcorr         = 1;

% fill parser object
addRequired(p,'ntviews',@isnumeric);
addRequired(p,'nx',@isnumeric);
addOptional(p,'ne',defNe,@isnumeric);
addOptional(p,'blipAngleDeg',defBlipAngleDeg,@isnumeric);
addParameter(p,'normalize',defNormalize,@isnumeric);
addParameter(p,'flip',defFlip);
addParameter(p,'gradcorr',defGradcorr,@isnumeric);
end