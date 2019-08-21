function out = modelbased_fw(par,kdata)

%% Extract parameters from struct
Nread           = par.Nread;
Nproj           = par.Nproj;
Ncoil           = par.Ncoil;
Neco            = par.Neco;

%% Scale kdata
fprintf('\t Scaling k-space data \n')
par.factor      = 100 / sqrt(kdata(:)'*kdata(:));
kdata           = kdata * par.factor;

%% Calculate trajectory
fprintf('\t Calculate trajectory \n')
phi     = zeros(Nread,Nproj,Neco);

GA = 111.246117975/180*pi;
for ip = 1:Nproj
    for ie = 1:Neco
        phi(:,ip,ie) = pi/2 + (ip-1)*GA + (ie-1)*par.blipangle;
    end
end

rho     = linspace(0,Nread-1,Nread)' - (Nread-1)/2;
rho     = repmat(rho,[1,Nproj,Neco]);
rho     = rho/Nread;

% Generate vector with k-space coordinates (as complex values kx + i ky)
traj    = double(rho.*exp(1j*phi));

% Calculate DCF
dcf     = ones(Nread,1);

for i=1:Nread
    dcf(i) = abs(Nread/2 - (i - 0.5));
end

dcf     = pi/Nproj*dcf;
dcf     = repmat(dcf,1,Nproj);

%% Calculate NUFFT operator
fprintf('\t Calculating NUFFT operators \n')
for ie = 1:Neco
    par.FT{ie} = NUFFT(traj(:,:,ie), 1, 1, [0,0], par.imsize, 2);
end

%% Gridding reconstruction
img_uncombined = zeros([par.imsize, Ncoil, Neco],par.prec);
fprintf('\t Gridding of k-space data \n')

for ie = 1:Neco
    for ic = 1:Ncoil
        if (par.verbose == 1)
            fprintf('\t \t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)
        end
        kdata_tmp = dcf .* cast(squeeze(kdata(:,:,ic,ie)),par.prec);
        img_uncombined(:,:,ic,ie) = par.FT{ie}' * kdata_tmp;
    end
end

img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));

%% Scale lambda according to maximum of gridded images
par.FD1Weight  = par.FD1Weight * max(abs(img_gridding_sos(:)));

%% Coil combination
img_coilcombined = squeeze(sum(bsxfun(@times,img_uncombined,conj(par.coilmaps)),3));

%% Compute frequencies in Hz
gyro            = 42.58;
par.f_wf        = [];
par.rel_amp     = [];
for i = 1:length(par.species)
    par.f_wf    = [par.f_wf; gyro * par.fieldStrength * par.species(i).frequency(:)];
    par.rel_amp = [par.rel_amp; par.species(i).relAmps(:)];
end

%% Calculate fat modulation frequencies
fprintf('Calculating fat modulation frequencies \n');
par.fatModulation = calculateFatModulation(par);

if par.verbose
    figure(574), imagesc(abs(squeeze(par.fatModulation))), colorbar, colormap gray
end

%% Estimate initial field map
if (par.estfieldmap == 1)
    fprintf('\t Calculating fieldmap estimate \n')
    % Uses method described in Berglund J et al. Magn Reson Med. 2010;63(6):1659-68
    % http://www.onlinelibrary.wiley.com/doi/10.1002/mrm.22385/abstract
    
    % Code taken from ISMRM FatWater Toolbox
    % http://www.ismrm.org/workshops/FatWater12/data.htm
    
    f_wf_fmap       = par.f_wf(2:end) * 2*pi;
    rel_amp_fmap    = par.rel_amp(2:end) / sum(par.rel_amp(2:end));
    a               = exp(complex(0,par.TE(1:3)' * f_wf_fmap')) * rel_amp_fmap;
    A               = [ones(3,1) a];
    S               = permute(img_coilcombined,[1 2 4 3]);
    
    if ~isa(S,'single')
        S = single(S); % Necessary for execution of C++ code
    end
    
    [bA, bB] = getPhasorCandidates(S,A);
    mw              = getMagnitudeWeight(S);
    phasorexp       = regionGrow(A,bA,bB,par.c1,par.c2,mw,par.voxelSize);
    fieldmap_init   = phasorexp/(2*pi*par.dTE);
    
    if (par.smoothfm == 1)
        fprintf('\t Smoothing of fieldmap estimate \n')
        mask = [0 1 0; 1 4 1; 0 1 0]/8;
        for j = 1:par.smoothfmiter
            fieldmap_init = conv2(fieldmap_init,mask,'same');
        end
    end
else
    fieldmap_init   = zeros(par.imsize);
end

if (par.verbose == 1)
    figure(575), imagesc(fieldmap_init); axis image; axis off; colorbar
end

%% Define further optimization parameters
par.alpha         = 0.01;
par.beta          = 0.6;
par.maxIter_ls    = 100;
par.maxIter_cg    = 10;
par.threshold     = 0;

par.pnorm_im      = 1;
par.pnorm_off     = 2;

par.FD1           = FD1OP();
par.dataWeight    = 1;
par.mu            = 1e-15; % Smoothing parameter of L1 norm

par.alpha_n       = 1; % Regularization parameter for update dx (see Uecker MRM60:674-682 (2008))
par.q             = 2/3; % Reduction factor

%% Initialization
% x contains the water image, fat image and field map
x0(:,:,1:2) = complex(zeros([par.imsize,2]));
x0(:,:,3)   = fieldmap_init;

%% Perform minimization
clearvars -except kdata x0 par img_coilcombined
fprintf('\n \t Perform CS-WF iterations \n')
[water,fat,fieldmap,optinfo] = wfcs(kdata, x0, par);

out.water           = water;
out.fat             = fat;
out.fieldmap        = real(fieldmap);
out.images          = img_coilcombined;
out.optinfo         = optinfo;