function fieldmap_init = calc3Dfieldmap(kdata,par)

%% Extract parameters from struct
Nread           = par.Nread;
Nproj           = par.Nproj;
Ncoil           = par.Ncoil;
Neco            = par.Neco;
Npar            = par.Npar;

%% PCA for coil compression
% Should be performed partition-wise to avoid inconsistencies
kdata_new = zeros(Nread,Nproj,par.compressedcoils,Neco,Npar);

if (par.coilcompression == 1)
    fprintf('Coil compression, using %d virtual coil(s) \n',par.compressedcoils)
    for ip = 1:Npar
        fprintf('\t Npar: %d/%d \n',ip,Npar)
        
        kdata_part = kdata(:,:,:,:,ip);
        kdata_sum = sum(kdata_part,4);
        kdata_sum = reshape(kdata_sum,[Nread * Nproj, Ncoil]);
        
        % Calculate coefficients
        coeffs = pca(kdata_sum);
        kdata_tmp = permute(kdata_part,[1 2 4 3]);
        kdata_tmp = reshape(kdata_tmp,[Nread * Nproj * Neco, Ncoil]);
        
        % Perform change of basis
        kdata_tmp = kdata_tmp * coeffs;
        
        % Crop data
        Ncoil_new = par.compressedcoils;
        kdata_tmp = kdata_tmp(:,1:Ncoil_new);
        
        kdata_tmp = reshape(kdata_tmp,[Nread, Nproj, Neco, Ncoil_new]);
        kdata_tmp = permute(kdata_tmp,[1 2 4 3]);
        
        kdata_new(:,:,:,:,ip) = kdata_tmp;
    end
    
    kdata = kdata_new;
    clear kdata_tmp kdata_sum kdata_part;
    Ncoil = Ncoil_new;
    par.Ncoil = Ncoil_new;
    
end

%% Calculate trajectory data
fprintf('Calculate trajectory \n')
[traj, dcf] = calcTrajectory(par);

%% Calculate NUFFT operator
shift   = [0,0];
fprintf('Calculating NUFFT operators \n')
for ie = 1:Neco
    par.FT{ie} = NUFFT(traj(:,:,ie), 1, 1, shift, par.imsize, 2);
end

%% Perform gridding
img_uncombined = zeros([par.imsize, Ncoil, Npar, Neco]);

for ip = 1:Npar
    for ie = 1:Neco
        for ic = 1:Ncoil
            fprintf('Npar: %d/%d Neco: %d/%d Ncoil: %d/%d \n',ip,Npar,ie,Neco,ic,Ncoil)
            kdata_tmp = dcf.*squeeze(double(kdata(:,:,ic,ie,ip)));
            img_uncombined(:,:,ic,ip,ie) = par.FT{ie}' * kdata_tmp;
        end
    end
end

img_gridding_sos = makesos(img_uncombined,3);

%% Adaptive combine
img_sum = squeeze(sum(squeeze(img_uncombined),5));
img_sum = permute(img_sum,[3 1 2 4]);
[~,~,wmap] = openadapt(img_sum,1);
image = squeeze(sum(bsxfun(@times,squeeze(img_uncombined),permute(conj(wmap),[2 3 1 4])),3));

%% Compute frequencies in Hz
gyro            = 42.58;
par.f_wf        = [];
par.rel_amp     = [];
for i = 1:length(par.species)
    par.f_wf    = [par.f_wf; gyro * par.fieldStrength * par.species(i).frequency(:)];
    par.rel_amp = [par.rel_amp; par.species(i).relAmps(:)];
end

%% Use Berglunds region growing algorithm
fprintf('\t Use Berglunds RG algorithm \n')

f_wf_fmap = par.f_wf(2:end) * 2*pi;
rel_amp_fmap = par.rel_amp(2:end) / sum(par.rel_amp(2:end)); % normalize amplitudes

a = exp(complex(0,par.TE(1:3)' * f_wf_fmap')) * rel_amp_fmap;
A = [ones(3,1) a];

S = image; % size[nx,ny,nz,3]

if (~strcmp(class(S),'single'))
    S = single(S); % necessary for execution of c++ code
end

[bA, bB] = getPhasorCandidates(S,A);

mw = getMagnitudeWeight(S);
phasorexp = regionGrow(A,bA,bB,par.c1,par.c2,mw,par.voxelSize);

fieldmap_init = phasorexp/(2*pi*par.dTE);

