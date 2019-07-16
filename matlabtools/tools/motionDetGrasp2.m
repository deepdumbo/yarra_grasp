function resSignal = motionDetGrasp2(rawdata, twix)
% Function for extracting respiratory motion signal from GRASP data.
% It contains several methods & options for performing this task, for
% optimization purposes. The effect of each of these on respiratory motion
% detection still needs to be evaluated.
% The option doAlign enables the alignment of k-space spokes prior to taking 
% the Fourier Transform along Z. This is because the location of the peak
% intensity along each spoke can vary between spokes; i.e., the max is not
% always observed at nx/2. Alignment compensates for these variations.
% The option doNormalize controls whether to normalize each profile in each 
% coil element to unit sum along z, possibly reducing profile-to-profile
% variations in total signal, e.g. due to cardiac pulsation
% The option doSmoothing controls whether temporal smoothing should be
% applied. There are two methods for this, controlled by the parameter
% smthMethod:
% 1: temporal smoothing, with a smoothing kernel type and width defined by 
% smthType and smthSpan
% 2: temporal low-pass filtering with a cutoff frequency (in Hz) defined by
% the parameter cutoff. This part uses the twix data to determine the
% sampling frequency of the temporal profiles (i.e. TR times number of kz
% samples per spoke)
% 
% Input parameters:
% - rawdata: contains the raw k-space data. Must have dimensions
% [nx,nc,ntviews,nkz], which is the same order they are read in by mapVBVD.
% - twix: twix header information
% 
% Written by Marnix Maas (Marnix.Maas@radboudumc.nl), April 2019

%% Set parameters
% Switches
doAlign         = 1;
doNormalize     = 1;
doSmoothing     = 1;

% Other parameters
smthMethod      = 1;                % Switch for selecting method for smoothing signals: 0 = no smoothing; 1 = smoothing; 2 = low-pass filtering
smthType        = 'lowess';         % Type of smoothing filter
smthSpan        = 7;                % Span of smoothing filter
cutoff          = 0.5;              % Cutoff frequency of low-pass filter in Hz
firOrder        = 16;               % Order of FIR low-pass filter
rfLim           = [0.1 0.5];        % Limits of 'likely' respiratory frequencies in Hz
nPCs            = 10;               % Number of principal components to analyze for respiratory frequency components

%% Get data dimensions
[nx, nc, ntviews, nkz] = size(rawdata);


%% Define time and frequency axes
% Time
TR              = twix.hdr.Config.TR;
nParMeas        = twix.hdr.Config.NParMeas;
dt              = 1E-6 * TR  * nParMeas;
t               = dt * (0:ntviews-1);

% Frequency
df              = 1/max(t);
f               = df*(0:length(t)-1);


%% Get central k-space lines along kz
if doAlign
    % Choosing the same x-coordinate on each spoke results in modulations,
    % because the peak signal intensity along x is not always at the same
    % point. This may be due to eddy currents, which lead to small shifts in
    % the echo location.
    % Therefore, the echoes should be 'aligned' first. We do this by selecting
    % the point along x for which the sum along z is maximum.
    % First, concatenate z and coil dimensions
    rawdata     = permute(rawdata,[1,3,4,2]);           % Dimensions now [x, views, partitions, coils]
    rawdata     = reshape(rawdata,[nx, ntviews, nkz*nc]);
    
    % Then, calculate sum of absolute values along z*c for each point x in each spoke
    sumzc       = sum(abs(rawdata),3);
    
    % Find the index of the maximum in each spoke
    [~, idx]    = max(sumzc,[],1);
    
    % Select data at this index to create the profile:
    % - Reshape raw data into nx x ntviews*nkz*nc matrix
    % - Replicate indices vector npar*nc times (same x-point should be 
    % selected in each partition and coil for a given view)
    % - Calculate linear indices of points to be used
    % - Select these points
    % - Reshape back to views x partitions x coils
    rawdata     = reshape(rawdata,[nx, ntviews*nkz*nc]);
    idx         = repmat(idx, [1, nkz*nc]);
    idx2        = sub2ind(size(rawdata), idx, 1:size(rawdata,2));
    kzProf      = rawdata(idx2);
    kzProf      = reshape(kzProf,[ntviews, nkz, nc]);
else
    kzProf      = squeeze(rawdata(nx/2+2,:,:,:));           % MCM: This was originally used in the NYU demo. Why nx/2 + 2, and not +1?
end


%% Create z-profiles (FFT along z)
zProf           = abs(fftshift(ifft(kzProf,[],2),2));     % FFT along the kz dimension
nz              = size(zProf,2);

% Note that NYU demo code used interpolation here, ie to 400 points per
% profile. It is unclear what the advantage of this would be, other than
% for plotting purposes.

%% Normalize z-profiles
if doNormalize
    % Normalize each profile in each coil element to unit sum along z:
    % This should eliminate profile-to-profile variations in total signal,
    % which could e.g. arise due to cardiac pulsation.
    % Perhaps there are other reasons for this as well, but could this occur
    % because of *respiratory* motion? 
    
    % Calculate sums
    sumz            = squeeze(sum(zProf,2));
    
    % Divide each profile by its sum
    sumz            = reshape(sumz,[1, ntviews*nc]);
    zProf           = permute(zProf,[2 1 3]);
    zProf           = reshape(zProf,[nz,ntviews*nc]);
    zProf           = bsxfun(@rdivide, zProf, sumz);
    
    % Reshape back to original dimensions    
    zProf           = reshape(zProf,[nz,ntviews,nc]);
    zProf           = permute(zProf,[2 1 3]);
end

%% Perform temporal smoothing

% Reshape profiles
zProf               = reshape(zProf,[ntviews, nz*nc]);
    
if doSmoothing
    % Temporal smoothing: ensure that only frequency components in the range of
    % respiratory motion remain. E.g. using low-pass filtering, smoothing or
    % spline fitting    
    switch smthMethod
        case 0
            % No smoothing or filtering
            zProfSmooth = zProf;
        case 1
            % Temporal smoothing
            zProfSmooth = zeros(size(zProf));
            for i=1:size(zProf,2)
                zProfSmooth(:,i) = smooth(zProf(:,i),smthSpan,smthType);
            end
            
        case 2
            % Temporal filtering
            fs              = 1/dt;                                 % Sampling frequency
            normCutoff      = cutoff/fs;                            % Normalized cutoff frequency
            lpf             = fir1(firOrder,normCutoff,'low');      % Low pass filter
            zProfSmooth     = filter(lpf,1,zProf,[],1);
            zProfSmooth     = [zProfSmooth(ceil(length(lpf)/2):end,:);
                               zeros(floor(length(lpf)/2), size(zProfSmooth,2))];
    end
    
end

%% Perform PCA to maximize information content
% Create linear combination of z,c-vs-t profiles that maximizes the
% variation along t, using PCA.
cmat                = cov(zProfSmooth);
[REV, lambda]       = eig(cmat);                    % REV = right eigenvectors; lambda = eigenvalues
lambda              = diag(lambda);
[~, rindices]       = sort(-1*lambda);
% lambda              = lambda(rindices);
REV                 = REV(:,rindices);
pcProfiles          = (REV' * zProfSmooth')';

% The first few principal components contain contributions from both
% respiratory and bulk motion. We want to select the one
% with the highest contribution from resp motion. We do this by selecting
% the one with the highest amplitude in the corresponding frequency range.
% - create freq axis & find indices where frequencies are in resp frequency
% range
[~, rfind]      = find(f>=rfLim(1) & f<rfLim(2));

% Fourier transform first few PCs along time
specProf        = fft(pcProfiles(:,1:nPCs),[],1);
% Calculate maxima within the resp frequency range
rfMax           = max(abs(specProf(rfind,:)),[],1);
% Get index of highest of these maxima
[~, rfMaxInd]   = max(rfMax);

% We use the corresponding principal component as respiratory signal
fprintf('Using PC %d\n', rfMaxInd);
resSignal       = pcProfiles(:,rfMaxInd);

%% Perform baseline correction
% To be implemented