function [data] = yarra_GRASP_xdme(work_path, meas_file, output_path, ~, mode_file, i_echo, echoparams)

% RECONSTRUCT ONE ECHO. SELECT ECHO THROUGH INPUT PARAMETERS!
% WATERFIETS: ALSO CORRECT FOR MOTION WITH XD-GRASP.

%ivom: select the echo index to reconstruct. Echoes sorted in
%increasing time!
    
    traj_flip    = [0 0];
    blipAngleDeg = 2;
    nresp        = 4;
  
    version='0.00611_S0L3R0';

    fprintf('\n');
    fprintf('Yarra Basic GRASP Recon %s\n',version);
    fprintf('----------------------------\n');
    fprintf('\n');
    
    if verLessThan('matlab','9.0')
        disp('ERROR: This module requires Matlab 2016a or newer.');
        return;
   end

    pct = ver('distcomp');
    pct_available = true;
   
    if isempty(pct)
        disp('Parallel Computing Toolbox not available.');
        disp('Reconstruction will run slowly.');
        pct_available=false;
    end    
    
    %reconStart=tic;    
    addpath(genpath('code'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reconstruction parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    rawfile=[work_path filesep meas_file];
    
    % Set default values for the reconstruction parameters
    recoparams=struct;
    recoparams.inner_iterations    = 5;
    recoparams.outer_iterations    = 3;
    recoparams.tv_lambda           = 0.02;
    recoparams.spokes_per_frame    = 73%144;
    recoparams.sliceFrom           = 0;
    recoparams.sliceTo             = 0;    
    recoparams.cores_to_use        = feature('numcores');
    recoparams.apply_kx_filter     = true;
    recoparams.compressed_channels = 10;
    %recoparams.nresp               = 4;

    % Read the module settings from the mode file
    settings=yarra_read_mode_section(mode_file, 'GRASP');
    
    % Overwrite the defaults with settings from the mode file
    if (isfield(settings,'spokesperframe'))
        recoparams.spokes_per_frame=settings.spokesperframe;
    end
    if (isfield(settings,'lambda'))
        recoparams.tv_lambda=settings.lambda;
    end
    if (isfield(settings,'slicefrom'))
        recoparams.sliceFrom=settings.slicefrom;
    end    
    if (isfield(settings,'sliceto'))
        recoparams.sliceFrom=settings.sliceto;
    end
    if (isfield(settings,'sliceto'))
        recoparams.sliceFrom=settings.sliceto;
    end    
    if (isfield(settings,'maxthreads'))
        recoparams.cores_to_use=settings.maxthreads;
    end        
    if (isfield(settings,'applykzfilter'))
        recoparams.apply_kx_filter=settings.applykzfilter;
    end        
    if (isfield(settings,'compressedchannels'))
        recoparams.recoparams.compressed_channels=settings.compressedchannels;
    end      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read rawdata 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Reading file %s...\n',rawfile);       
    
    twixfile=mapVBVD(rawfile);
    if length(twixfile)>1
        twixfile = twixfile{end};
    end
    
    centerpar    =max(twixfile.image.centerPar); %center of 
    partitions   =twixfile.image.NPar; %measured partitions with OS w/o PF
    imagesPerSlab=twixfile.hdr.Meas.lImagesPerSlab; %final number of reconstructed images. 
    nparPF       =twixfile.hdr.Meas.Partitions; %partitions with OS with PF
    TA           =twixfile.hdr.Meas.lTotalScanTimeSec;
    
    disp('Sorting data...');
    
    % Get the kspace data and convert to double now to avoid repetitive
    % recasting into at later points
    kdata = double(twixfile.image{''});
    
    [nx, nc, ntviews, ~, ne] = size(kdata);
    %kdata as [nx,nc,ntviews,npar]
    
    
    if i_echo > ne
        error('Echo index exceeds the number of echoes of dataset.');
    end
    
    kdata = squeeze(kdata(:, :, :, :, i_echo));
    %Res_Signal = motionDetGrasp2(kdata, twixfile);
    Res_Signal = motionDetGrasp(kdata, 0, 0);

    % Free up the data again
    clear twixfile;    


    kdata=permute(kdata,[1 3 4 2]);
    kdata=zerofill_Kz(kdata,3,partitions,imagesPerSlab,centerpar,recoparams.apply_kx_filter);

    
    disp('Finished reading data');          
    disp('Performing Partial Fourier...');
    %input as: [nx, nc, ntviews, nz] 
    %kdata = pfGoldenAngle0611(kdata, nparPF);
    %output same as input: permute. 
 
    Traj_ME = Trajectory_GoldenAngle_ME(ntviews, nx, ne, blipAngleDeg, 'flip', traj_flip, 'gradcorr', 1); 
    Traj = Traj_ME(:, :, i_echo);
    clear Traj_ME
    %disp('Zero-padding...')
    %kdata = permute(kdata, [1 3 4 2]);
    %kdata=zerofill_Kz(kdata,3,partitions,imagesPerSlab,centerpar,recoparams.apply_kx_filter);

    [nx,ntviews,nz,nc]=size(kdata);
    
    DensityComp = abs(Traj);
    DensityComp = DensityComp/(max(DensityComp(:)));
    
    %[Traj,DensityComp]=Trajectory_GoldenAngle(ntviews,nx);

    % FFT along the kz dimension
    kdata1=fftshift(ifft(fftshift(kdata,3),[],3),3); 
    clear kdata
    
    figure()
    imagesc(squeeze(abs(kdata1(256, :, :, 18)))');
    colormap('gray');
    title('Projection before coil compression')
    hold on
    plot(Res_Signal*40+40, 'r');
    
    %nline=recoparams.spokes_per_frame;  % Number of spokes per frame
    nline=floor(ntviews/nresp);            % Number of frames   
    nt = nresp;
    
    disp(' ');
    disp('== Information ================================');
    fprintf('Number of slices = %d\n',   imagesPerSlab);
    fprintf('Frames           = %d\n',   nt);
    fprintf('Spokes per frame = %d\n',   nline);
    fprintf('Time per frame   = %f s\n', TA/ntviews*nline);    
    fprintf('TV lambda        = %f\n',   recoparams.tv_lambda);
    disp('===============================================');
    disp(' ');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coil compression
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    if (recoparams.compressed_channels > 0) && (nc>recoparams.compressed_channels)
        fprintf('Performing coil compression (to %d channels)...\n',recoparams.compressed_channels);
        ncc=recoparams.compressed_channels;
        D=reshape(kdata1,nx*nz*ntviews,nc);
        [~,~,V]=svd(D,'econ');
        kdata1=reshape(D*V(:,1:ncc),nx,ntviews,nz,ncc);
        [nx,ntviews,nz,nc]=size(kdata1);
    end
   
    figure()
    imagesc(squeeze(abs(kdata1(256, :, :, 5)))');
    colormap('gray');
    title('Projection after coil compression')
    hold on
    plot(Res_Signal*40+40, 'r');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coil sensitivity calculation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    disp('Estimating coil profiles...');
    
    ref=ReconForCoilMaps(kdata1,Traj,DensityComp, nz,0);
    b1=adapt_array_3d(ref);
    b1=b1/max(abs(b1(:)));
           
    %IVOM%IVOM%IVOM%IVOM
    % SORT FOR MOTION CORR. 
    %IVOM%IVOMAATMAN$IVOM
    
    [~, index] = sort(Res_Signal, 'descend'); 
    Traj = Traj(:, index);
    DensityComp = DensityComp(:, index);
    kdata1 = kdata1(:, index, :, :); 
    
% % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRASP reconstruction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    disp('Preparing GRASP...');
    fprintf('correct version...');
    
    % Start the thread pool for parallel processing
    if pct_available
        poolobj=parpool(recoparams.cores_to_use);
        fprintf('Using pool with %d threads.\n',poolobj.NumWorkers);
    end       
    
    % Sort the data into time frames with desired number of spokes
    clear kdata_nt Traj_nt DensityComp_nt;
    for ii=1:nt
        kdata_nt(:,:,:,:,ii)=kdata1(:,(ii-1)*nline+1:ii*nline,:,:);
        Traj_nt(:,:,ii)=Traj(:,(ii-1)*nline+1:ii*nline,:);
        DensityComp_nt(:,:,ii)=DensityComp(:,(ii-1)*nline+1:ii*nline);
    end
    
    clear kdata1
    
    % Apply the density compensation function to the k-space data
    kdata_nt=kdata_nt.*repmat(sqrt(permute(DensityComp_nt,[1 2 4,5,3])),[1,1,nz,nc,1]);
    nt=size(kdata_nt,5);
           
    % Preinitialize the results array for multithreading use
    data=zeros(nx/2,nx/2,nz,nt);

    % Determine loop limit for the multi-slice reconstruction.
    % Limiting the number of spokes can be helpful when testing the recon.
    minSlice=1;
    maxSlice=nz;    
    if (recoparams.sliceFrom>0)
        disp('Limiting number of slices...');
        minSlice=min(recoparams.sliceFrom,nz);
        maxSlice=min(recoparams.sliceTo,nz);
    end

    % Calculate gridding solutions as initial guess and to properly scale 
    % the TV weight
    initialEstimate=zeros(nx/2,nx/2,nz,nt);
    parfor zz = 1:nz 
        tempparam=[];
        tempparam.y=double(kdata_nt(:,:,zz,:,:));
        tempparam.SG=1;
        tempparam.E=MCNUFFT3D(double(Traj_nt),double(DensityComp_nt),double(b1(:,:,zz,:)));

        initialEstimate(:,:,zz,:)=tempparam.E'*tempparam.y;
        tempparam=[];
    end    
    globalTVWeight=max(abs(initialEstimate(:)))*recoparams.tv_lambda;
    
    % Load into local variables to avoid communication overhead
    innerIterations=recoparams.inner_iterations;
    outerIterations=recoparams.outer_iterations;    
    
    disp('Running GRASP reconstruction now...');
        
    %fullstart=tic;
    parfor zz = minSlice:maxSlice %zz = minSlice:maxSlice % loop through all slices    
        fprintf('Calculating slice %d...\n', zz);
        param=[];
        
        % No softweighting is used here
        param.SG=1; 
        
        % Limit data for MCNUFFT3D operator to one slice for enabling
        % multithreaded calculation of slices.
        param.y=double(kdata_nt(:,:,zz,:,:));
        
        % Prepare the operator for the forward operation
        param.E=MCNUFFT3D(double(Traj_nt),double(DensityComp_nt),double(b1(:,:,zz,:)));
        param.TV=TV_Temp3D;
        param.TVWeight=globalTVWeight;
        param.TVWeightRes=0;
        param.L1Weight=0;
        param.nite=innerIterations;
        param.display=1;
        param.slice=zz;

        recon_cs=initialEstimate(:,:,zz,:);
        
        for n=1:outerIterations
            recon_cs = CSL1NlCg_4DRadial_SG(recon_cs,param);
        end
        data(:,:,zz,:)=abs(recon_cs);    

        param=[];
    end
    
    % Normalize the result
    data=data/max(data(:));

    %toc(fullstart);

    %Shutdown the thread pool
    if pct_available
        delete(gcp('nocreate'));
    end    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Storage of images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %disp('Saving the results...');
    
    %data=flip(data,3); % Flip slice order
    %yarra_dicom_timeseries(data,output_path);

    %disp('Complete duration: ');
    %toc(reconStart);
    