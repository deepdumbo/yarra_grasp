function [out_img_TVxy, out_img_mcgridding, t] = reconMeGrasp(kdata)
% Function for performing a reconstruction of multi-echo golden-angle 
% radial sparse (GRASP) k-space data using compressed sensing. 
% NOTE: it is assumed that a fourier transform along the slice direction
% has already been performed!
% Input parameters:
% kdata: contains the raw k-space data. Must have dimensions
% [nx,nc,ntviews,nz], which is the same order they are read in by mapVBVD.
% Res-Signal: respiratory motion signal (estimated in a previous step)
% 
% Adapted from Demo_GRASP_DCE and Demo_XDGRASP_NonContrast 
% (NYU Demo provided by Li Feng) 
% and mostly from the demo file 'iterative_multicoil.m' by Tobias Block
% by Camille Van Speybroeck (Camille.vanSpeybroeck@radboudumc.nl), 
% August 2018
    
% Get data dimensions.
[nx,nc,ntviews,nz,ne]=size(kdata);

% Generate trajectory.
Traj=Trajectory_GoldenAngle_ME(ntviews, nx, ne, 2, 'normalize',1,'flip',[0 0]);

% Calculate density compensation.
dcf = dcfGridding(ntviews, nx);

param.dcf = dcf;
dcf = repmat(dcf, [1,1,3]);

% Set NUFFT parameters.
fftSize=[nx,nx];


% Create empty images to sum the reconstructed channels (for sum-of-squares
% combination).
slices = 1:nz;
out_img_cs   = single(zeros(nx,nx,length(slices),ne));
out_img_TVxy = out_img_cs;
out_img_mcgridding = out_img_cs;

% Set parameters.
param.nvar=2*nx^2; 
param.br=nx;
param.rawdata=kdata;
param.spokes=ntviews;
param.channels=nc;

% Coil profile parameters.
param.lambdaCoils=10;       %iterativeradial version = 10
param.stopTolCoils=1e-10;   %iterativeradial version = 1e-5; 
                            % This seems like a good value
param.iterationsCoils=1000; %iterativeradial version = 1000

%Image parameters.
%param.lambdaImage=2e-7;    %iterativeradial version = 0.00002 (=2e-5)
%This parameter is obtained by evaluating the outcome without penalty terms.
param.iterationsImage=60;   %iterativeradial version = 60
param.stopTolImage=1e-11;   %iterativeradial version = 4e-8
                            % This seems like a good value

%lbfgs parameters
param.MaxFuncEvals=10000;   %poblano version = 10000
param.RelFuncTol=1e-6;      %poblano version = 1e-6
param.LineSearch_maxfev=20; %poblano version = 20


% Global counter for displaying the number of evaluations
global iterationCounter
iterationCounter=0;

%%% Motion resolved gridding and CS-recon using MCNUFFT.
%spokesOffset = 0;
b1 = zeros(nx, nx, nc, ne);
%load('b1.mat');            To skip the coil profile estimation.
recon_cs = zeros(nx, nx, ne);
temp_recon_cs = zeros(nx, nx);
t = [];
tic;

    
for sl=slices
    % Select slice to be reconstructed.
    kdata1=squeeze(kdata(:,:,:,sl,:));
    for ec=1:ne
        FT = NUFFT(Traj(:,:,ec), 1, 1, [0,0], fftSize, 2);
        param.FT = FT;

        fprintf('Starting GRASP recon of echo %d, slice %d / %d\n', ec, sl, length(slices));

        % ## Step 1: Estimate coil profiles.
        % Loop over the channels to create "smooth" image for each coil.
        for ic=1:nc
            % Read k-space data for the channel.
            param.y=double(squeeze(kdata1(:,ic,:,ec)));

            % Initialize optimizer with empty image.
            x0=zeros(param.nvar,1);
            iterationCounter=0;

            % Set stop criterion depending on the value range of the raw data.
            options.StopTol=param.stopTolCoils;
            options.RelFuncTol=param.RelFuncTol;
            options.MaxFuncEvals=param.MaxFuncEvals;
            options.LineSearch_maxfev = param.LineSearch_maxfev;

            % Run the optimizer for 10 iteration without penalty terms.
            param.enablePenalty=0;
            options.MaxIters=10;
            out = lbfgs(@(x) costfunction_coils(x,param), x0, options);

            % Now enable the penalty terms and run the optimizer for the remaining
            % iterations.
            param.enablePenalty=1;
            options.MaxIters=param.iterationsCoils-10; 
            out = lbfgs(@(x) costfunction_coils(x,param), out.X, options);

            % Reshape the result vector into a 2D image and store it.
            b1(:,:,ic,ec)=vec_to_img(out.X,param.br);
        end

        % Sum-of-squares calculation of the coil profiles.
        ssqcoil=zeros(nx,nx);
        for ic=1:nc
            ssqcoil=ssqcoil+abs(b1(:,:,ic,ec).*b1(:,:,ic,ec));
        end

        ssqcoil=sqrt(ssqcoil);
        for ic=1:nc
            b1(:,:,ic,ec)=b1(:,:,ic,ec)./ssqcoil;
        end
       
   
        % ## Step 2: TV in xy dimension.
        % Initialize the optimizer.
        param.RelFuncTol=1e-20;
        param.coilprofile=b1;
        param.y=double(squeeze(kdata1(:,:,:,ec)));
        x0=zeros(param.nvar,1);
        iterationCounter=0;
        options.StopTol=param.stopTolImage;
        options.RelFuncTol=param.RelFuncTol;
        options.MaxFuncEvals=param.MaxFuncEvals;
        options.LineSearch_maxfev = param.LineSearch_maxfev;
        
        
        % First run the optimizer for 5 iterations without penalty terms.
        param.enablePenalty=0;
        options.MaxIters=5;
        out = lbfgs(@(x) costfunction_image(x,param), x0, options);
        
        temp_recon_cs(:,:)=vec_to_img(out.X,param.br);
        param.lambdaImage=mean(abs(temp_recon_cs(:)))*0.0008;       % Determine parameter with outcome without penalty terms
%        fprintf('TV_xy weightingfactor for echo %d = %g\n', ec, param.lambdaImage);
        
        % Now enable the penalty terms and run the optimizer for the remaining
        % iterations.
        param.enablePenalty=1;
        options.MaxIters=param.iterationsImage-5;
        out = lbfgs(@(x) costfunction_image(x,param), out.X, options);
        fprintf(1, 'Repeat... ExitFlag = %d; %s\n', out.ExitFlag, out.ExitDescription);
        
        % Reshape the result vector into 2D image format.
        recon_cs(:,:,ec)=vec_to_img(out.X,param.br);

    end  %for ec=1:ne
    
%     save('/nis_home/mrs/camilles/Software/ImageReconstruction/ImageReconstruction/b1.mat', 'b1');
    
    % ## Step 3: TV in echo dimension.
    % Perform multi-coil Gridding recon as starting point for CS.

    param.E=MCNUFFT(double(Traj),double(dcf),double(b1(:,:,:,1)));  % MCNUFFT does not seem to work well.
    kdata1 = permute(kdata1, [1,3,2,4]);
    param.y=double(kdata1);
%    recon_cs = param.E'*param.y;
    out_img_mcgridding(:,:,sl,:)=param.E'*param.y;
    out_img_TVxy(:,:,sl,:)=recon_cs;                                % Gives very smooth outcome

    % Perform compressed sensing recon in extra dimension (e.g. echo)
    % TV in echo dimension does not work yet
%    param.TV_dim1=TV_Temp;
%    param.TVWeight_dim1=max(abs(recon_cs(:)))*0.0; %Was 0.02
%    param.TVWeight_dim2=0;
%    param.nite = 4; %Was 4
%    param.display=1;
    %     
    % tic
%    for n=1:2
%      recon_cs = CSL1NlCg_XDGRASP_Mx(recon_cs,param);
%    end
    % time=toc;
    % time=time/60
%    out_img_cs(:,:,sl,:)=recon_cs;
%    t(sl) = toc;
end
    
% Normalize images
out_img_mcgridding = out_img_mcgridding/max(abs(out_img_mcgridding(:)));
out_img_TVxy       = out_img_TVxy/max(abs(out_img_TVxy(:)));
%out_img_cs         = out_img_cs/max(abs(out_img_cs(:)));           TV in echo dimension does not work yet


figure, imshow3(abs(out_img_mcgridding(:,:,:,:)),[0 .5]);title('MCGrdding');

figure, imshow3(abs(out_img_TVxy(:,:,:,:)),[0 .5]);title('TVxy');   % This image is very smooth

%figure, imshow3(abs(out_img_cs(:,:,:,:)),[0 .5]);title('CS');      TV in echo dimension does not work yet


disp('...done.');