%% Demonstration code for Dixon-RAVE
% See readme.txt for more details

clear;
clc;
close all;

fprintf('Dixon-RAVE demo reconstruction \n \n')
addpath(genpath(pwd));

%% Load measurement
fprintf('\t Loading sample data \n')
% load 'data_liver.mat' % Load liver data
load 'data_breast.mat' % Load breast data

%% Set optimization parameters
par.FD1Weight           = 0.00; % TV regularization parameter for fat and water

par.iter_gn             = 3;
par.iter_cg             = 3;
par.prec                = 'single'; % Precision for fitting ('single'/'double')
par.threshold           = 0.4; % stopping criterion

par.verbose             = 1;

%% Set fat/water/field map parameters
par.species(1).frequency = 0;
par.species(1).relAmps   = 1;

% 6-peak model
par.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
par.species(2).relAmps   = [0.087 0.693 0.128 0.004 0.039 0.048];

par.estfieldmap          = 1; % Perform fieldmap estimation for initialization
par.c1                   = 0.75; % Magnitude weight threshold for seed points
par.c2                   = 0.25; % Determines number of used seeds
par.smoothfm             = 1; % Smooth fieldmap after estimation
par.smoothfmiter         = 10;

par.timemap.usetimemap   = 1; % Use time map instead of constant TEs
par.coilmaps             = coilmaps; % Save coilmaps in par

%% Undersample dataset retrospectively
%{
Nprojcut            = 64;
fprintf('\t Use only first %d projections \n',Nprojcut)
kdata               = kdata(:,1:Nprojcut,:,:);
par.Nproj           = Nprojcut;
%}

%% Perform modelbased reconstruction
par.recotime            = tic;
out                     = modelbased_fw(par,kdata);

water                   = out.water;
fat                     = out.fat;
fieldmap                = out.fieldmap;
images                  = out.images;
optinfo                 = out.optinfo; % convergence information

fprintf('\t Optimization completed \n')
par.recotime            = toc(par.recotime);
fprintf('\t Time for optimization: %.2f min\n',par.recotime/60)

%% Show results
water_disp = imresize(water,[size(water,1)*4 size(water,2)*4],'bilinear');
water_disp = water_disp./max(abs(water_disp(:)));
fat_disp = imresize(fat,[size(fat,1)*4 size(fat,2)*4],'bilinear');
fat_disp = fat_disp./max(abs(fat_disp(:)));

figure(1)
imshow([abs(water_disp) abs(fat_disp)],[],'Border','tight'); title('Water');
