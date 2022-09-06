% This is the main script for reconstruction of 7T SMS dMRI data.
% Original scripts were authored by Ziyi Pan.
% Modified by Xiaodong Ma: combing the main scripts of 5 steps into one
% script.
%%
clear; clc;
%% add utils
addpath(genpath('funcs'))
addpath(genpath('scripts'))
%% Acq and Recon parameters. Specify them mannually!!!
% filepath: the folder containing the raw data files
% savepath_root: any folder you want to save the reconstruction results
filepath = '/home/range6-raid5/xpwu-data/dMRI_SMS_7T/20220722-ST001-Xiaodong_Nova_Denoising_LowSNR/';
savepath_root = '/home/range6-raid5/xpwu-data/dMRI_SMS_7T/20220722-ST001-Xiaodong_Nova_Denoising_LowSNR/';

% file names for multi-band imaging and refs raw data; no extension.
filename = 'raw_fidsort_1mm_GreACS_PeShift2_AP_b2000_0p9mm';
filename_refs = 'raw_refs_1mm_GreACS_PeShift2_AP_b2000_0p9mm';

% name the folder where the reconstruction results will be saved.
seq_name = 'b2000_0p9mm';

% Acq parameters
ACS = 'Gre'; % 'GRE' or 'Fleet'
PEShift = 2; % FOV shift; 2 if shift=1/2
AP_direction = 'AP'; %'AP' or 'PA'

MB = 2; % MB factor
R_PE = 3; % in-plane acceleration factor
PF = 0.75; % partial fourier factor
par.kx_r = 234; % readout dimension without acceleration (r:real)
par.ky_r = 234; % PE dimension without acceleration 
par.ky_s = 59; % readout dimension after acceleration: round(ky_r/R_PE*PF) (s:sampled)
par.nSL = 64; % number of slices acquired (after MB acceleration)
nVol = 53; % number of volumes (or diffusion directions)

par.nCH = 32; % number of receive channels
par.nDYN = 1; % number of dynamics (usually 1)
par.nSHOT = 1;% number of shots (usually 1)

% Recon Par
GRAPPA = 0; GRAPPA_sb = 1; SPIRIT = 0; % skip the recon method by set to 0

%% Create the folder for saving the results. No need to specify mannually.
savepath = [savepath_root,seq_name,'/PEShift',num2str(PEShift),'_',AP_direction,'/'];
if ~exist(savepath,'dir')
    mkdir(savepath)
end

%% Step1: load ACS data
SMS_20210817_1get_ACS_MB2

%% Step2
if AP_direction(1) == 'A' % for display only
    rotation_90 = 90;
else
    rotation_90 = -90;
end
savepath0 = savepath;
savepath = [savepath0,'single_band/'];
if ~exist(savepath,'dir')
    mkdir(savepath)
end

par.nSL = par.nSL * MB;

SMS_20210817_2get_single_band_MB2

%% Step3: load imaging (multi-band) data
savepath = savepath0;
DWI = true;  % for DWI 
par.nSL = par.nSL / MB;

SMS_20210817_3get_multi_band_MB2

%% Step4: reconstruct images using specified 3D GRAPPA or SPIRiT;
which_ghostcorr_method = {'Ref','SVD','Ref-Ent'}; 

filepath = [savepath0,'single_band/'];
savepath = savepath0;
savepath2 = [savepath,'img_SMS_PF_SoS/'];
if ~exist(savepath2,'dir')
    mkdir(savepath2)
end

% switching of variable names is confusing but should be correct
slice_gap = par.nSL; 
par.ky_sr = par.ky_s; 
par.ky_s = par.ky_sr * R_PE;
par.nSL0 = par.nSL;
par.nSL = MB; 

ACS_name = {'Fleet','Gre'};

SMS_20210817_4GRAPPA_SPIRIT_sbmb_MB2_GCC

%% Step5: re-organize and combine reconstructed images into one mat file 
filepath = savepath2;
savepath = [savepath_root,seq_name,'/PEShift',num2str(PEShift),'_',AP_direction,'/AResult/'];
if ~exist(savepath,'dir')
    mkdir(savepath)
end
which_ghost_correction_method = 3;
GRAPPA_plus = 'sb';%'sb',''

% Show GRAPPA Img
IMG_GRAPPA = zeros(par.kx_r,par.ky_r,par.nSL0*MB,nVol);

for which_slice = 1:size(IMG_GRAPPA,3) %140
    disp(['Slice ',num2str(which_slice),' Loading....' ])
   
    Img1 = load([filepath,ACS,'_img_GRAPPA_PF_ACC_slice',GRAPPA_plus,'_',num2str(which_slice),...
    '_GhostMethod',num2str(which_ghost_correction_method),'_1']);

    IMG_GRAPPA(:,:,which_slice,:) = Img1.img;
end
IMG_GRAPPA = single((IMG_GRAPPA));

img_name = ['MB',num2str(MB),'_',AP_direction,'_',ACS,'_img_GRAPPA_PF_ACC',GRAPPA_plus,'_GhostMethod',num2str(which_ghost_correction_method)];
save([savepath,img_name,'.mat'],'IMG_GRAPPA','-v7.3');

%% delete intermediate files (NOT TESTED YET!!!)
% step1
delete([savepath0, 'ksp_*']);
% step2
rmdir([savepath0,'single_band'],'s');
% step3
delete([savepath0, 'ksp_*']);
% step4
rmdir([savepath2,'single_band'],'s');


