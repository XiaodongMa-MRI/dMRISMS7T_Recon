clear; clc;
addpath(genpath('funcs'))
addpath(genpath('scripts'))
%% Common parameters
filepath = '/home/range6-raid5/xpwu-data/dMRI_SMS_7T/20220722-ST001-Xiaodong_Nova_Denoising_LowSNR/';
savepath_root = '/home/range6-raid5/xpwu-data/dMRI_SMS_7T/20220722-ST001-Xiaodong_Nova_Denoising_LowSNR/';

filename = 'raw_fidsort_1mm_GreACS_PeShift2_AP_b2000_0p9mm';
filename_refs = 'raw_refs_1mm_GreACS_PeShift2_AP_b2000_0p9mm';
seq_name = 'b2000_0p9mm';

ACS = 'Gre'; % GRE
PEShift = 2;
AP_direction = 'AP';

savepath = [savepath_root,seq_name,'/PEShift',num2str(PEShift),'_',AP_direction,'/'];
if ~exist(savepath,'dir')
    mkdir(savepath)
end

MB = 2;
R_PE = 3; % acceleration factor
PF = 0.75;
par.kx_r = 234; % r:real
par.ky_r = 234; 
par.ky_s = 59; % s:sampled
par.nCH = 32;
par.nSL = 64;

%
%% Step1
% Para 1

% Para 2

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

par.nSL = 128;
par.nDYN = 1;
par.nSHOT = 1;

SMS_20210817_2get_single_band_MB2
%% Step3
savepath = savepath0;
DWI = true;  % for DWI 
par.nSL = 64;

SMS_20210817_3get_multi_band_MB2
%% Step4
which_ghostcorr_method = {'Ref','SVD','Ref-Ent'}; 

filepath = [savepath0,'single_band/'];
savepath = savepath0;
savepath2 = [savepath,'img_SMS_PF_SoS/'];
if ~exist(savepath2,'dir')
    mkdir(savepath2)
end
slice_gap = 64;
par.ky_s = 177;
par.ky_sr = 59;
par.nSL = 2;
par.nSL0 = 64;

% Recon Par
GRAPPA = 0; GRAPPA_sb = 1; SPIRIT = 0; 

ACS_name = {'Fleet','Gre'};

SMS_20210817_4GRAPPA_SPIRIT_sbmb_MB2_GCC
%% Step5
filepath = savepath2;
savepath = [savepath_root,seq_name,'/PEShift',num2str(PEShift),'_',AP_direction,'/AResult/'];
if ~exist(savepath,'dir')
    mkdir(savepath)
end
which_ghost_correction_method = 3;
GRAPPA_plus = 'sb';%'sb',''

% Show GRAPPA Img
IMG_GRAPPA = zeros(234,234,128,53);

for which_slice = 1:size(IMG_GRAPPA,3) %140
    disp(['Slice ',num2str(which_slice),' Loading....' ])
   
    Img1 = load([filepath,ACS,'_img_GRAPPA_PF_ACC_slice',GRAPPA_plus,'_',num2str(which_slice),...
    '_GhostMethod',num2str(which_ghost_correction_method),'_1']);

    IMG_GRAPPA(:,:,which_slice,:) = Img1.img;
end
IMG_GRAPPA = single((IMG_GRAPPA));

img_name = ['MB',num2str(MB),'_',AP_direction,'_',ACS,'_img_GRAPPA_PF_ACC',GRAPPA_plus,'_GhostMethod',num2str(which_ghost_correction_method)];
save([savepath,img_name,'.mat'],'IMG_GRAPPA','-v7.3');