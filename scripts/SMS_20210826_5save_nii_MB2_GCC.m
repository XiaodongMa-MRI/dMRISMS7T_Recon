clear; clc;
addpath(genpath('./funcs/'))

datapath = '/home/naxos2-raid1/maxiao/Projects/dMRISMSRecon/Recon/data/Subj5_2021-618';
%% Para 
PEShift = 2;
MB = 2;
AP_direction = 'AP';
ACS = 'Gre';

if AP_direction(1) == 'A'
    rotation_90 = 90;
else
    rotation_90 = -90;
end

filepath = [datapath,'/PEShift',num2str(PEShift),'_',AP_direction,'/img_SMS_PF_SoS/'];
mkdir([datapath,'/PEShift',num2str(PEShift),'_',AP_direction,'/AResult/'])
savepath = [datapath,'/PEShift',num2str(PEShift),'_',AP_direction,'/AResult/'];

which_ghost_correction_method = 3;
GRAPPA_plus = '';%'sb',''

%% Show GRAPPA Img
IMG_GRAPPA = zeros(200,200,140,71);

for which_slice = 1:size(IMG_GRAPPA,3) %140
    disp(['Slice ',num2str(which_slice),' Loading....' ])
   
    Img1 = load([filepath,ACS,'_img_GRAPPA_PF_ACC_slice',GRAPPA_plus,'_',num2str(which_slice),...
    '_GhostMethod',num2str(which_ghost_correction_method),'_1']);
    Img2 = load([filepath,ACS,'_img_GRAPPA_PF_ACC_slice',GRAPPA_plus,'_',num2str(which_slice),...
    '_GhostMethod',num2str(which_ghost_correction_method),'_2']);
    Img3 = load([filepath,ACS,'_img_GRAPPA_PF_ACC_slice',GRAPPA_plus,'_',num2str(which_slice),...
    '_GhostMethod',num2str(which_ghost_correction_method),'_3']);

    IMG_GRAPPA(:,:,which_slice,1:30) = Img1.img;
    IMG_GRAPPA(:,:,which_slice,31:60) = Img2.img;
    IMG_GRAPPA(:,:,which_slice,61:71) = Img3.img(:,:,1:11);
end
IMG_GRAPPA = single(abs(IMG_GRAPPA));

img_name = ['MB',num2str(MB),'_',AP_direction,'_',ACS,'_img_GRAPPA_PF_ACC',GRAPPA_plus,'_GhostMethod',num2str(which_ghost_correction_method)];
save([savepath,img_name,'.mat'],'IMG_GRAPPA');
%% Save Nii
% img_name = ['MB',num2str(MB),'_',AP_direction,'_',ACS,'_img_GRAPPA_PF_ACC',GRAPPA_plus,'_GhostMethod',num2str(which_ghost_correction_method)];
% if AP_direction(1) == 'A'
% 	IMG_GRAPPA2 = abs(IMG_GRAPPA)*1e6;
% else
%     IMG_GRAPPA2 = abs(flip(flip(IMG_GRAPPA,1),2))*1e6;
% end
% IMG_GRAPPA2 = flipud(IMG_GRAPPA2);
% nii = make_nii(IMG_GRAPPA2,[1.05 1.05 1.05]); 
% save_nii(nii,[savepath,img_name,'.nii']);

