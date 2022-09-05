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

%% Show GRAPPA Img
IMG_SPIRiT = zeros(200,200,140,71);
Img1 = zeros(200,200,30);
Img2 = zeros(200,200,30);
Img3 = zeros(200,200,11);

for which_slice = 1:size(IMG_SPIRiT,3) %140
    disp(['Slice ',num2str(which_slice),' Loading....' ])
    for nbi = 1:size(Img1,3)   
        Imgtmp = load([filepath,ACS,'_img_L1-SPIRiT_PF_ACC_',num2str(which_slice),'_nb',num2str(nbi)...
        '_GhostMethod',num2str(which_ghost_correction_method),'_1_Lambda5e-08']);
        Img1(:,:,nbi) = Imgtmp.img;
    end
    for nbi = 1:size(Img2,3)   
        Imgtmp = load([filepath,ACS,'_img_L1-SPIRiT_PF_ACC_',num2str(which_slice),'_nb',num2str(nbi)...
        '_GhostMethod',num2str(which_ghost_correction_method),'_2_Lambda5e-08']);
        Img2(:,:,nbi) = Imgtmp.img;
    end
    for nbi = 1:size(Img3,3)   
        Imgtmp = load([filepath,ACS,'_img_L1-SPIRiT_PF_ACC_',num2str(which_slice),'_nb',num2str(nbi)...
        '_GhostMethod',num2str(which_ghost_correction_method),'_3_Lambda5e-08']);
        Img3(:,:,nbi) = Imgtmp.img;
    end
    IMG_SPIRiT(:,:,which_slice,1:30) = Img1;
    IMG_SPIRiT(:,:,which_slice,31:60) = Img2;
    IMG_SPIRiT(:,:,which_slice,61:71) = Img3;
end
IMG_SPIRiT = single(abs(IMG_SPIRiT));
%% Save Nii
img_name = ['MB',num2str(MB),'_',AP_direction,'_',ACS,'_img_L1-SPIRiT_PF_ACC','_GhostMethod',num2str(which_ghost_correction_method)];
save([savepath,img_name,'.mat'],'IMG_SPIRiT');

