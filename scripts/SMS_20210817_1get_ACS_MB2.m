% clear; clc;
% %addpath(genpath('C:/Users/GuoLab/Desktop/MR code/'))

%% Para 1
% PEShift = 2;
% AP_direction = 'AP';
% filepath = '/home/naxos2-raid1/maxiao/Projects/dMRISMSRecon/Recon/data/Subj5_2021-618/';
% savepath = [filepath,'/PEShift',num2str(PEShift),'_',AP_direction,'/'];
% mkdir( savepath )
% 
%% Para 2
% MB = 2;
% R_PE = 3; % acceleration factor
% PF = 0.75;
% par.kx_r = 200; % r:real
% par.ky_r = 200; 
% par.ky_s = 50; % s:sampled
% par.nCH = 32;
% par.nSL = 70;

%% 1. Fleet ref data for GRAPPA recon
if strcmp(ACS,'Fleet')
    % ksp_ref_FLEET.mat: ksp_ref ( 200kx  66ky   32coil   66slice   2kz )
    % >>>>>>>>>>>>>>>>>>>>> 1. ref data for GRAPPA recon
    %load ([filepath, 'raw_refs_1mm_FleetACS_PeShift',num2str(PEShift),'_',AP_direction,'.mat']);
%     load ([filepath, 'raw_refs_1mm_FleetACS_PeShift',num2str(PEShift),'_',AP_direction,'.mat']);
    load ([filepath, filename_refs]);
    % load ([filepath, 'raw_refs_1mm_FleetACS_',AP_direction,'.mat']);

    ksp_ref_temp = zeros(size(refs.sbpatref.fid));
    ksp_ref_temp(:,1:3:end,:,:,:,:) = refs.sbpatref.fid(:,1:22,:,:,:,:);
    ksp_ref_temp(:,2:3:end,:,:,:,:) = refs.sbpatref.fid(:,23:44,:,:,:,:);
    ksp_ref_temp(:,3:3:end,:,:,:,:) = refs.sbpatref.fid(:,45:66,:,:,:,:);
    ksp_ref = permute(squeeze(ksp_ref_temp),[1 2 4 3]); % ref in [x y ch nsl]
    % myshow3(squeeze(sos(ifft2c(ksp_ref),3))); %imcontrast

    % >>>>>>>>>>>>>>>>>>>> 2. Phase Correction for ref data
    phase0 = epighost_fast_low_rank_dep_v2(ksp_ref(:, :, :,70),3);%center_slice_ind,nshot
    % save([savepath,'phase_Fleet_ref_slice70'],'phase0','-v7.3')
    % load ([savepath, 'phase_Fleet_ref_slice70']);
    ksp_ref_temp = apply_phase(ksp_ref,phase0, 3);% Apply phase
    % myshow3(squeeze(sos(ifft2c(ksp_ref_temp),3))); %imcontrast
    % save([savepath, 'ksp_ref_FLEET_temp'],'ksp_ref_temp');

    % >>>>>>>>>>>>>>>>>>> 3. z-->kz
    % load([savepath, 'ksp_ref_FLEET_temp']);
    ksp_ref_temp2 = zeros(size(ksp_ref_temp,1),size(ksp_ref_temp,2),size(ksp_ref_temp,3),size(ksp_ref_temp,4)/MB,MB);
    for i = 1:size(ksp_ref_temp,4)/MB
    %     ksp_ref_temp2(:,:,:,i,:) = ksp_ref_temp(:,:,:,[i,i+size(ksp_ref_temp,4)/MB,i+size(ksp_ref_temp,4)/MB*2]);     
        ksp_ref_temp2(:,:,:,i,:) = ksp_ref_temp(:,:,:,[i,i+size(ksp_ref_temp,4)/MB]);
    end
    ksp_ref = fftc(ksp_ref_temp2,5);
%     myshow3(squeeze(sos(ifftc(ifft2c(ksp_ref),5),3))); %imcontrast
    clear refs; clear ksp_ref_temp2; clear ksp_ref_temp
    save([savepath, 'ksp_ref_Fleet'],'ksp_ref');

%% 1. Gre ref data for GRAPPA recon
% ksp_ref_GRE.mat: ksp_ref ( 200kx  66ky   32coil   70slice   2kz )
% >>>>>>>>>>>>>>>>>>>>> 1. ref data for GRAPPA recon
else
%     load ([filepath, 'raw_refs_1mm_GreACS_PeShift',num2str(PEShift),'_',AP_direction,'.mat']);
    load ([filepath, filename_refs]);
    % load ([filepath, 'raw_refs_1mm_GreACS_PeShift',num2str(PEShift),'_',AP_direction,'.mat']);
    ksp_ref_temp = permute(squeeze(refs.sbpatref.fid),[1 2 4 3]);
    % ksp_ref_temp = [x200 y66 ch32 nsl140]
    % myshow3(sos(ifft2c(ksp_ref_temp),3))
    ksp_ref_temp2 = zeros(size(ksp_ref_temp,1),size(ksp_ref_temp,2),size(ksp_ref_temp,3),size(ksp_ref_temp,4)/MB,MB);
    for i = 1:size(ksp_ref_temp,4)/MB
    %     ksp_ref_temp2(:,:,:,i,:) = ksp_ref_temp(:,:,:,[i,i+size(ksp_ref_temp,4)/MB,i+size(ksp_ref_temp,4)/MB*2]);
        ksp_ref_temp2(:,:,:,i,:) = ksp_ref_temp(:,:,:,[i,i+size(ksp_ref_temp,4)/MB]);
    end
    ksp_ref = fftc(ksp_ref_temp2,5);
    clear refs; clear ksp_ref_temp2; clear ksp_ref_temp
%     myshow3(squeeze(sos(ifftc(ifft2c(ksp_ref),5),3)));
    save([savepath, 'ksp_ref_Gre'],'ksp_ref');
    %ksp_ref = 200 66 32 70 2
end
