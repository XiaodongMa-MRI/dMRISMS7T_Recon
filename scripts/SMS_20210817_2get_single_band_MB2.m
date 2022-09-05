%% Data: 7TSMS_Subj2 MB2 
% clear; clc;
%addpath(genpath('C:/Users/GuoLab/Desktop/MR code/'))
set(0,'defaultfigurecolor','w')

%% Para 1
% ACS = 'Gre'; % GRE
% PEShift = 2;
% AP_direction = 'AP';
% if AP_direction(1) == 'A'
%     rotation_90 = 90;
% else
%     rotation_90 = -90;
% end
% filepath = '/home/naxos2-raid1/maxiao/Projects/dMRISMSRecon/Recon/data/Subj5_2021-618/';
% savepath0 = [filepath,'/PEShift',num2str(PEShift),'_',AP_direction,'/'];
% savepath = [savepath0,'single_band/'];
% mkdir(savepath);

%% Para 2
% MB = 2;
% R_PE = 3; % acceleration factor
% PF = 0.75;
% par.kx_r = 200; % r:real
% par.ky_r = 200; 
% par.ky_s = 50; % s:sampled
% par.nCH = 32;
% par.nSL = 140;
% par.nDYN = 1;
% par.nSHOT = 1;

%% Load single band data 
% load ([filepath, 'raw_refs_1mm_',ACS,'ACS_PeShift',num2str(PEShift),'_',AP_direction,'.mat']);
load ([filepath, filename_refs]);

% single band data
kspace0_temp = squeeze(refs.sbfid);
kspace0 = permute(kspace0_temp,[1 2 4 3]); % [x0 y0 ch nsl]
kspace0_bp = kspace0;
S0 = size(kspace0);

% 3 line ref 
ref_temp = squeeze(refs.sbpcor_fid);
ref_ghost_correction =  permute(ref_temp,[1 2 4 3]); % [x0 y0 ch nsl]
disp('Single band data Loaded.')

%% DO phase correction: 3 Line Ref 
if exist([savepath, 'ksp_',ACS,'_sb_3LineRef.mat'],'file')==0
    edata_ref = mean(ref_ghost_correction(:,2:2:end,:,:),2); 
    odata_ref = mean(ref_ghost_correction(:,1:2:end,:,:),2); 
    eImg_ref =  ifftc(edata_ref,1); %x-ky
    oImg_ref =  ifftc(odata_ref,1); %x-ky
    theta_x = angle(eImg_ref./(oImg_ref+eps)); 
    phase_ref = repmat(theta_x/2,[1 par.ky_s 1 1]);

    idx1 = phase_ref>0; idx1(76:end,:,:,:) = 0; % phase unwrap
    idx2 = phase_ref<0; idx2(1:125,:,:,:) = 0;
    phase_ref_unwrap = phase_ref;
    phase_ref_unwrap(idx1) = phase_ref(idx1)-pi;
    phase_ref_unwrap(idx2) = phase_ref(idx2)+pi;

    phase_ref_coil_mean = repmat(mean(phase_ref_unwrap,3),[1 1 S0(3) 1]);
    idx = abs(phase_ref_coil_mean-phase_ref_unwrap)>0.3;
    phase_ref = phase_ref_unwrap;
    phase_ref(idx) = phase_ref_coil_mean(idx); % spike removal

%     figure; 
%     for i = 1:S0(3) 
%         plot(phase_ref(:,1,i,80)); hold on
%     %     plot(phase_Ref_Ent(:,1,i)); hold on
%     end
%     title('Different Coil')
% 
%     figure;
%     for i = 1:5:S0(4) 
%         plot(phase_ref_coil_mean(:,1,1,i)); hold on
%     end
%     title('Different Slice')
    save([savepath, 'ksp_',ACS,'_sb_3LineRef'],'phase_ref','-v7.3');
end

%% DO phase correction: referenceless SVD
if exist([savepath, 'ksp_',ACS,'_sb_SVD.mat'],'file')==0
    tic
    phase_SVD = epighost_fast_low_rank_dep_v2(kspace0(:,:,:,65),par.nSHOT);
    toc
    phase_SVD = repmat(phase_SVD(:,:,1),[1 1 S0(3) S0(4)]);
    save([savepath, 'ksp_',ACS,'_sb_SVD'],'phase_SVD','-v7.3');
end

%% Load Phase --------------------------
fprintf('Single band data ghost correction...it may take long time...\n');  
load([savepath, 'ksp_',ACS,'_sb_SVD']);
load([savepath, 'ksp_',ACS,'_sb_3LineRef']);
phase_ref_cor = medfilt1(phase_ref,5);

% figure; 
% for i = 1:32
%     plot(phase_ref_cor(:,1,i,30));hold on
% end
% plot(mean(phase_ref_cor(:,1,:,30),3));hold on
% plot(phase_SVD(:,1,1,30))

phase_ref_cor_coil_mean = repmat(mean(phase_ref_cor,3),[1 1 32 1]);
phase_ref_abs = phase_ref_cor_coil_mean - phase_SVD;
idx = abs(phase_ref_abs)>0.3;
phase_ref_cor = phase_SVD.*idx + phase_ref_cor.*(1-idx);
phase_ref_cor = medfilt1(phase_ref_cor,5);
% phase_ref_cor = 0.5*phase_ref_cor + 0.5*phase_SVD;

%% Slice Loop --------------------------
for which_slice = 1:size(kspace0_bp,4)
    disp(['Slice ',num2str(which_slice)])
    % 1: Apply 3line-Ref Phase
    kspace0_cor0 = apply_phase(kspace0_bp(:,:,:,which_slice),phase_ref_cor(:,:,:,which_slice),par.nSHOT);    

    % 2: Apply SVD Phase
    kspace0_cor1 = apply_phase(kspace0_bp(:,:,:,which_slice),phase_SVD(:,:,:,which_slice),par.nSHOT);    
    
    % 3: Apply Ref-Ent Phase
    phase_Ref_Ent = zeros([S0(1),S0(2),S0(3)]);
    xunits = (1:S0(1)) - (S0(1)/2) - 0.5;
    for coil_i = 1:S0(3)
        [~,kappaEstSimplex0, phiEstSimplex0] = ghost_correction(kspace0_cor0(:,:,coil_i),'ent',0.05,0,0,0);
%         disp(['Coil ',num2str(coil_i),', kappaEstSimplex=',num2str(kappaEstSimplex0),', phiEstSimplex=',num2str(phiEstSimplex0)])
        phase_tmp(:,1) = (2*pi*kappaEstSimplex0 * xunits./(S0(1))+phiEstSimplex0); 
        phase_tmp = repmat(phase_tmp(:,1), [1,S0(2)]);
        if abs(kappaEstSimplex0)>0.2
            phase_Ref_Ent(:,:,coil_i) = 0;
        elseif abs(kappaEstSimplex0)>0.15
            phase_Ref_Ent(:,:,coil_i) = phase_tmp/3;
        elseif abs(kappaEstSimplex0)>0.1
            phase_Ref_Ent(:,:,coil_i) = phase_tmp/2;
        else
            phase_Ref_Ent(:,:,coil_i) = phase_tmp;        
        end
    end
    phase_Ref_Ent = phase_Ref_Ent + phase_ref_cor(:,:,:,which_slice);
    kspace0_cor2 = apply_phase(kspace0_bp(:,:,:,which_slice),phase_Ref_Ent,par.nSHOT);    
    
%     % Show Result
%     Img(:,:,1) = squeeze(rsos(ifft2c(kspace0_cor0),3));
%     Img(:,:,2) = squeeze(rsos(ifft2c(kspace0_cor1),3));
%     Img(:,:,3) = squeeze(rsos(ifft2c(kspace0_cor2),3));
%     mosaic(imrotate(Img, rotation_90), 1, 3, 2, ['Slice ',num2str(which_slice)], 0.6*genCaxis(Img))
% 
%     mosaic(imrotate(ifft2c(kspace0_cor0), rotation_90), 8, 4, 21, ['Nonlinear, Slice ',num2str(which_slice)], 0.8*genCaxis(ifft2c(kspace0_cor0)))
%     mosaic(imrotate(ifft2c(kspace0_cor1), rotation_90), 8, 4, 22, ['Nonlinear, Slice ',num2str(which_slice)], 0.8*genCaxis(ifft2c(kspace0_cor0)))
%     mosaic(imrotate(ifft2c(kspace0_cor2), rotation_90), 8, 4, 23, ['Nonlinear, Slice ',num2str(which_slice)], 0.8*genCaxis(ifft2c(kspace0_cor0)))
    
%     save results
%     kspace0 = kspace0_cor0;
%     save([savepath, 'ksp_',ACS,'_sb_Ref_Slice',num2str(which_slice),],'kspace0','-v7.3');
    
%     kspace0 = kspace0_cor1;
%     save([savepath, 'ksp_',ACS,'_sb_SVD_Slice',num2str(which_slice),],'kspace0','-v7.3');

    kspace0 = kspace0_cor2;
    save([savepath, 'ksp_',ACS,'_sb_Ref-Ent_Slice',num2str(which_slice),],'kspace0','-v7.3');
    
end

%% Show GRAPPA Results
% load([savepath0,'ksp_ref_Fleet']);
% 
% % Slice Loop -------------------
% for which_slice = 1:140
%     which_slice_ = which_slice; flag_ = 1;
%     if which_slice>70
%         which_slice_ = which_slice - 70;
%         flag_ = 2;
%     end
%     
%     ref0_tmp=squeeze(permute(ksp_ref(:,:,:,which_slice_,:),[1 2 5 3 4]));  %ref0_tmp [200 60 2 32]
%     ref0_tmp_sb = ifftc(ref0_tmp,3); 
%     [P1,P2] = GRAPPA_calibration(ref0_tmp_sb(:,:,flag_,:), 1, 3, 3, 4, -0.03); % P1 is related to MB
% 
%     
%     disp(['Slice ',num2str(which_slice),' Loading....' ])
%     img2 = zeros(200,150,4);
%     for method_i = 3
%         switch method_i
%             case 1
%                 Img1 = load([savepath,'ksp_',ACS,'_sb_Ref_Slice',num2str(which_slice)]);
%             case 2
%                 Img1 = load([savepath,'ksp_',ACS,'_sb_SVD_Slice',num2str(which_slice)]);
%             case 3
%                 Img1 = load([savepath,'ksp_',ACS,'_sb_Ref-Ent_Slice',num2str(which_slice)]);
%         end
%         if which_slice>70
%             Img1.kspace0(:,2:2:end) = Img1.kspace0(:,2:2:end).*exp(1i*pi);
%         end
%         K0 = zeros(200,150,1,32);
%         K0(:,1:3:end,:,:) = Img1.kspace0;
%         k2 =  GRAPPA_recovery(P1,P2,K0,1,3,3,4,0);
%         img2(:,:,method_i) = squeeze(rsos(ifft2c(k2),4));    
%   
%     end
% %     mosaic(imrotate(img2(:,:,method_i), rotation_90), 1, 1, method_i, ['Slice',num2str(which_slice)], 1*genCaxis(img2))
% end
