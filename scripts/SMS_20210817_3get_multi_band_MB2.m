% clear; clc;
%addpath(genpath('C:/Users/GuoLab/Desktop/MR code/'))
set(0,'defaultfigurecolor','w')

%% Para 1
% ACS = 'Gre'; % 
% PEShift = 2;
% AP_direction = 'AP';
% if AP_direction(1) == 'A'
%     rotation_90 = 90;
% else
%     rotation_90 = -90;
% end
% 
% filepath = '/home/naxos2-raid1/maxiao/Projects/dMRISMSRecon/Recon/data/Subj5_2021-618/';
% savepath = [filepath,'/PEShift',num2str(PEShift),'_',AP_direction,'/'];

%% Para 2
% MB = 2;
% R_PE = 3; % acceleration factor
% PF = 0.75;
% DWI = true;  % for DWI 
% par.kx_r = 200; % r:real
% par.ky_r = 200; 
% par.ky_s = 50; % s:sampled
% par.nCH = 32;
% par.nSL = 70;
% par.nDYN = 1;
% par.nSHOT = 1;

%% DO phase correction: referenceless SVD
if exist([savepath, 'ksp_Gre_SVD.mat'],'file')==0
    load ([filepath, filename]);
    kspace0 = permute(squeeze(fidsort),[1 2 4 5 3]); % in [x y ch nb nsl]
    tic
    phase_SVD = epighost_fast_low_rank_dep_v2(kspace0(:,:,:,1,35),par.nSHOT);
    toc
    phase_SVD = repmat(phase_SVD(:,:,1),[1 1 par.nCH  par.nSL]);
    save([savepath, 'ksp_Gre_SVD'],'phase_SVD','-v7.3');
end

%% DO phase correction: 3 Line Ref 
if exist([savepath, 'ksp_',ACS,'_3LineRef.mat'],'file')==0
    load ([filepath, filename]);
    % 3 line ref 
    ref_temp = squeeze(refs.pcor_fid);
    ref_ghost_correction =  permute(ref_temp,[1 2 4 3 5]); % [x0 y0 ch nsl]
    ref_ghost_correction = squeeze(mean(ref_ghost_correction,5));
    edata_ref = mean(ref_ghost_correction(:,2:2:end,:,:),2); 
    odata_ref = mean(ref_ghost_correction(:,1:2:end,:,:),2); 
    eImg_ref =  ifftc(edata_ref,1); %x-ky
    oImg_ref =  ifftc(odata_ref,1); %x-ky
    theta_x = angle(eImg_ref./(oImg_ref+eps)); 
    phase_ref = repmat(theta_x/2,[1 par.ky_s 1 1]);
    S0 = size(phase_ref);

    idx1 = phase_ref>0; idx1(76:end,:,:,:) = 0; % phase unwrap
    idx2 = phase_ref<0; idx2(1:125,:,:,:) = 0;
    phase_ref_unwrap = phase_ref;
    phase_ref_unwrap(idx1) = phase_ref(idx1)-pi;
    phase_ref_unwrap(idx2) = phase_ref(idx2)+pi;

    phase_ref_coil_mean = repmat(mean(phase_ref_unwrap,3),[1 1 S0(3) 1]);
    idx = abs(phase_ref_coil_mean-phase_ref_unwrap)>0.3;
    phase_ref = phase_ref_unwrap;
    phase_ref(idx) = phase_ref_coil_mean(idx);% spike removal

%     figure; 
%     for i = 1:S0(3) 
%         plot(phase_ref(:,1,i,35)); hold on
%     end
%     title('Different Coil')
% 
%     figure;
%     for i = 1:5:S0(4) 
%         plot(phase_ref_coil_mean(:,1,1,i)); hold on
%     end
%     title('Different Slice')

    save([savepath, 'ksp_',ACS,'_3LineRef'],'phase_ref','-v7.3');
end
%% Ghost Correction 
for plus_number = 1
    plus_ = ['_',num2str(plus_number)];
    
    %% Load Data
    switch plus_number
        case 1
            load ([filepath, filename]);
%         case 2
%             load ([filepath, 'raw_fidsort_1mm_GreACS_PeShift',num2str(PEShift),'_',AP_direction,'_slice31to60.mat']);
%         case 3
%             load ([filepath, 'raw_fidsort_1mm_GreACS_PeShift',num2str(PEShift),'_',AP_direction,'_slice61to71.mat']);
    end

    
    % fid = 200 50 70slice 1 1 32 30nb(every 10 light 1)
    kspace0 = permute(squeeze(fidsort),[1 2 4 5 3]); % in [x y ch nb nsl]
    kspace0_bp = kspace0;
    sz_GRAPPA = size(kspace0);
    clear fidsort
    % myshow3(squeeze(sos(ifft2c(kspace0(:,:,:,1,:)),3))); 
    % kspace0 = double(kspace0);
    
    %% Load data
    fprintf('ghost correction...it may take long time...\n');  
    load([savepath, 'ksp_Gre_SVD']);
    load([savepath, 'ksp_',ACS,'_3LineRef']);
    phase_ref_cor = medfilt1(phase_ref,5);
%     figure; 
%     for i = 1:32
%         plot(phase_ref_cor(:,1,i,30));hold on
%     end
%     plot(mean(phase_ref_cor(:,1,:,30),3));hold on
%     plot(phase_SVD(:,1,1,30))

    phase_ref_cor_coil_mean = repmat(mean(phase_ref_cor,3),[1 1 32 1]);
    phase_ref_abs = phase_ref_cor_coil_mean - phase_SVD;
    idx = abs(phase_ref_abs)>0.3;
    phase_ref_cor = phase_SVD.*idx + phase_ref_cor.*(1-idx);
    phase_ref_cor = medfilt1(phase_ref_cor,5);
%     phase_ref_cor = 0.5*phase_ref_cor + 0.5*phase_SVD;

    
        %% DO phase correction
    for which_slice = 1:par.nSL
        disp(['slice ',num2str(which_slice), ' processing'])
        % 1: Apply Ref Phase
        kspace0_cor0 = squeeze(kspace0_bp(:,:,:,:,which_slice));
        for b_flag = 1:size(kspace0_bp,4)
            kspace0_cor0(:,:,:,b_flag) = apply_phase(kspace0_bp(:,:,:,b_flag,which_slice),phase_ref_cor(:,:,:,which_slice),par.nSHOT);
        end

        % 2: Apply SVD Phase
%         kspace0_cor1 = apply_phase(kspace0_bp(:,:,:,:,which_slice),phase_SVD(:,:,:,which_slice),par.nSHOT);            

        % 3: Apply Ref-Ent Phase
        phase_Ref_Ent = zeros([sz_GRAPPA(1),sz_GRAPPA(2),sz_GRAPPA(3)]);
        xunits = (1:sz_GRAPPA(1)) - (sz_GRAPPA(1)/2) - 0.5;
        for coil_i = 1:sz_GRAPPA(3)
            [~,kappaEstSimplex0, phiEstSimplex0] = ghost_correction(kspace0_cor0(:,:,coil_i,1),'ent',0.05,0,0,0);
    %             disp(['Coil ',num2str(coil_i),', kappaEstSimplex=',num2str(kappaEstSimplex0),', phiEstSimplex=',num2str(phiEstSimplex0)])
            phase_tmp(:,1) = (2*pi*kappaEstSimplex0 * xunits./(sz_GRAPPA(1))+phiEstSimplex0); 
            phase_tmp = repmat(phase_tmp(:,1), [1,sz_GRAPPA(2)]);
            if abs(kappaEstSimplex0)>0.15
                phase_Ref_Ent(:,:,coil_i) = 0;
            elseif abs(kappaEstSimplex0)>0.1
                phase_Ref_Ent(:,:,coil_i) = phase_tmp/2;
            else
                phase_Ref_Ent(:,:,coil_i) = phase_tmp;        
            end

        end
        phase_Ref_Ent = phase_Ref_Ent + phase_ref_cor(:,:,:,which_slice);
        kspace0_cor2 = squeeze(kspace0_bp(:,:,:,:,which_slice));
        for b_flag = 1:size(kspace0_bp,4)
            kspace0_cor2(:,:,:,b_flag) = apply_phase(kspace0_bp(:,:,:,b_flag,which_slice),phase_Ref_Ent,par.nSHOT);
        end


        % Show Result
%         Img(:,:,1) = squeeze(rsos(ifft2c(kspace0_cor0(:,:,:,1)),3));
%         Img(:,:,2) = squeeze(rsos(ifft2c(kspace0_cor1(:,:,:,1)),3));
%         Img(:,:,3) = squeeze(rsos(ifft2c(kspace0_cor2(:,:,:,1)),3));
        %     mosaic(imrotate(ifft2c(kspace0_cor0(:,:,:,1)), rotation_90), 8, 4, 91, ['Slice ',num2str(which_slice)], 0.6*genCaxis(ifft2c(kspace0_cor0(:,:,:,1))))
        %     mosaic(imrotate(ifft2c(kspace0_cor1(:,:,:,1)), rotation_90), 8, 4, 92, ['Slice ',num2str(which_slice)], 0.6*genCaxis(ifft2c(kspace0_cor1(:,:,:,1))))
        %     mosaic(imrotate(ifft2c(kspace0_cor2(:,:,:,1)), rotation_90), 8, 4, 93, ['Slice ',num2str(which_slice)], 0.6*genCaxis(ifft2c(kspace0_cor2(:,:,:,1))))

        % save results
        %     kspace0 = kspace0_cor0;
        %     save([savepath, 'ksp_',ACS,'_Ref_Slice',num2str(which_slice),plus_],'kspace0','-v7.3');
        %     kspace0 = kspace0_cor1;
        %     save([savepath, 'ksp_',ACS,'_SVD_Slice',num2str(which_slice),plus_],'kspace0','-v7.3');
        kspace0 = kspace0_cor2;
        save([savepath, 'ksp_',ACS,'_Ref-Ent_Slice',num2str(which_slice),plus_],'kspace0','-v7.3');

    end
end

