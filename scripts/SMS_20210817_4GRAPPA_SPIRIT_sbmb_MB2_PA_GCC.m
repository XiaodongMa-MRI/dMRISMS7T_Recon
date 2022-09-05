clear; clc;
%addpath(genpath('C:/Users/GuoLab/Desktop/MR code/'))
set(0,'defaultfigurecolor','w')

%% Para 1
PEShift = 2; % only 2
which_ghostcorr_method = {'Ref','SVD','Ref-Ent'}; 
AP_direction = 'PA';
if AP_direction(1) == 'A'
    rotation_90 = 90;
else
    rotation_90 = -90;
end

filepath = ['/home/naxos2-raid1/maxiao/Projects/dMRISMSRecon/Recon/data/Subj5_2021-618/','PEShift',num2str(PEShift),'_',AP_direction,'/single_band/'];
savepath = ['/home/naxos2-raid1/maxiao/Projects/dMRISMSRecon/Recon/data/Subj5_2021-618/','PEShift',num2str(PEShift),'_',AP_direction,'/'];
savepath2 = ['/home/naxos2-raid1/maxiao/Projects/dMRISMSRecon/Recon/data/Subj5_2021-618/PEShift',num2str(PEShift),'_',AP_direction,'/img_SMS_PF_SoS/'];
mkdir(savepath2)

%% Para 2
MB = 2;
R_PE = 3; % acceleration factor
PF = 0.75;
par.kx_r = 200;
par.nCH = 32;
slice_gap = 70;
par.ky_r = 200; 
par.ky_s = 150;
par.ky_sr = 50;
par.nSL = 2;

%% Recon Par
GRAPPA = 1; GRAPPA_sb = 1; SPIRIT = 1; 

%% Big Loop
ACS_name = {'Fleet','Gre'};
%ACS_name = {'Gre'};
for i = 2:length(ACS_name)
    ACS = ACS_name{i};
    load([savepath, 'ksp_ref_',ACS]);   
    for plus_number = 1:3
        plus_ = ['_',num2str(plus_number)];
        %% slice loop
%        for which_slice = 50%56:10:66%56:10:66%:10:66
        for which_slice = 1:70%56:10:66%56:10:66%:10:66
            for which_ghost_correction_method = 3  % which ghost correction method
                   %% load single band data [kx200 ky50  ch32  kz2]
                par.nCH = 32;
                if GRAPPA_sb 
                    single_band_ref = zeros([par.kx_r,par.ky_sr,par.nCH,MB]);
                    load([filepath,'ksp_',ACS,'_sb_',which_ghostcorr_method{which_ghost_correction_method},'_Slice',num2str(which_slice)]);
                    single_band_ref(:,:,:,1) = kspace0;          
                    load([filepath,'ksp_',ACS,'_sb_',which_ghostcorr_method{which_ghost_correction_method},'_Slice',num2str(which_slice + slice_gap)]);
                    kspace0(:,2:2:end,:,:,:)= kspace0(:,2:2:end,:,:,:).*exp(1i*(1)*pi);
                    single_band_ref(:,:,:,2) = kspace0;        
                    single_band_ref = permute(single_band_ref,[1 2 4 3]);
                %             myshow3(squeeze(rsos((ifft2c(single_band_ref)),4)));
                end

                    %% Multi-band Data  [kx200 ky150 kz2 ch32 nb3]
                load([savepath, 'ksp_','Gre','_',which_ghostcorr_method{which_ghost_correction_method},'_Slice',num2str(which_slice),plus_]);
                kspace0(:,2:2:end,:,:,:)= kspace0(:,2:2:end,:,:,:).*exp(1i*pi);

                kspace = zeros(par.kx_r,par.ky_sr*R_PE,par.nCH,size(kspace0,4),1,MB);
                for mb=1:MB
                    kspace(:,(mb-1)*R_PE+1:MB*R_PE:end,:,:,:,mb) = kspace0(:,mb:MB:end,:,:,:);
                end
                kspace = permute(kspace,[1 2 6 3 4 5]); 
                % kspace in [kx200 ky150 kz2 ch32 nb30]

              %% Prepare the 3D GRAPPA ACS  [ref0_tmp + single_band_ref2]
                % get single_band_ref (In-Plane GRAPPA)
                if MB>1
                    ref0_tmp = squeeze(permute(ksp_ref(:,:,:,which_slice,:),[1 2 5 3 4]));  %ref0_tmp [200 60 2 32]
                else
                    ref0_tmp = ksp_ref(:,:,which_slice,:);
                end
                %  myshow3(squeeze(sos(ifftc(ifft2c(ref0_tmp),3),4)));

               if GRAPPA_sb
                    ref0_tmp_sb = ifftc(ref0_tmp,3); 
                    S1 = size(single_band_ref);
                    single_band_ref_sb_kspace = zeros(S1(1),S1(2)*R_PE,S1(3),S1(4));
                    single_band_ref_sb_kspace(:,1:R_PE:end,:,:) = single_band_ref;
					clear single_band_ref2
                    [P1,P2] = GRAPPA_calibration(ref0_tmp_sb(:,:,1,:), 1, R_PE, 3, 4, 0); % P1 is related to MB
                    single_band_ref2(:,:,1,:) = GRAPPA_recovery(P1,P2,single_band_ref_sb_kspace(:,:,1,:),1,R_PE,3,4,0); % dkspace1 in [x y z ch]           
                    [P1,P2] = GRAPPA_calibration(ref0_tmp_sb(:,:,2,:), 1, R_PE, 3, 4, 0); % P1 is related to MB
                    single_band_ref2(:,:,2,:) = GRAPPA_recovery(P1,P2,single_band_ref_sb_kspace(:,:,2,:),1,R_PE,3,4,0); % dkspace1 in [x y z ch]
                    single_band_ref2 = fftc(single_band_ref2,3);
    %                 img_temp = (sos(ifftc(ifft2c(single_band_ref2),3),4));  
    %                 mosaic(imrotate(img_temp, rotation_90), 1, 2, which_ghost_correction_method, ['single band ',num2str(which_slice)], 0.6*genCaxis(img_temp))
               end
               
               %% Coil combination
               par.nCH = 16;

               if GRAPPA_sb
                    [ref0_GCC,single_band_ref_GCC,kspace0_GCC] = coil_combination_for_MB_data_sb(ref0_tmp,single_band_ref2,kspace,par.nCH);                    
                    ref0_tmp = ref0_GCC;
                    kspace = kspace0_GCC;
                    single_band_ref2 = single_band_ref_GCC;
               else
                    [ref0_GCC,kspace0_GCC] = coil_combination_for_MB_data(ref0_tmp,kspace,par.nCH );
                    ref0_tmp = ref0_GCC;
                    kspace = kspace0_GCC;
                end
                
              %% GRAPPA
                if GRAPPA
                    ker_x=3; ker_y=4; lamb0=0.1; % ker_y is chosen according to the R_PE 
                    fprintf('Reconstruction slice %d ...\n',which_slice)
                    tic;
                    [P11,P21] = GRAPPA_calibration(ref0_tmp, MB, R_PE, ker_x, ker_y, lamb0); % use GRE as ref only 
                    toc

                    fprintf('Recovering data...\n');
                    tic
                    dkspace1_all_tmp1 = kspace;
                    for nb = 1:size(kspace0,4)
                        fprintf(['nb = ',num2str(nb),'...\n']);
                        dkspace1_tmp1 = GRAPPA_recovery(P11,P21, kspace(:,:,:,:,nb),MB,R_PE,ker_x,ker_y,lamb0); % dkspace1 in [x y z ch]
                        dkspace1_all_tmp1(:,:,:,:,nb) = dkspace1_tmp1;
                    end                   
                    toc
                    
                    tic
                    complex_data1 =  PF_coil_combination_GRAPPA(dkspace1_all_tmp1,par);   
                    toc
                    
                    img = single(squeeze(complex_data1(:,:,1,1,:)));
                    img_name1 = [ACS,'_img_GRAPPA_PF_ACC_','slice_',num2str(which_slice),'_GhostMethod',num2str(which_ghost_correction_method),plus_];
                    save([savepath2,img_name1],'img','-v7.3')
                    img = single(squeeze(complex_data1(:,:,2,1,:)));
                    img_name2 = [ACS,'_img_GRAPPA_PF_ACC_','slice_',num2str(which_slice + slice_gap),'_GhostMethod',num2str(which_ghost_correction_method),plus_];
                    save([savepath2,img_name2,],'img','-v7.3')
                    
                    % Show Results
                    img_GRAPPA1 = reshape(squeeze(complex_data1(:,:,:,:,1:2)),[size(complex_data1,1),size(complex_data1,2),2*MB]);
                    %mosaic(imrotate(img_GRAPPA1, rotation_90), 1, MB, 70+which_ghost_correction_method, [''], 1.2*genCaxis(img_GRAPPA1)) % ref ACS                    
                    
                end
              %% GRAPPA_sb
                if GRAPPA_sb
                    ker_x = 3; ker_y = 4; lamb0 = 0.1; % ker_y is chosen according to the R_PE 
                    fprintf('Reconstruction slice %d ...\n',which_slice)
                    tic;
                    [P12,P22] = GRAPPA_calibration(single_band_ref2, MB, R_PE, ker_x, ker_y, lamb0); %  use recon single-band as ref
                    toc

                    fprintf('Recovering data...\n');
                    tic
                    dkspace1_all_tmp2 = kspace;
                    for nb = 1:size(kspace0,4)
                        fprintf(['nb = ',num2str(nb),'...\n']);
                        dkspace1_tmp2 = GRAPPA_recovery(P12,P22,kspace(:,:,:,:,nb),MB,R_PE,ker_x,ker_y,lamb0); % dkspace1 in [x y z ch]
                        dkspace1_all_tmp2(:,:,:,:,nb) = dkspace1_tmp2;
                    end                   
                    toc
                    
                    tic
                    complex_data2 =  PF_coil_combination_GRAPPA(dkspace1_all_tmp2,par);   
                    toc
                    
                    img = single(squeeze(complex_data2(:,:,1,1,:)));
                    img_name1 = [ACS,'_img_GRAPPA_PF_ACC_','slicesb_',num2str(which_slice),'_GhostMethod',num2str(which_ghost_correction_method),plus_];
                    save([savepath2,img_name1],'img','-v7.3')
                    img = single(squeeze(complex_data2(:,:,2,1,:)));
                    img_name2 = [ACS,'_img_GRAPPA_PF_ACC_','slicesb_',num2str(which_slice + slice_gap),'_GhostMethod',num2str(which_ghost_correction_method),plus_];
                    save([savepath2,img_name2,],'img','-v7.3')
                    
                    % Show Results
                    img_GRAPPA2 = reshape(squeeze(complex_data2(:,:,:,:,1:2)),[size(complex_data2,1),size(complex_data2,2),2*MB]);
                    %mosaic(imrotate(img_GRAPPA2, rotation_90),1, MB, 80+which_ghost_correction_method, [''], 1.2*genCaxis(img_GRAPPA2)) % ref ACS                    
                    
                end


                %% L1 Spirit Recon
                if SPIRIT        
                    ACS_DATA = permute(ref0_tmp,[1 2 4 3]);  
                    Ry = 6;                     % accl factor
                    del_step = 3;               % amount of k-space staggering between phase-cycles
                    N = [200 150];
                    num_acsX = size(ACS_DATA,1);            % acs size in readout
                    num_acsY = size(ACS_DATA,2);              % acs size in PE
                    num_cycle = 2;
                    num_chan = par.nCH;
                    CalibSize = [num_acsX, num_acsY];     
                    mask = zeros([N, num_cycle]);
                    del = mod((0:num_cycle-1) * del_step, Ry);       % staggering amount

                    for t = 1:num_cycle
                        mask(:,1+del(t):Ry:end,t) = 1;
                    end

                    Mask_spirit = ( mask ) ~= 0;
                    Mask_spirit = repmat(permute(Mask_spirit, [1,2,4,3]), [1,1,num_chan,1]);
                    % myshow3(rsos(ifft2c(ACS_DATA),3))
                    ACS_DATA = reshape( ACS_DATA, [CalibSize, num_chan * num_cycle] );

                    % Perform Calibration   
                    kSize = [3,3];              % SPIRiT kernel size
                    CalibTyk = 2.5e-11;%8e-15;            % Tykhonov regularization in the calibration
                    kCalib = ACS_DATA;
                    tic
                        kernel = calibSPIRiT(kCalib, kSize, num_chan*num_cycle, CalibTyk);
                    toc
                    % save([savepath, 'kernel_3x5'],'kernel','-v7.3' );
                    % load([savepath, 'kernel_3x5'])
                    GOP = SPIRiT(kernel, 'fft', N);

                    lambda = 5e-10 %5e-13;              % TV penalty
                    Lambda_ = '_Lambda5e-08';
                    L = 1;                      % set "L=1" for TV penalty on each coil separately
                                                % set "L=12" for joint TV penalty across all coils/cycles
                    iter_inner = 10;            % outer cg loops
                    iter_outer = 6;            % inner cg loops
                    for nb = 1 :size(kspace0,4)
                        % input slice: kx200 ky150 ch32 kz2
                        kspace0_slice = permute(squeeze((kspace(:,:,:,:,nb))),[1 2 4 3]);

                        % undersampled data
                        DATA_sampled = (kspace0_slice).*Mask_spirit;
                        DATA_Sampled = reshape( DATA_sampled, [N, num_chan * num_cycle] );
                        % myshow3(rsos(ifft2c(DATA_sampled),3))

                        idx_mis = find(abs(DATA_Sampled)==0);
                        param = init;
                        param.dsp = 1;
                        param.GOP = GOP;            % (G-I) spirit operator
                        param.TV = TViFFT;          % (Gradient * IFFT) operator
                        param.data = -(param.GOP * DATA_Sampled);   
                        param.grad = param.TV * DATA_Sampled;
                        param.N = N;
                        param.num_chan = num_chan * num_cycle;
                        param.TVWeight = lambda;     
                        param.Itnlim = iter_inner;
                        param.idx_mis = idx_mis;
                        param.L = L;
                        res = zeros(size(idx_mis));     % unknown: missing k-space samples
                        Res = DATA_Sampled;             % acquired k-space samples

                        tic
                        for n = 1:iter_outer
                            disp(['iter outer ',num2str(n)])

                            res = fnlCg_l1spirit(res, param);

                            Res(idx_mis) = res;
                            img_jspirit = ifftc( reshape(ifft2c(Res), [N, num_chan, num_cycle]),4);

                            %mosaic(imrotate(rsos(img_jspirit, 3), rotation_90), 1, num_cycle, 99, ['Joint L1-SPIRiT R=', num2str(Ry)], 0.8*genCaxis(rsos(img_jspirit, 3)))

                        end
                        toc
                        dkspace1_tmp = reshape(Res,[N, num_chan, num_cycle]);

                         % PF+ACC+Save
                        tic
                            complex_data =  PF_coil_combination_L1SPIRiT(dkspace1_tmp,par,0);   
                        toc
                        
                        img = single(squeeze(complex_data(:,:,1,1,:)));
                        img_name1 = [ACS,'_img_L1-SPIRiT_PF_ACC_',num2str(which_slice),'_nb',num2str(nb),'_GhostMethod',num2str(which_ghost_correction_method),plus_,Lambda_];
                        save([savepath2,img_name1],'img','-v7.3')

                        img = single(squeeze(complex_data(:,:,2,1,:)));
                        img_name2 = [ACS,'_img_L1-SPIRiT_PF_ACC_',num2str(which_slice + slice_gap),'_nb',num2str(nb),'_GhostMethod',num2str(which_ghost_correction_method),plus_,Lambda_];
                        save([savepath2,img_name2,],'img','-v7.3')


                        %mosaic(imrotate(complex_data(:,:,:,1), rotation_90), 1, par.nSL, 90+which_ghost_correction_method, ['Slice',num2str(which_slice),', nb',num2str(nb)], 1.0*genCaxis(complex_data(:,:,:,1)))
                       
                    end

                end

            end



        end
    end
end