% 20151127 modified for in-vivo data processing
clear all;close all;clc;
format compact;warning off;
addpath('E:\Simin\GuoLab\GuoLab他人工作\吴鹏\Wu Peng paper\code\Multichannel WFS_20151125\simulation\support_function')
% path_pub = genpath('E:\Simin\GuoLab\Software&Toolbox\吴鹏public-routines');
% addpath(path_pub);
% path_wfs = genpath('E:\work\Fatty_liver\codes\fwtoolbox_v1_code\fwtoolbox_v1_code');
% addpath(path_wfs);
% % Data acquired from phase-array receiver coils
% file_data(1) = {'G:\data\LIVER_FAT\2yzh\CPX\cpx_006'}; % 1 degree
% file_data(2) = {'G:\data\LIVER_FAT\2yzh\CPX\cpx_007'}; % 2 degrees
% file_data(3) = {'G:\data\LIVER_FAT\2yzh\CPX\cpx_005'}; % 5 degrees
% file_data(4) = {'G:\data\LIVER_FAT\2yzh\CPX\cpx_003'}; % 10 degrees
% file_data(5) = {'G:\data\LIVER_FAT\2yzh\CPX\cpx_008'}; % 20 degrees
%
% % Data acquired from body coil
% file_body(1) = {'G:\data\LIVER_FAT\2yzh\CPX\cpx_004'}; % 10 degrees
% % file_body(2) = {'G:\data\channel combination\20151126_LIVER_WUPENG\CPX\cpx_006'};

% Data acquired from phase-array receiver coils
file_data(1) = {'E:\Simin\GuoLab\GuoLab他人工作\吴鹏\Wu Peng paper\data\LIVER_FAT\1xg\raw\raw_005'}; % 1 degree
file_data(2) = {'E:\Simin\GuoLab\GuoLab他人工作\吴鹏\Wu Peng paper\data\LIVER_FAT\1xg\raw\raw_007'}; % 2 degrees
file_data(3) = {'E:\Simin\GuoLab\GuoLab他人工作\吴鹏\Wu Peng paper\data\LIVER_FAT\1xg\raw\raw_006'}; % 5 degrees
file_data(4) = {'E:\Simin\GuoLab\GuoLab他人工作\吴鹏\Wu Peng paper\data\LIVER_FAT\1xg\raw\raw_003'}; % 10 degrees
file_data(5) = {'E:\Simin\GuoLab\GuoLab他人工作\吴鹏\Wu Peng paper\data\LIVER_FAT\1xg\raw\raw_008'}; % 20 degrees

% Data acquired from body coil
file_body(1) = {'E:\Simin\GuoLab\GuoLab他人工作\吴鹏\Wu Peng paper\data\LIVER_FAT\1xg\raw\raw_004'}; % 10 degrees
% file_body(2) = {'G:\data\channel combination\20151126_LIVER_WUPENG\CPX\cpx_006'};

% Parameters
TE = 1.47:0.5:3.97;
chemifre = -434.3;
imDataParams.TE = TE*1e-3;
imDataParams.PrecessionIsClockwise = 1;
imDataParams.FieldStrength = 3;

for file_idx = 4
    % load imaging data
    disp([num2str(file_idx),'.  ',cell2mat(file_data(file_idx))]);

    tic;
    if ~exist ([cell2mat(file_data(file_idx)),'_raw_images.mat'],'file')
        
        [ksig, dim, coil, kRange, kOSFac, kSENSE, slicenum, echonum] = RAWOpener_dyn(cell2mat(file_data(file_idx)));
        ksig = permute(ksig,[1 5 3 2 4 6]);
        recda0 = zeros(size(ksig));
        [x_len, y_len, sli_no, coil_no, te_no, dyn_no]=size(ksig);
        for sli_idx = 1:sli_no
            for coil_idx = 1:coil_no
                for te_idx = 1:te_no
                    for dyn_idx = 1:dyn_no
                        %                      hehe = ifftshift(fft2(fftshift(ksig(:,:,1,5,1,1)),400,222),1);
                        recda0(:,:,sli_idx,coil_idx,te_idx,dyn_idx) = ifftshift(fft2(ifftshift(ksig(:,:,sli_idx,coil_idx,te_idx,dyn_idx))),1);
                    end
                end
            end
        end
        save ([cell2mat(file_data(file_idx)),'_raw_images.mat'],'recda0')
    else
        load ([cell2mat(file_data(file_idx)),'_raw_images.mat'])
    end
    recda0 = recda0(57:344,:,:,:,:,:);
    toc
    disp('Finished loading imaging data!')

    
    % load qbody data
    disp([num2str(file_idx),'.  ',cell2mat(file_body(1))]);

    tic;
    if ~exist ([cell2mat(file_body(1)),'_raw_images.mat'],'file')
        
        [ksig, dim, coil, kRange, kOSFac, kSENSE, slicenum, echonum] = RAWOpener_dyn(cell2mat(file_body(1)));
        ksig = permute(ksig,[1 5 3 2 4 6]);
        qbodyrecda = zeros(size(ksig));
        [x_len, y_len, sli_no, coil_no, te_no, dyn_no]=size(ksig);
        for sli_idx = 1:sli_no
            for coil_idx = 1:coil_no
                for te_idx = 1:te_no
                    for dyn_idx = 1:dyn_no
                        %                      hehe = ifftshift(fft2(fftshift(ksig(:,:,1,5,1,1)),400,222),1);
                        qbodyrecda(:,:,sli_idx,coil_idx,te_idx,dyn_idx) = ifftshift(fft2(ifftshift(ksig(:,:,sli_idx,coil_idx,te_idx,dyn_idx))),1);
                    end
                end
            end
        end
        save ([cell2mat(file_body(1)),'_raw_images.mat'],'qbodyrecda')
    else
        load ([cell2mat(file_body(1)),'_raw_images.mat'])
    end
    qbodyrecda = qbodyrecda(57:344,:,:,:,:,:);
    toc
    disp('Finished loading data!')
    
%     qbodyrecda = read_cpx_ordered(cell2mat(file_body(1)));    
    
    

%     qbodyrecda = qbodyrecda(11:end-10,:,:,:);
%     recda0 = recda0(11:end-10,:,:,:,:,:);% Nx,Ny,Nz,Nte,Ncoil,Ndyn->Nx,Ny,Nz,Ncoil,Nte,Ndyn
    [x_len, y_len, sli_no, coil_no, te_no, dyn_no]=size(recda0);
    
    qbodyrecda = squeeze(mean(qbodyrecda,4));%Nx,Ny,Nz,Ncoil,Nte,Ndyn->Nx,Ny,Nz,Nte,Ndyn
    
    dyn_idx = 1;
    recda = recda0(:,:,:,:,:,dyn_idx);
    qbodyrecda = qbodyrecda(:,:,:,:,dyn_idx);
    sense_map = recda./(repmat(permute(qbodyrecda,[1 2 3 5 4]),[1 1 1 coil_no 1])+1e-6);
    sense_map = mean(sense_map,5);
    sens_m = FitSensMap(sense_map,recda);
    
    save raw_data_sense_mask.mat sens_m
    
end