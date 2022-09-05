%--------------------------------------------------------------------------
%% load fully sampled data
%--------------------------------------------------------------------------
 
set(0,'DefaultFigureWindowStyle','docked')
% 使得matlab的figures不弹出，而是在主窗口显示


% fully-sampled bSSFP data with four phase cycles
% coil compressed to 12 channels using GCC (T Zhang et al MRM 2013)
load IMG_gcc
% size = 160(readout) x 160(phase encoding) x 12(coils) x 4(phase)
% IMG_gcc = IMG_gcc(:,:,:,1);
[N(1), N(2), num_chan, num_cycle] = size(IMG_gcc);


mosaic(imrotate(rsos(IMG_gcc, 3), -90), 1, num_cycle, 1, 'GCC R=1', genCaxis(rsos(IMG_gcc, 3)))
% 一个好用的画图函数mosaic；sos函数rsos；genCaxis：寻找合适的图像范围
% mosaic( imgs, row_num, col_num, fig_num, title_str, disp_range )

mip_full = max(rsos(IMG_gcc, 3), [], 4);        % MIP from fully-sampled data


%%------------------------------------------------------------------------- 
%% Grappa with 1D acceleration
%--------------------------------------------------------------------------


Ry = 6;                     % acceleration factor
 
num_acsX = N(1);            % acs size in readout
num_acsY = 32;              % acs size in phase encoding
 

lambda_tik = 1e-3;          % Tikhonov reg parameter for GRAPPA kernel estimation


% sampling and acs masks
mask = zeros(N);
mask_acs = zeros(N);

mask(:,1:Ry:end) = 1;
mask_acs(1+end/2-num_acsX/2:end/2+num_acsX/2, 1+end/2-num_acsY/2:end/2+num_acsY/2) = 1;



kernel_size = [13,3];       % assume odd kernel size
kernel_hsize = (kernel_size-1)/2;


pad_size = kernel_hsize .* [1,Ry];
N_pad = N + 2*pad_size;
% k-sapce是有pad（补边缘）的，对应后边的padarray

% k-space limits for training:
ky_begin = 1 + Ry * kernel_hsize(2);       % first kernel center point that fits acs region 

ky_end = num_acsY - Ry * kernel_hsize(2);  % last kernel center point that fits acs region 


kx_begin = 1 + kernel_hsize(1);            % first kernel center point that fits acs region 

kx_end = num_acsX - kernel_hsize(1);       % last kernel center point that fits acs region 

 


% k-space limits for recon:
Ky_begin = 1 + Ry * kernel_hsize(2);       % first kernel center point that fits acs region 

Ky_end = N_pad(2) - Ry * kernel_hsize(2);  % last kernel center point that fits acs region 


Kx_begin = 1 + kernel_hsize(1);            % first kernel center point that fits acs region 

Kx_end = N_pad(1) - kernel_hsize(1);       % last kernel center point that fits acs region 

 


% count the no of kernels that fit in acs 
% 计算acs有多少个kernel
ind = 1;

for ky = ky_begin : ky_end
    for kx = kx_begin : kx_end
        ind = ind + 1;        
    end
end


num_ind = ind;



Img_recon = zeros([N, num_chan, num_cycle]);

 

for pc = 1:num_cycle

    kspace_full = fft2c(IMG_gcc(:,:,:,pc));

    kspace_sampled = kspace_full .* repmat(mask, [1,1,num_chan]);

    kspace_acs = kspace_full .* repmat(mask_acs, [1,1,num_chan]);



    % train kernel
    kspace_acs_crop = kspace_acs(1+end/2-num_acsX/2:end/2+num_acsX/2, 1+end/2-num_acsY/2:end/2+num_acsY/2, :);
    % size(160 x 32 x 12)

    Rhs = zeros([num_ind, num_chan, Ry-1]);
    Acs = zeros([num_ind, prod(kernel_size) * num_chan]);
    
    size(Acs)
    

    ind = 1;

    for ky = ky_begin : ky_end
        for kx = kx_begin : kx_end

            acs = kspace_acs_crop(kx-kernel_hsize(1):kx+kernel_hsize(1), ky-kernel_hsize(2)*Ry:Ry:ky+kernel_hsize(2)*Ry, :);

            Acs(ind,:) = acs(:);
 
            for ry = 1:Ry-1
                Rhs(ind,:,ry) = kspace_acs_crop(kx, ky-ry, :);
            end    

            ind = ind + 1;
        end
    end
 
    
    if lambda_tik
        [u,s,v] = svd(Acs, 'econ');
        
        s_inv = diag(s); 
        s_inv = conj(s_inv) ./ (abs(s_inv).^2 + lambda_tik);
        
        Acs_inv = v * diag(s_inv) * u';
    end
    

    % estimate kernel weights
    
    weights = zeros([prod(kernel_size) * num_chan, num_chan, Ry-1]);

    for r = 1:Ry-1
        disp(['Kernel group : ', num2str(r)])

        for c = 1:num_chan
            
            if ~lambda_tik
                weights(:,c,r) = Acs \ Rhs(:,c,r);
            else
                weights(:,c,r) = Acs_inv * Rhs(:,c,r);
            end
            
        end
    end



    % recon undersampled data

    Weights = permute(weights, [2,1,3]);

    kspace_recon = padarray(kspace_sampled, [pad_size, 0]);


    for ky = Ky_begin : Ry : Ky_end
        for kx = Kx_begin : Kx_end

            data = kspace_recon(kx-kernel_hsize(1):kx+kernel_hsize(1), ky-kernel_hsize(2)*Ry:Ry:ky+kernel_hsize(2)*Ry, :);                
 
            for ry = 1:Ry-1
                 kspace_recon(kx, ky-ry, :) = Weights(:,:,ry) * data(:);
            end

        end
    end
    
    kspace_recon = kspace_recon(1+pad_size(1):end-pad_size(1), 1+pad_size(2):end-pad_size(2), :);

    % subsititute sampled & acs data
    kspace_recon = kspace_recon .* repmat((~mask & ~mask_acs), [1,1,num_chan]) + kspace_sampled .* repmat(~mask_acs, [1,1,num_chan]) + kspace_acs;
 
    Img_recon(:,:,:,pc) = ifft2c(kspace_recon);
   
end    
 

 

% combine coils and cycles  

rmse_grappa_mip = 100 * norm2( mip_full - max(rsos(Img_recon, 3), [], 4) ) / norm2( mip_full );


mosaic(imrotate(rsos(Img_recon, 3), -90), 1, num_cycle, 2, ['Grappa R=', num2str(Ry), ' RMSE=', num2str(rmse_grappa_mip), '%'], genCaxis(rsos(IMG_gcc, 3)))
 



%--------------------------------------------------------------------------
%% Joint Grappa with different sampling
%--------------------------------------------------------------------------
 

Ry = 6;                     % accl factor
Lambda_tik = 3e-3;          % Tikhonov reg parameter
del_step = 3;               % amount of k-space staggering between phase-cycles
kernel_size = [9,3];        % assume odd kernel size


num_acsX = N(1);            % acs size in readout
num_acsY = 32;              % acs size in phase encoding



mask = zeros([N, num_cycle]);
mask_acs = zeros(N);


mask_acs(1+end/2-num_acsX/2:end/2+num_acsX/2, 1+end/2-num_acsY/2:end/2+num_acsY/2) = 1;
Mask_acs = repmat(mask_acs, [1,1,num_cycle]);


% create sampling patterns

del = mod((0:num_cycle-1) * del_step, Ry)


for t = 1:num_cycle
    mask(:,1+del(t):Ry:end,t) = 1;
end

mosaic( squeeze(mask(1,:,:)).', 1, 1, 10, 'sampling patterns' ), set(gcf, 'color', [1,1,1]*0.2)



kernel_hsize = (kernel_size-1)/2;


pad_size = kernel_hsize .* [1,Ry];
N_pad = N + 2 * pad_size; 


% k-space limits for training:
ky_begin = 1 + Ry * kernel_hsize(2);        % first kernel center point that fits acs region 

ky_end = num_acsY - Ry * kernel_hsize(2);   % last kernel center point that fits acs region 

ky_end = ky_end - max(del);                 % make sure other cycles remain within acs


kx_begin = 1 + kernel_hsize(1);             % first kernel center point that fits acs region 

kx_end = num_acsX - kernel_hsize(1);        % last kernel center point that fits acs region 




% k-space limits for recon:
Ky_begin = 1 + Ry * kernel_hsize(2);        % first kernel center point that fits acs region 

Ky_end = N_pad(2) - Ry * kernel_hsize(2);   % last kernel center point that fits acs region 

Ky_end = Ky_end - max(del);                 % make sure data from other images remain in matrix size


Kx_begin = 1 + kernel_hsize(1);             % first kernel center point that fits acs region 

Kx_end = N_pad(1) - kernel_hsize(1);        % last kernel center point that fits acs region 



% count the no of kernels that fit in acs 
ind = 1;

for ky = ky_begin : ky_end
    for kx = kx_begin : kx_end
        ind = ind + 1;        
    end
end

num_ind = ind;

 
Img_Rec = zeros([N, num_chan, num_cycle]);
    

kspace_sampled = zeros(size(Img_Rec));
kspace_acs = zeros(size(Img_Rec));


for pc = 1:num_cycle
    kspace_full = fft2c( IMG_gcc(:,:,:,pc) );

    kspace_sampled(:,:,:,pc) = kspace_full .* repmat(mask(:,:,pc), [1,1,num_chan]);

    kspace_acs(:,:,:,pc) = kspace_full .* repmat(mask_acs, [1,1,num_chan]);  
end



% train kernel
kspace_acs_crop = kspace_acs(1+end/2-num_acsX/2:end/2+num_acsX/2, 1+end/2-num_acsY/2:end/2+num_acsY/2, :, :);



Rhs = zeros([num_ind, num_chan, Ry-1, num_cycle]);
Acs = zeros([num_ind, prod(kernel_size) * num_chan * num_cycle]);


disp(size(Acs))


ind = 1;

for ky = ky_begin : ky_end
    for kx = kx_begin : kx_end

        acs = zeros( prod(kernel_size) * num_chan * num_cycle, 1 );

        for ry = 0 : num_cycle - 1

            tmp = kspace_acs_crop(kx-kernel_hsize(1):kx+kernel_hsize(1), del(ry+1) + ky-kernel_hsize(2)*Ry : Ry : del(ry+1) + ky+kernel_hsize(2)*Ry, :, ry + 1);

            acs( 1 + prod(kernel_size) * num_chan * ry : prod(kernel_size) * num_chan * (ry + 1) ) = tmp(:);

        end

        Acs(ind,:) = acs;


        for phs_cyc = 1:num_cycle
            for ry = 1:Ry-1
            
                Rhs(ind,:,ry,phs_cyc) = kspace_acs_crop(kx, del(phs_cyc) + ky-ry, :, phs_cyc);
                
            end
        end 

        ind = ind + 1;

    end
end

  


if Lambda_tik
    [u,s,v] = svd(Acs, 'econ');

    s_inv = diag(s); 
    s_inv = conj(s_inv) ./ (abs(s_inv).^2 + Lambda_tik);

    Acs_inv = v * diag(s_inv) * u';
end
    

% estimate kernel weights
weights = zeros([prod(kernel_size) * num_chan * num_cycle, num_chan, Ry-1, num_cycle]);

for phs_cyc = 1:num_cycle
    disp(['Phase cycle: ', num2str(phs_cyc)])

    for r = 1:Ry-1
        for c = 1:num_chan

            if ~Lambda_tik
                weights(:,c,r,phs_cyc) = Acs \ Rhs(:,c,r,phs_cyc);
            else
                weights(:,c,r,phs_cyc) = Acs_inv * Rhs(:,c,r,phs_cyc);
            end
        end
    end
end



% recon undersampled data
Weights = permute(weights, [2,1,3,4]);

kspace_recon = padarray(kspace_sampled, [pad_size, 0, 0]);


for ky = Ky_begin : Ry : Ky_end
    for kx = Kx_begin : Kx_end

        data = zeros( prod(kernel_size) * num_chan * num_cycle, 1 );

        for ry = 0:num_cycle-1

            dt = kspace_recon(kx-kernel_hsize(1):kx+kernel_hsize(1), del(ry+1) + ky-kernel_hsize(2)*Ry : Ry : del(ry+1) + ky+kernel_hsize(2)*Ry, :, ry + 1);

            data( 1 + prod(kernel_size) * num_chan * ry : prod(kernel_size) * num_chan * (ry + 1) ) = dt(:);

        end


        for phs_cyc = 1:num_cycle
            for ry = 1:Ry-1
                
                kspace_recon(kx, del(phs_cyc) + ky-ry, :, phs_cyc) = Weights(:,:,ry,phs_cyc) * data(:);
                 
            end
        end

    end
end


kspace_recon = kspace_recon(1+pad_size(1):end-pad_size(1), 1+pad_size(2):end-pad_size(2), :, :);


% subsititute sampled data
kspace_recon = kspace_recon .* repmat( permute( (~mask & ~Mask_acs), [1,2,4,3]), [1,1,num_chan,1]) + ...
    kspace_sampled .* repmat( permute( ~Mask_acs, [1,2,4,3]), [1,1,num_chan,1]) + kspace_acs;


Img_Rec = zeros(size(kspace_recon));

for phs_cyc = 1:num_cycle
    for c = 1:num_chan
        Img_Rec(:,:,c,phs_cyc) = ifft2c( kspace_recon(:,:,c,phs_cyc) );
    end
end



rmse_jgrappa_mip = 100 * norm2( mip_full - max(rsos(Img_Rec, 3), [], 4) ) / norm2( mip_full );


mosaic(imrotate(rsos(Img_Rec, 3), -90), 1, num_cycle, 3, ['Joint Grappa R=', num2str(Ry), ' RMSE=', num2str(rmse_jgrappa_mip), '%'], genCaxis(rsos(IMG_gcc, 3)))



 
%--------------------------------------------------------------------------
%% Joint L1-Spirit with different sampling
%--------------------------------------------------------------------------


Ry = 6;                     % accl factor
del_step = 3;               % amount of k-space staggering between phase-cycles


num_acsX = N(1);            % acs size in readout
num_acsY = 32;              % acs size in PE


mask = zeros([N, num_cycle]);
mask_acs = zeros(N);


mask_acs(1+end/2-num_acsX/2:end/2+num_acsX/2, 1+end/2-num_acsY/2:end/2+num_acsY/2) = 1;
 

% create sampling patterns
del = mod((0:num_cycle-1) * del_step, Ry)       % staggering amount

for t = 1:num_cycle
    mask(:,1+del(t):Ry:end,t) = 1;
end


Mask_spirit = ( repmat(mask_acs, [1,1,num_cycle]) + mask ) ~= 0;
Mask_spirit = repmat(permute(Mask_spirit, [1,2,4,3]), [1,1,num_chan,1]);

DATA = fft2c2( IMG_gcc );


% undersampled data
DATA_sampled = DATA .* Mask_spirit;
myshow3(rsos(ifft2c(DATA_sampled),3))
DATA_Sampled = reshape( DATA_sampled, [N, num_chan * num_cycle] );

 

% Perform Calibration   
kSize = [7,7];              % SPIRiT kernel size
CalibTyk = 8e-3;            % Tykhonov regularization in the calibration

CalibSize = [num_acsX, num_acsY];     

kCalib = crop(DATA_Sampled, [CalibSize, num_chan*num_cycle]);
myshow3(rsos(ifft2c(kCalib),3))

tic
    kernel = calibSPIRiT(kCalib, kSize, num_chan*num_cycle, CalibTyk);
toc


GOP = SPIRiT(kernel, 'fft', N);



% Perform Reconstruction
lambda = 5e-5;              % TV penalty
L = 1;                      % set "L=1" for TV penalty on each coil separately
                            % set "L=12" for joint TV penalty across all coils/cycles

iter_inner = 10;            % outer cg loops
iter_outer = 10;            % inner cg loops

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
    
	res = fnlCg_l1spirit(res, param);
    
    Res(idx_mis) = res;
    
    img_jspirit = reshape(ifft2c(Res), [N, num_chan, num_cycle]);
        
    rmse_jspirit_mip = 100 * norm2( mip_full - max(rsos(img_jspirit, 3), [], 4) ) / norm2(mip_full)
    
    mosaic(imrotate(rsos(img_jspirit, 3), -90), 1, num_cycle, 4, ['Joint L1-SPIRiT R=', num2str(Ry), ' RMSE=', num2str(rmse_jspirit_mip), '%'], genCaxis(rsos(IMG_gcc, 3)))
	      
end
toc
 



%--------------------------------------------------------------------------
%% compare MIP images from 3 different reconstructions
%--------------------------------------------------------------------------

mosaic(imrotate(max(rsos(Img_recon, 3), [], 4), -90), 1, 1, 10, ['GRAPPA R=', num2str(Ry), ' RMSE=', num2str(rmse_grappa_mip), '%'], genCaxis(rsos(IMG_gcc, 3)))
    
mosaic(imrotate(max(rsos(Img_Rec, 3), [], 4), -90), 1, 1, 11, ['Joint GRAPPA R=', num2str(Ry), ' RMSE=', num2str(rmse_jgrappa_mip), '%'], genCaxis(rsos(IMG_gcc, 3)))

mosaic(imrotate(max(rsos(img_jspirit, 3), [], 4), -90), 1, 1, 12, ['Joint L1-SPIRiT R=', num2str(Ry), ' RMSE=', num2str(rmse_jspirit_mip), '%'], genCaxis(rsos(IMG_gcc, 3)))



  