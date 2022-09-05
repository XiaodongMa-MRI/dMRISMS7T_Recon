function complex_data =  PF_coil_combination_L1SPIRiT(dkspace1_tmp,par,move_flag1)
            
    tempimg = ifftc(ifft2c(dkspace1_tmp),4);
    tempimgbp = tempimg;
    if move_flag1
        tempimg(:,1:move_flag1,:,1) = tempimgbp(:,(end-move_flag1+1):end,:,1);
        tempimg(:,(move_flag1+1):end,:,1) = tempimgbp(:,1:end-move_flag1,:,1);
        tempimg(:,(end-move_flag1+1):end,:,3) = tempimgbp(:,1:move_flag1,:,3); %why move?
        tempimg(:,1:end-move_flag1,:,3) = tempimgbp(:,(move_flag1+1):end,:,3);
    end

    %% PF recon
    dkspace1_all_tmp = permute(fftc(fft2c(tempimg),4),[1 2 3 4]);
    ksp_PFrecon = zeros([200 200 par.nCH par.nSL]);

    for b_flag = 1
        disp(['PF Recon, nb = ',num2str(b_flag)])
        ksp_PFrecon(:,:,:,:,b_flag) = MB_partial_fourier_dep_v2(squeeze(dkspace1_all_tmp(:,:,:,:,b_flag)),par);
    end

    img_PFrecon = ifftc(ifft2c(ksp_PFrecon),4);

    %% cpx_output channel combination
    img_PFrecon = permute(img_PFrecon,[1 2 4 3]); % ref in [x y z ch]
    complex_data = coilCombine3D( img_PFrecon ); %[x y z ch]
    
end
