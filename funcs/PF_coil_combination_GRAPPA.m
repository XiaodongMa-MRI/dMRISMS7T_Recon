function complex_data = PF_coil_combination_GRAPPA(temp_kespace,par)            
    dkspace1_all_tmp = permute(temp_kespace,[1 2 4 3 5]);
    sz_dkspace1_all_tmp = size(dkspace1_all_tmp);
%     ksp_PFrecon = zeros([200 200 par.nCH sz_dkspace1_all_tmp(4) sz_dkspace1_all_tmp(5)]);
    ksp_PFrecon = zeros([size(dkspace1_all_tmp,1) size(dkspace1_all_tmp,1) par.nCH sz_dkspace1_all_tmp(4) sz_dkspace1_all_tmp(5)]);
    for b_flag = 1:sz_dkspace1_all_tmp(5)
        disp(['PF Recon, nb = ',num2str(b_flag)])
        ksp_PFrecon(:,:,:,:,b_flag) = MB_partial_fourier_dep_v2(squeeze(dkspace1_all_tmp(:,:,:,:,b_flag)),par);
    end

           
    %% cpx_output channel combination
    img_PFrecon = ifftc(ifft2c(ksp_PFrecon),4);
    img_PFrecon = permute(img_PFrecon,[1 2 4 3 5]); % ref in [x y z ch]

    complex_data = coilCombine3D( img_PFrecon ); %[x y z ch]
end


