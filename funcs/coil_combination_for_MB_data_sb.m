function [ref0_GCC,single_band_ref_GCC,kspace0_GCC] = coil_combination_for_MB_data_sb(ref0,single_band_ref,kspace0,ncc)
    Calib_size = [56 56];% the calibration matrix to calculate compression matrix
    Method = 1;% compression method option: 0->SCC; 1->GCC 
   
    [f0,~,MB,c0] = size(ref0);
    if Method == 0  % SCC
%         dim = 0;
%         calibkdata = crop(squeeze(ref0(:,:,:,1)),Calib_size(1),Calib_size(2),(c0));
%         sccmtx = calcSCCMtx(calibkdata);
%         ccmtx = sccmtx(:,1:ncc);

    else % GCC
        slwin = 9;
        dim = 1;
        calib = crop(ref0,[f0 Calib_size(2) MB c0]);
        gccmtx = calcGCCMtx(calib,dim,slwin);
        ccmtx = alignCCMtx(gccmtx(:,1:ncc,:));
    end
    
    %% Coil compression
    for mb=1:MB % Reference
        ref0_GCC(:,:,mb,:)=CC(ref0(:,:,mb,:),ccmtx,dim);
        single_band_ref_GCC(:,:,mb,:)=CC(single_band_ref(:,:,mb,:),ccmtx,dim);
        for nbi = 1:size(kspace0,5)
            kspace0_GCC(:,:,mb,:,nbi)=CC(kspace0(:,:,mb,:,nbi),ccmtx,dim);
        end
    end


