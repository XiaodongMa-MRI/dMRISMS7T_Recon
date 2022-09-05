%% Func: recon partial fourier data
% Input: b0ksp2, f0*p0*c0*s0; par, parameter matrix
% Input: b0ksp2, f0*p0*s0 ---modified by Simin June 2018
% Output: recovered full kspace
% Note: for half fourier matraix, no zero is allowed in the sampled matrix.
% 150202, daiep, use POCS_PF instead of pocs_mute in case of possible
% scaling

function [b0ksp3]=MB_partial_fourier_dep_v2(b0ksp2,par)
for ch=1:par.nCH
    for n=1:par.nSL
        ksptemp=single(zeros(par.kx_r,par.ky_r));
        ksptemp(:,par.ky_r-par.ky_s+1:par.ky_r)=b0ksp2(:,:,ch,n);
        b0ksp3(:,:,ch,n)=single(POCS_PF_v3(ksptemp,par.ky_r,[par.ky_r-par.ky_s+1:par.ky_r]));%daiep 160816
    end
end
end
