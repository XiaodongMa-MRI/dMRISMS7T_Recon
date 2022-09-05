% zijing original
% daiep 160105. Perform coil compression for SMS
% xiaodong change
function [b0ksp2,b1ksp2,b0ksp2_nav,b1ksp2_nav,par] = CoilCompression_mxd(b0ksp2,b1ksp2,b0ksp2_nav,b1ksp2_nav,par,Par_SYM);
%%
ncc = Par_SYM.ncc;
% Calib_ysize = 32;
%% Calculate the coil compression matrix
calib_data = squeeze(mean(b0ksp2,4));
[f0_ref,p0_ref,c0,s0]=size(calib_data);
Calib_ysize = 32;
% if Method==0  % SCC
%     ccmtx = zeros(c0,ncc,s0);
%     dim=0;
%     for s_i=1:s0
%         calibkdata=crop(squeeze(calib_data(:,:,:,s_i)),f0_ref,double(Calib_ysize),double(c0));
%         sccmtx=calcSCCMtx(calibkdata);
%         ccmtx(:,:,s_i)=sccmtx(:,1:ncc);
%     end
% else % GCC
    ccmtx = zeros(c0,ncc,f0_ref,s0);
    slwin=5;
    dim=1;
    for s_i=1:s0
        gccmtx=calcGCCMtx(crop(calib_data(:,:,:,s_i),f0_ref,double(Calib_ysize),double(c0)),dim,slwin);
        % alignment
        gccmtx_aligned=alignCCMtx(gccmtx(:,1:ncc,:));
        ccmtx(:,:,:,s_i)=gccmtx_aligned;
    end
% end
%% Coil compression
[f0, p0, c0, Nex_b0, s0] = size(b0ksp2);
f0_nav = size(b0ksp2_nav,1);
p0_nav = size(b0ksp2_nav,2);
ind1 = round( ( (1:f0_nav)-1 )*f0_ref/f0_nav+1 );
ccmtx1 = ccmtx(:,:,ind1,:);

b0ksp2_cc = zeros(f0, p0, ncc, par.nNEX0, s0);
b1ksp2_cc = zeros(f0, p0, ncc, par.nB, par.nNEX, s0);
b0ksp2_nav_cc = zeros(f0_nav, p0_nav, ncc, par.nSHOT, par.nNEX0, s0);
b1ksp2_nav_cc = zeros(f0_nav, p0_nav, ncc, par.nSHOT, par.nB, par.nNEX, s0);
for s_i=1:s0
    for nex0=1:par.nNEX0
        b0ksp2_cc(:,:,:,nex0,s_i)=CC(b0ksp2(:,:,:,nex0,s_i),ccmtx(:,:,:,s_i),dim);
        for shot=1:par.nSHOT
            b0ksp2_nav_cc(:,:,:,shot,nex0,s_i)=CC(b0ksp2_nav(:,:,:,shot,nex0,s_i),ccmtx1(:,:,:,s_i),dim);
        end
    end
    
    for nex=1:par.nNEX
        for nb=1:par.nB
            b1ksp2_cc(:,:,:,nb,nex,s_i)=CC(b1ksp2(:,:,:,nb,nex,s_i),ccmtx(:,:,:,s_i),dim);
            for shot=1:par.nSHOT
                b1ksp2_nav_cc(:,:,:,shot,nb,nex,s_i)=CC(b1ksp2_nav(:,:,:,shot,nb,nex,s_i),ccmtx1(:,:,:,s_i),dim);
            end
        end
    end
end

%% output
b0ksp2 = single(b0ksp2_cc);
b1ksp2 = single(b1ksp2_cc);
b0ksp2_nav = single(b0ksp2_nav_cc);
b1ksp2_nav = single(b1ksp2_nav_cc);
par.nCH = ncc;