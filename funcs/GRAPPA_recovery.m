%%

%% setting Parameter
function [dkspace0]=GRAPPA_recovery(P1,P2,kspace0,MB, R_PE, ker_x, ker_y, lamb);
vir_R=R_PE*MB;

[f0, p0, z0, c0]=size(kspace0);
dkspace0=kspace0;
ker_y_off=(ker_y-R_PE+1)/2-1;
% the kernel is symmetric about the missing data in ky direction
ker_x_off = floor(ker_x/2);

for mb=1:MB    
    dkspace_t=zeros([f0+ker_x-1, p0+ker_y-1, z0, c0]);
    dkspace_t(ker_x_off+1:ker_x_off+f0,1:p0,:,:)=dkspace0;
    %% step1
    if MB>1
       MB_ker_y=5;
       ind_z=mod(mb+floor(MB_ker_y/2)-1,MB)+1;
        for i_y=(mb-1)*R_PE+1:vir_R:p0-(R_PE*MB_ker_y-1)
            k_block=[];
            bb=[];
            tmp=[];
            for xx=1:ker_x
                for yy=1:R_PE:R_PE*MB_ker_y
                    tmp=dkspace_t(xx:f0+xx-1,i_y+yy-1,mod(mb+(yy-1)/R_PE-1,MB)+1,:);
                    k_block=[k_block,reshape(tmp,f0,[])];
                end
            end
            bb =k_block * P1(:,:,mb);
            bb=reshape(bb,[f0,1,MB-1,c0]);
            % dkspace_t(ker_x_off+1:ker_x_off+f0,i_y+R_PE,[1:ind_z-1, ind_z+1:MB],:)=bb;
            dkspace0(1:f0,i_y+R_PE*floor(MB_ker_y/2),[1:ind_z-1, ind_z+1:MB],:)=bb;
        end
    end
end

for mb=1:MB    
    dkspace_t=zeros([f0+ker_x-1, p0+ker_y-1, z0, c0]);
    dkspace_t(ker_x_off+1:ker_x_off+f0,1:p0,:,:)=dkspace0;
    %% step 2
    if R_PE>1
        for i_y=(mb-1)*R_PE+1:vir_R:p0-floor(ker_y/2)
            k_block=[];
            bb=[];
            tmp=[];
            for xx=1:ker_x
                for yy=1:R_PE:ker_y
                    tmp=dkspace_t(xx:f0+xx-1,i_y+yy-1,:,:);
                    k_block=[k_block,reshape(tmp,f0,[])];
                end
            end
            bb =k_block * P2(:,:,mb);
            bb=reshape(bb,[f0,R_PE-1,MB,c0]);
            dkspace0(1:f0,i_y+ker_y_off+1:i_y+ker_y_off+R_PE-1,:,:)=bb;
        end
    end
end