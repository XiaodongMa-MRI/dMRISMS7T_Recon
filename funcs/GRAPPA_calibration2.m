function [P1, P2]=GRAPPA_calibration2(ref0, MB, R_PE, ker_x, ker_y, lamb);
vir_R=R_PE*MB;

[f1, acs_y, ~, ~]=size(ref0);
ker_y_off=(ker_y-R_PE+1)/2-1;
% the kernel is symmetric about the missing data in ky direction
ker_x_off = floor(ker_x/2);
P1=[];
P2=[];
for mb=1:MB    
    %% step1
    if MB>1
        MB_ker_y=3;
        k_acs0 = [];
        k_block0=[];
        for i_y=1:acs_y-(R_PE*MB_ker_y-1)
            k_acs = [];
            k_block=[];
            tmp=[];
            for xx=1:ker_x
                for yy=1:R_PE:R_PE*MB_ker_y
                    tmp=ref0(xx:f1-(ker_x-xx),i_y+yy-1,mod(mb+(yy-1)/R_PE-1,MB)+1,:);
                    k_block=[k_block,reshape(tmp,f1-ker_x+1,[])];
                end
            end
            k_block0=[k_block0;k_block];
            
            ind_z=mod(mb+floor(MB_ker_y/2)-1,MB)+1;
            k_acs=ref0((ker_x+1)/2:f1-(ker_x-1)/2,i_y+R_PE*floor(MB_ker_y/2),[1:ind_z-1, ind_z+1:MB],:);
            k_acs0=[k_acs0;reshape(k_acs,f1-ker_x+1,[])];
        end
        I=ones(1,size(k_block,2));
        I=diag(I);
        AtA = k_block0'*k_block0;
        lamb_new = norm(AtA,'fro')/size(AtA,1)*lamb;
        P1(:,:,mb)=(k_block0'*k_block0+lamb_new*I) \k_block0'* k_acs0;
        
    end
end

for mb=1:MB    
    %% step 2
    if R_PE>1
        k_acs0 = [];
        k_block0=[];
        for i_y=1:acs_y-(ker_y-1)
            k_acs = [];
            k_block=[];
            tmp=[];
            for xx=1:ker_x
                for yy=1:R_PE:ker_y
                    tmp=ref0(xx:f1-(ker_x-xx),i_y+yy-1,:,:);
                    k_block=[k_block,reshape(tmp,f1-ker_x+1,[])];
                end
            end
            k_block0=[k_block0;k_block];
            k_acs=ref0((ker_x+1)/2:f1-(ker_x-1)/2,i_y+ker_y_off+1:i_y+ker_y_off+R_PE-1,:,:);
            k_acs0=[k_acs0;reshape(k_acs,f1-ker_x+1,[])];
        end
        I=ones(1,size(k_block,2));
        I=diag(I);
        AtA = k_block0'*k_block0;
        lamb_new = norm(AtA,'fro')/size(AtA,1)*lamb;
        P2(:,:,mb)=(k_block0'*k_block0+lamb_new*I) \k_block0'* k_acs0;
    end
end