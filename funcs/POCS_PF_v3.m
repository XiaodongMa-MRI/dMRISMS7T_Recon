function res=POCS_PF_v3(ksp,ky,acqky)
% v3 daiep 160816
% First ifft in kx direction, then PF recon in ky direction.

% zhangzhe 2015.7.10
% input: kspace, zero filling with non acquired part, [kx,ky,..]
%       number of max ky, should be 2nd dim of input kspace
%       acqky: array of acquired ky lines, for example [32:128] when ky=128
% output:reconed kspace

if max(acqky)>max(ky)
    	error('check input')
else
    [nX,nY,n1,n2,n3,n4]=size(ksp);
    if nY~=ky
        error('check input,k-space should be alined by [ky,ky,...]')
    else
        N=n1*n2*n3*n4;
        ksp=reshape(ksp,[nX,nY,N]);
        ksp=ifftc(ksp,1);% daiep 160816
        res=single(zeros(size(ksp)));
        for n=1:N 
            hannwin=hann(18);
            win_width=2*length(acqky)-ky;
            filt0=ones(win_width,1);
            filt0(1:8)=hannwin(2:9);
            filt0(end-8+1:end)=hannwin(10:17);
            filt0=zpad(filt0,[nY,1]);
            filt0=repmat(filt0',[nX,1]);
            % phi=ifft2c(ksp(:,:,n).*filt0);
            phi=ifftc(ksp(:,:,n).*filt0,2);
            phi=phi./abs(phi); %œ‡Œªπ¿º∆

            if acqky(1)~=1 % 011 PF style
                filt1=filt0;
                filt1(:,round(nY/2):nY)=1;
                % I0=ifft2c(ksp(:,:,n).*filt1);
                I0=ifftc(ksp(:,:,n).*filt1,2);

                for iter=1:4
                    I1=abs(I0).*phi;
                    % k1=fft2c(I1);
                    k1=fftc(I1,2);
                    k2=ksp(:,:,n);
                    k2(:,1:acqky(1)-1)=k1(:,1:acqky(1)-1);
                    filt2=repmat(hannwin',[nX,1]);
                    k2(:,acqky(1):acqky(1)+8)=k1(:,acqky(1):acqky(1)+8).*filt2(:,10:18)+ksp(:,acqky(1):acqky(1)+8).*filt2(:,1:9);% daiep 160816
                    % I0=ifft2c(k2);
                    I0=ifftc(k2,2);
                end
                
            else % 110 PF style
                filt1=filt0;
                filt1(:,1:round(nY/2))=1;
                % I0=ifft2c(ksp(:,:,n).*filt1);
                I0=ifftc(ksp(:,:,n).*filt1,2);

                for iter=1:6
                    I1=abs(I0).*phi;
                    % k1=fft2c(I1);
                    k1=fftc(I1,2);
                    k2=ksp(:,:,n);
                    k2(:,acqky(end)+1:end)=k1(:,acqky(end)+1:end);
                    filt2=repmat(hannwin',[nX,1]);
                    k2(:,acqky(end)-8:acqky(end))=k1(:,acqky(end)-8:acqky(end)).*filt2(:,1:9)+ksp(:,acqky(end)-8:acqky(end)).*filt2(:,10:18);% daiep 160816
                    % I0=ifft2c(k2);
                    I0=ifftc(k2,2);
                end
            end
            
            res(:,:,n)=k2;
        end
        res=fftc(res,1);% daiep 160816
        
    end
end
res=reshape(res,[nX,nY,n1,n2,n3,n4]);
end
