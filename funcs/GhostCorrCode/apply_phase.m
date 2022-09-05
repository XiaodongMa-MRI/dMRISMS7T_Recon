function [ksp_res]=apply_phase(b0k,phase,nSHOT)
    [nFE,nPE,nCH,nNEX]=size(b0k);
    nEX=nPE/nSHOT;
    oddk=zeros(size(b0k));
    evenk=zeros(size(b0k));
    for i=1:nEX
        if mod(i,2)==1
            oddk(:,(1:nSHOT)+(i-1)*nSHOT,:,:,:)=b0k(:,(1:nSHOT)+(i-1)*nSHOT,:,:);
        end
        if mod(i,2)==0
            evenk(:,(1:nSHOT)+(i-1)*nSHOT,:,:,:)=b0k(:,(1:nSHOT)+(i-1)*nSHOT,:,:);
        end
    end

    oddxky=ifftshift((fftshift(ifft(ifftshift(oddk,1),[],1),1).*exp(1i*repmat(phase,[1,1,1,nNEX]))));
    evenxky=ifftshift((fftshift(ifft(ifftshift(evenk,1),[],1),1).*exp(-1i*repmat(phase,[1,1,1,nNEX]))));
    im=fftshift(ifft(oddxky+evenxky,[],2));
    ksp_res=fftshift(fftshift(fft(fft(ifftshift(ifftshift(im,1),2),[],1),[],2),1),2);
end
