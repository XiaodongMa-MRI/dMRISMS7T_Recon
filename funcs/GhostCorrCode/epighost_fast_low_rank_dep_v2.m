% function [phase, E0, E1]=epighost_fast_low_rank_dep(b0k,nSHOT)
function [phase, E0, E1]=epighost_fast_low_rank_dep_v2(b0k,nSHOT,nav_meas_sign)
% Zhang Zhe 2015.9.16
% low rank for EPI ghost corr
% ref: Proc. Intl. Soc. Mag. Reson. Med. 23 (2015) 0075
% Acquisition-free Nyquist ghost correction for parallel imaging accelerated EPI
% Based on epighost_fast_low_rank_v2, correct possible aliasing in PE
% directino for navigator.

%SAKE paramater
[nFE,nPE,nCH]=size(b0k);
nEX=nPE/nSHOT;
kSize=[3,4];
nFE_r=round(nFE/4);
nPE_r=nSHOT*12;
if nSHOT==1
%     kSize=[3,6];
    if nPE>48 % Indicate ss-EPI data is not navigator
        nPE_r=48;
        nFE_r=nFE/2;
    else % ss-EPI is likely to be a navigator
        nPE_r=nPE;
        nFE_r=nFE;
    end
end
wnRank=1.8;
%% First search
% phiCON=linspace(-pi/2,pi/2,27);% original
% phiLIN=linspace(-3,3,41);% original
% if nSHOT == 1 && exist('nav_meas_sign','var')
if nSHOT == 1
    phiCON=linspace(-pi/4,pi/4,29);% v2
    phiLIN=linspace(-1.5,1.5,25);% v2
else
    phiCON=linspace(-pi/2,pi/2,29);% The total number should be an odd number. For philips
    phiLIN=linspace(-1.5,1.5,13);
end
% phiLIN=linspace(-3,3,11);% The envelope change smoothly, so loop twice.
% phiCON=phiCON(2:end);

oddk=zeros(size(b0k));
evenk=zeros(size(b0k));
for i=1:nEX
    if mod(i,2)==1
        oddk(:,[1:nSHOT]+(i-1)*nSHOT,:)=b0k(:,[1:nSHOT]+(i-1)*nSHOT,:);
    end
    if mod(i,2)==0
        evenk(:,[1:nSHOT]+(i-1)*nSHOT,:)=b0k(:,[1:nSHOT]+(i-1)*nSHOT,:);
    end
end
% E is for minimization
E0=zeros(length(phiCON),length(phiLIN));
E1=[];
hwait=waitbar(0,'ghost correction using SVD...');
for i=1:length(phiCON)
    waitbar(i/length(phiCON),hwait,'ghost correction using SVD...');
    tic
    for j=1:length(phiLIN)
        % phi=repmat((phiCON(i)+2*pi*(0:(nFE-1))/nFE*phiLIN(j))',[1,nPE,nCH]);
        phi=repmat((phiCON(i)+2*pi*((0:(nFE-1))-floor(nFE/2))/nFE*phiLIN(j))',[1,nPE,nCH]); % daiep 151202
        oddxky=fftshift(ifft(ifftshift(oddk,1),[],1),1).*exp(1i*phi);
        evenxky=fftshift(ifft(ifftshift(evenk,1),[],1),1).*exp(-1i*phi);
        ksptemp=crop(fftshift(fft(ifftshift(oddxky+evenxky,1),[],1),1),nFE_r,nPE_r,nCH);
        % ksptemp=crop(fftshift(fft(ifftshift(oddxky+evenxky,1),[],1),1),nFE_r,nPE,nCH);
        tmp = im2row(ksptemp,kSize);
        tmp = reshape(tmp,size(tmp,1),size(tmp,2)*size(tmp,3));
        [~,S,~] = svd(tmp,'econ');
        s=diag(S);
        % E(i,j)=sum(s(round(length(s)/5):end));
        E0(i,j)=sum(s(round(wnRank*prod(kSize)):end));
    end
    toc
end
[~, minind] = min(E0(:));
[x,y] = ind2sub(size(E0),minind);
close(hwait);
% phase=repmat((phiCON(x)+2*pi*(0:(nFE-1))/nFE*phiLIN(y))',[1,nPE,nCH]);

%% Second search
if y==1
    phiLIN=linspace(phiLIN(y),phiLIN(y+1),5);
elseif y==size(phiLIN,2)
    phiLIN=linspace(phiLIN(y-1),phiLIN(y),5);
else
    phiLIN=linspace(phiLIN(y-1),phiLIN(y+1),9);% The envelope change smoothly, so loop twice.
end
% E is for minimization
E1=zeros(length(phiCON),length(phiLIN));
hwait=waitbar(0,'ghost correction using SVD...');
for i=1:length(phiCON)
    waitbar(i/length(phiCON),hwait,'ghost correction using SVD...');
    for j=1:length(phiLIN)
        % phi=repmat((phiCON(i)+2*pi*(0:(nFE-1))/nFE*phiLIN(j))',[1,nPE,nCH]);
        phi=repmat((phiCON(i)+2*pi*((0:(nFE-1))-floor(nFE/2))/nFE*phiLIN(j))',[1,nPE,nCH]); % daiep 151202
        oddxky=fftshift(ifft(ifftshift(oddk,1),[],1),1).*exp(1i*phi);
        evenxky=fftshift(ifft(ifftshift(evenk,1),[],1),1).*exp(-1i*phi);
        ksptemp=crop(fftshift(fft(ifftshift(oddxky+evenxky,1),[],1),1),nFE_r,nPE_r,nCH);
        tmp = im2row(ksptemp,kSize);
        tmp = reshape(tmp,size(tmp,1),size(tmp,2)*size(tmp,3));
        [~,S,~] = svd(tmp,'econ');
        s=diag(S);
        E1(i,j)=sum(s(round(wnRank*prod(kSize)):end));
    end
end
[~, minind] = min(E1(:));
[x,y] = ind2sub(size(E1),minind);
close(hwait);
% phase=repmat((phiCON(x)+2*pi*(0:(nFE-1))/nFE*phiLIN(y))',[1,nPE,nCH]);
phase=repmat((phiCON(x)+2*pi*((0:(nFE-1))-floor(nFE/2))/nFE*phiLIN(y))',[1,nPE,nCH]); % daiep 151202
%% Show results
oddk=zeros(size(b0k));
evenk=zeros(size(b0k));
for i=1:nEX
    if mod(i,2)==1
        oddk(:,[1:nSHOT]+(i-1)*nSHOT,:)=b0k(:,[1:nSHOT]+(i-1)*nSHOT,:);
    end
    if mod(i,2)==0
        evenk(:,[1:nSHOT]+(i-1)*nSHOT,:)=b0k(:,[1:nSHOT]+(i-1)*nSHOT,:);
    end
end

oddxky=ifftshift((fftshift(ifft(ifftshift(oddk,1),[],1),1).*exp(1i*phase)));
evenxky=ifftshift((fftshift(ifft(ifftshift(evenk,1),[],1),1).*exp(-1i*phase)));
if nCH>1
    myshow3(sos(fftshift(ifft(oddxky+evenxky,[],2))))
end
if nCH==1
    myshow3((fftshift(ifft(oddxky+evenxky,[],2))))
end

end