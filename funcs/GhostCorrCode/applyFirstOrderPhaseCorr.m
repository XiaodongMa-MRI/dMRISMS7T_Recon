%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify an EPI K-space according to constant and linear order phase
% Assumes first dimension is RO, second dimension is PE
%
%  Inputs: 
%       ksp     [nRO,nPE] - uncorrected k-space
%       kappa   shift of each line from ideal. The total shift between even
%               and odd is 2*kappa. Units are pixel.
%       phi     phase modulation of each line from ideal. Total pahse 
%               difference between even and odd is 2*phi. Units are radians
%      
%%%  Copyright  Jessica A McKay and Patrick J Bolan, Ph.D. 
%%%     University of Minnesota , September 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function kspmod = applyFirstOrderPhaseCorr(ksp, kappa, phi)


if mod(size(ksp,2),2); % If odd
    kspmod = zeros(size(ksp,1),size(ksp,2)+1);
end
    
kspmod(1:size(ksp,1),1:size(ksp,2),1:size(ksp,3)) = ksp;

kspmod = double(kspmod);
% img = fftshift(fft2(ifftshift(ksp)));

% First apply phase correction
% kspmod(:,1:2:end) = kspmod(:,1:2:end) .* exp(+1j*phi);
% kspmod(:,2:2:end) = kspmod(:,2:2:end) .* exp(-1j*phi);


%% 

% Then perform sub-pixel shifting
% I had done this with floatingCircShift, but this is incorrect, leading to
% a cosine modulation in image space.
%kspmod(:,1:2:end) = floatingCircShift(kspmod(:,1:2:end),[-kappa 0]); % shift
%kspmod(:,2:2:end) = floatingCircShift(kspmod(:,2:2:end),[kappa 0]); % shift

% Instead, FT to image space in RO, apply phase correction
phase_correction = zeros(size(kspmod));
[Nro, Npe, Ncoil] = size(kspmod);
% xunits = (1:Nro)  - 1;
xunits = (1:Nro) - (Nro/2) - 0.5;

% Specify the even and odd traces
phase_correction(:,1) = exp(+1j * (2*pi*kappa * xunits./(Nro)+phi)); 
phase_correction(:,2) = exp(-1j * (2*pi*kappa * xunits./(Nro)+phi)); 

% Replicate them over the full matrix 
phase_correction = repmat(phase_correction(:,1:2), [1, Npe/2]);
phase_correction = repmat(phase_correction, [1, 1, Ncoil]);

% Apply phase correction in partial image domain
% ftData is the FID, fourier transformed ONLY in the readout direction
% ftData = fftshift(fft(ifftshift(kspmod,1)),1);
% kspmod = ifftshift(ifft(ifftshift(ftData.* phase_correction,1),[],1),1);


%pzy
ftData = fftshift(ifft(ifftshift(kspmod,1)),1);
kspmod = fftshift(fft(ifftshift(ftData.* phase_correction,1),[],1),1);