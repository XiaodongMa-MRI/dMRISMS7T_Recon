%   This function calculates the cost function for various referenceless
%   ghost correction methods

% Note: This version does NOT handel undersampled data -- may need to
% modify for R>1

%   Inputs:
%       kspOrig                 [nRO, nPE] - uncorrected ksp
%       method                	options: 'ent', 'svd', 'Gh/Ob')
%       idxMap                  for SVD method (saves time to caculate once
%                               and pass it in). Use [] for other methods.  
%
%   Outputs:
%       metricOut               the calculated value for the cost function
%                               *This is what you want to minimize 
%
%
%%%  Copyright  Jessica A McKay and Patrick J Bolan, Ph.D. 
%%%     University of Minnesota , September 2018


%%

function [metricOut] = getMetric_v2(kspOrig,method,idxMap); 

switch method
    
    case 'entSmooth'
%         imgTmp = fftshift(fft2(ifftshift(kspOrig)));
        imgTmp = rsos(fftshift(ifft2(ifftshift(kspOrig))),3);
        imgSmooth = medfilt2(abs(imgTmp));
        metricOut = entropy(mat2gray(double(abs(imgSmooth))));
    
    case 'ent'
        if size(kspOrig,3)==1
            imgTmp = fftshift(fft2(ifftshift(kspOrig)));
        else
            imgTmp = rsos(fftshift(ifft2(ifftshift(kspOrig))),3);
        end
%         metricOut = entropy(mat2gray(double((imgTmp))));
        metricOut = entropy(mat2gray(double(abs(imgTmp))));

    case 'svd'
        minZone = 2; % first index of min zone
%         minZone = size(idxMap,2)/3; % first index of min zone -- JAM HACK 
        Mat = kspOrig(idxMap);
        S = svd(Mat);
        tail = S(minZone:end);
%         tail = S(minZone+1:end); % HACK
        metricOut = sum(tail); % The integratl of the "tail" of the SVD vector
        
    case 'GhOb'
        imgTmp = rsos(fftshift(ifft2(ifftshift(kspOrig))),3);
        shift = circshift(imgTmp, [0,size(imgTmp,2)/2]);
        
        metGhOb = abs(imgTmp./shift); % signal/shifted signal
        metGhOb = medfilt2(metGhOb);
        metric = nanmean(metGhOb(:));
        metricOut = 1./metric; % Invert for minizmation problem
        
    case 'energy'
        imgTmp = rsos(fftshift(ifft2(ifftshift(kspOrig))),3);
        tmp = graycoprops(im2uint8(mat2gray(abs(imgTmp))),'energy');
        metricOut = 1./ (tmp.Energy); 
end
