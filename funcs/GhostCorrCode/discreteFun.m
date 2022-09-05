%   This function finds the (kappa, phi) value that minimizes the phase
%   ghost using the a referenceless method. It takes a starting range and
%   returns the solution (kappa, phi). Uses multi-resolution discrete
%   search
%
%   Inputs:
%       kspOrig                 [nRO, nPE] - uncorrected ksp
%       startOpt                fields: .kRange, .pRange, .kCenter,
%                                       .pCenter, .nGrid
%       method                	options: 'ent', 'svd', 'Gh/Ob')
%
%   Outputs:
%       kappa, phi              the estimated solution
%       outputMap               the search space used including kappas,
%                               phis, and cost function 
%
%
%%%  Copyright  Jessica A McKay and Patrick J Bolan, Ph.D. 
%%%     University of Minnesota , September 2018




%%

function [kappaCenter, phiCenter, outputMap] = discreteFun(...
    kspOrig, startOpt, method)

nGrids = length(startOpt.nGrid);
[nRO,nPE,~,~] = size(kspOrig);

kappaRange = startOpt.kRange; 
phiRange = startOpt.pRange; 
kappaCenter = startOpt.kCenter; 
phiCenter = startOpt.pCenter;

% SVD Indeces
if strcmp(method,'svd');
    
    blockX = 3; % kernel size
    blockY = 3; % kernel size

    nBlock = blockX * blockY;
    
    % Precaculate the mapping from kspTmp to the Mat matrix
    idxMap = zeros((nRO-blockX+1) * (nPE-blockY+1), nBlock);
    
    % These matrices define the offsets for the kernel
    kernelOffsetX = repmat(0:blockX-1, [blockY 1])';
    kernelOffsetY = repmat(0:blockY-1, [blockX 1]);
    
    count = 0;
    for kx = 1:nRO-blockX+1;
        for ky = 1:nPE-blockY+1; 
            count = count + 1;
            kernelIndices = sub2ind([nRO,nPE], kernelOffsetX+kx, kernelOffsetY+ky); % JAM
            idxMap(count,:) = kernelIndices(:);
        end
    end
else 
    idxMap = [];
end

for iGrid = 1:nGrids;
    kappaSteps = startOpt.nGrid(iGrid); % number of steps to use
    phiSteps = startOpt.nGrid(iGrid); % number of steps to use
    
    kappaArray = linspace(kappaCenter - kappaRange, kappaCenter + kappaRange, kappaSteps);
    phiArray = linspace(phiCenter - phiRange, phiCenter + phiRange, phiSteps);
    
    % Preallocation:
    mapTmp = zeros(kappaSteps, phiSteps);
    
    for iKappa = 1:kappaSteps;
        kappaTmp = kappaArray(iKappa);
        for iPhi = 1:phiSteps;
            phiTmp = phiArray(iPhi);
            
            kspMod = applyFirstOrderPhaseCorr(kspOrig,kappaTmp,phiTmp);
            
            mapTmp(iKappa,iPhi) = getMetric_v2( kspMod, method, idxMap);
        end
    end
    ind = find(mapTmp == min(min(mapTmp)));
    [kapInd, phiInd] = ind2sub(size(mapTmp),ind);
    
    kappaCenter = mean(kappaArray(kapInd));
    phiCenter = mean(phiArray(phiInd));
    
    kappaRange = kappaRange / (iGrid+1);
    phiRange = phiRange/ (iGrid+1);
    
    outputMap.map{iGrid} = mapTmp; 
    outputMap.kappasUsed{iGrid} = kappaArray; 
    outputMap.phisUsed{iGrid} = phiArray; 
    
end

end






