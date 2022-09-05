%   This function finds the (kappa, phi) value that minimizes the phase
%   ghost using the Entropy method. It takes a starting points and
%   returns the solution (kappa, phi). Uses simplex search
%
%
%   Inputs:
%       kspOrig                 [nRO, nPE] - uncorrected ksp
%       startOpt                fields: .kCenter, .pCenter
%       method                	options: 'ent', 'svd', 'Gh/Ob')
%
%   Outputs:
%       kappa, phi              the estimated solution
%
%
%%%  Copyright  Jessica A McKay and Patrick J Bolan, Ph.D. 
%%%     University of Minnesota , September 2018




%%

function [kappa, phi] = simplexFun(kspOrig, startOpt, method)

[nRO,nPE,~,~] = size(kspOrig);

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



% Create an anonymous function that you can use to pass parameters to the
% optimization function
anonFun = @(x)calculateMetric(x, kspOrig, method, idxMap);

% Starting values
x0(1) = startOpt.kCenter;
x0(2) = startOpt.pCenter;

% Specify fitting options
options = optimset(...
    'Display', 'off', ...
    'TolFun', 1E-6, 'TolX', 1E-6, ...
    'MaxFunEvals', 2000);

% Begin Fitting
% start = tic;
% [x, fval, exitflag, output] = fminsearch(anonFun, x0, options);
[x] = fminsearch(anonFun, x0, options); % JAM
% dur = toc(start);
% fprintf('Completed in %.1f s\n', dur);

% Save outputs
kappa = x(1);
phi = x(2);

return;



% This calculates the Entropy-based metric, that needs to be minimized
% The optimization parameters kappa and phi are arrayed in the parameter
% matrix "x"
    function metric = calculateMetric(x, ksp, method, idxMap)
        kappaTmp = x(1);
        phiTmp = x(2);
                
        % Apply these corrections
        kspMod = applyFirstOrderPhaseCorr(ksp,kappaTmp,phiTmp);
        metric = getMetric_v2(kspMod, method, idxMap); % JAM
    end


end






