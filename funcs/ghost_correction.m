function [kspCorrRefLessSimplex,kappaEstSimplex, phiEstSimplex] = ghost_correction(kspGhost,method,kRange,pRange,kCenter,pCenter)
    %kspGhost is a 2D-kspace; the first dimension is the PE direction
%     method = 'Gh/Ob'; % options: 'ent','svd','Gh/Ob'

    startOptDiscrete.kRange = kRange; 
    startOptDiscrete.pRange = pRange;
    startOptDiscrete.kCenter = kCenter; % May need to check the polarity here!!!
    startOptDiscrete.pCenter = pCenter;
    startOptDiscrete.nGrid = [11, 11, 11]; % This means use three grids of 11x11 each
    % ^ I recommend using one larger grid at first so you can see the whole
    % metric space for debugging

%     fprintf('Starting DISCRETE SEARCH \r');

    kspOrigTmp = squeeze(kspGhost);
    [kappaEstDiscrete, phiEstDiscrete, outputMap] = discreteFun(kspOrigTmp, startOptDiscrete, method);


    % % Option 2 - use the solutions from discrete search:
    startOptSimplex.kCenter = kappaEstDiscrete; % you need to get the polarity right on this
    startOptSimplex.pCenter = phiEstDiscrete; % In my data phi ~= 0 

    [kappaEstSimplex, phiEstSimplex] = simplexFun(kspGhost, startOptSimplex, method);



%% Now apply them and see results

%     % Apply the first guess from the discrete search
%     kspCorrRefLessDiscrete = applyFirstOrderPhaseCorr( kspGhost, kappaEstDiscrete, phiEstDiscrete);   
%     imgRefLessDiscrete = fftshift(fft2(ifftshift( kspCorrRefLessDiscrete ))); 

    % Apply the refined guess via simplex search
    kspCorrRefLessSimplex = applyFirstOrderPhaseCorr( kspGhost, kappaEstSimplex, phiEstSimplex); 
%     imgRefLessSimplex = fftshift(fft2(ifftshift( kspCorrRefLessSimplex ))); 