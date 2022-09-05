
%  runScript_example.m
%  Calls simplexFun.m, discreteFun.m, and applyFirstOrderPhaseCorr.m
%  Load kspGhost (or simulate below):
%%% kspGhost    [nRO,nPE] ~ this is your kspace data with a 1st order Nyquist ghost
%
% Need to set: method, startOptDiscrete, and startOptSimplex, 
%
%%%  Copyright  Jessica A McKay and Patrick J Bolan, Ph.D. 
%%%     University of Minnesota , May 2019
%%%
%%%  Reference:  McKay, JA, Moeller, S, Zhang, L, Auerbach, EJ, Nelson, MT,
%%%              Bolan, PJ. Nyquist ghost correction of breast diffusion
%%%              weighted imaging using referenceless methods. Magn Reson
%%%              Med. 2019; 81: 2624–2631. https://doi.org/10.1002/mrm.27563

%% Start with Read in
 

clear all

% Make up a fake data set for testing... 

img = phantom(256,'Modified Shepp-logan'); 
kspTrue = fftshift(ifft2(ifftshift( img ))); 
kspTrue = kspTrue + randn(size(kspTrue)).* 0.0001; % add some noise - adding it here will affect the ghost some too!
imgTrue = fftshift(fft2(ifftshift( kspTrue ))); 

% In this simulated data I'll choose reasonable values for kappa and phi
kappaTrue = 0.623; % this is the 1st order parameter (i.e. slope of phase 
             % difference from navigator). It corresponts to a 0.6 shift
             % between RO+ and RO- in k-space
phiTrue = 0.0826; % This is the 0th order parameter (i.e. the intercept of phase
            % difference from navigator). 
            
kspGhost = applyFirstOrderPhaseCorr(kspTrue, kappaTrue, phiTrue); 
imgGhost = fftshift(fft2(ifftshift( kspGhost ))); 

figure(1); 
subplot(1,2,1); 
imagesc(abs(imgTrue)); axis equal; axis tight;  
title('Ideal')

subplot(1,2,2); 
imagesc(abs(imgGhost)); axis equal; axis tight;  
title('Original ghost')

method = 'Gh/Ob'; % options: 'ent','svd','Gh/Ob'


%% Option 1: Discrete search 

% Set up a multiresolution grid for searching -- may need to adjust
% depending on data!!

% Make the ranges larger if you don't know where to start, then take a look
% at the cost function search space to narrow down to speed up the search
startOptDiscrete.kRange = 1; 
startOptDiscrete.pRange = 0.3;
startOptDiscrete.kCenter = 0; % May need to check the polarity here!!!
startOptDiscrete.pCenter = 0;
startOptDiscrete.nGrid = [11, 11, 11]; % This means use three grids of 11x11 each
% ^ I recommend using one larger grid at first so you can see the whole
% metric space for debugging

fprintf('Starting DISCRETE SEARCH \r');

kspOrigTmp = squeeze(kspGhost);
[kappaEstDiscrete, phiEstDiscrete, outputMap] = discreteFun(kspOrigTmp, startOptDiscrete, method);
discreteMap = outputMap;

% We can look at the cost function to see if it seems to be getting closer
% to the minimum 

figure(2); 
subplot(1,3,1); 
imagesc(outputMap.map{1})
subplot(1,3,2); 
imagesc(outputMap.map{2})
title('Cost function');
subplot(1,3,3); 
imagesc(outputMap.map{3})

%% Option 2: Simplex search
% I recommend using a low-res discrete search first, then use that as your
% starting point for simplex search. This will refine your estimate

% Define the starting points here:
% % Option 1 - educated guess:
% startOptSimplex.kCenter = -1; % you need to get the polarity right on this
% startOptSimplex.pCenter = 0; % In my data phi ~= 0 

% % Option 2 - use the solutions from discrete search:
startOptSimplex.kCenter = kappaEstDiscrete; % you need to get the polarity right on this
startOptSimplex.pCenter = phiEstDiscrete; % In my data phi ~= 0 
            
[kappaEstSimplex, phiEstSimplex] = simplexFun(kspGhost, startOptSimplex, method);
                        
 

%% Now apply them and see results

% Apply the first guess from the discrete search
kspCorrRefLessDiscrete = applyFirstOrderPhaseCorr( kspGhost, kappaEstDiscrete, phiEstDiscrete);   
imgRefLessDiscrete = fftshift(fft2(ifftshift( kspCorrRefLessDiscrete ))); 

% Apply the refined guess via simplex search
kspCorrRefLessSimplex = applyFirstOrderPhaseCorr( kspGhost, kappaEstSimplex, phiEstSimplex); 
imgRefLessSimplex = fftshift(fft2(ifftshift( kspCorrRefLessSimplex ))); 

figure(3)
subplot(1,3,1); 
imagesc(abs( imgGhost )); axis equal; axis tight;
title({['True kappa=',num2str(kappaTrue)],['True phi=', num2str(phiTrue)]})
subplot(1,3,2);
imagesc(abs( imgRefLessDiscrete )); axis equal; axis tight;
title({['Discrete kappa=',num2str(kappaEstDiscrete)],['Discrete phi=', num2str(phiEstDiscrete)]})
subplot(1,3,3); 
imagesc(abs( imgRefLessSimplex )); axis equal; axis tight;
title({['Simplex kappa=',num2str(kappaEstSimplex)],['Simplex phi=', num2str(phiEstSimplex)]})



