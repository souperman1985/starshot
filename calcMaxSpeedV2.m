% Calculate maximum speed attainable
%
% Copyright (c) 2025 Matthew Campbell
% Permission is hereby granted, free of charge, to use, copy, modify, and distribute this software for any purpose, with or without modification, provided that the above copyright notice and this permission notice appear in all copies.
% THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
%
% Start fresh
clear all
close all
format short
format compact
clc

% Constants
cLight = 299792458; % m/s ... light speed
hPlanck = 6.62607015e-34; % J-s ... Planck's constant
hPlanckBar = 1.054571817e-34; % J-s ... Planck's constant with a bar... h/(2*pi)
kB = 1.380649e-23; % J/K ... Boltzmann constant
sigmaBB = (2*pi^5*kB^4)/(15*cLight^2*hPlanck^3); % W/m2-K4

% Inputs
maxPwr0 = 100e9; % W ... maximum power producible by light engine
dLaser0 = 30e3; % m ... laser array diameter on Earth 
beta0 = 1e-10; % start from rest... Beta = v/cLight so 1e-10 is ~ 0.03 m/s

% Material parameters
rhoAl2O3 = 3200; % kg/m3
rhoMoS2 = 5060; % kg/m3
rhoSi3N4 = 3200; % kg/m3
rhoSiO2 = 2030; % kg/m3
rhoSi = 2330; % kg/m3
rhoTiO2 = 3700; % kg/m3
sigYldAl2O3 = 2e9; % Pa ... was 5e9.  Jen 2011: interpolate to find yield tensile strain is about 1.3% for 50 nm films.  Multiply by E=170GPa. take only 1 significant figure. 
sigYldMoS2poly = 0.75e9; % Pa ... polycrystalline
sigYldMoS2cryst = 2.3e9; % Pa ... crystalline monolayer
sigYldSi3N4 = 14e9; % Pa ... was 10e9
sigYldSiO2 = 1.5e9; % Pa ... was 2e9
sigYldSi = 2e9; % Pa ... was 7e9.  new value taken from TSUCHIYA 2005. 
sigYldTiO2 = 1.1e9; % Pa ... estimate based on yield strength of 151 GPa (Borgese 2012) and tensile yield strain estimate of 0.75% (Tavares2008) 

% Calculate parameters for each sail
dataSelect = 2;
constraintSelect = 1; % for the paper, we're always using #1. 

% data:
% 1: our prototype
% 2: optimal
% 3: Brewer
% 4: Chang
% 5: Ilic
% 6: Lien
% 7: Salary
% 8: Santi
% 9: Taghavi

% This study
% Uses the reflectivity actually measured by Caltech (this is The January 25, 2022 sample).
if dataSelect == 1
    lambda0 = 1.2; % um ... laser wavelength at zero sail velocity.  
    fNameRef = 'octoberAvgPlusMinusWithWL.csv'; % First column is lambda (nm), second is reflectivity (0-1)
    refData = csvread(fNameRef,1,0); % read the CSV file into a matrix
    refData(:,1) = refData(:,1)*(1/1000); % nm --> um
    hexD = 77e-6; % m ... optical microscope view. Technically the diameter of the hexagons on the Si mold
    hexW = 15e-6; % m ... optical microscope view. Technically the width of the ribs on the Si mold
    hexH = 10e-6; % m ... 95 DRIE smooth loops
    tA1 = 21e-9; % m ... ellipsometry, bottom alumina thickness
    tM = 53e-9; % m ... fit to optical data, MoS2 thickness
    tA2 = 51e-9; % m ... ellipsometry, top alumina thickness
    rhoA1 = rhoAl2O3;
    rhoM = rhoMoS2;
    rhoA2 = rhoAl2O3; 
    rhoAnc = (rhoA1*tA1 + rhoM*tM + rhoA2*tA2)*(1000/1); % g/m2 ... Non-corrugated areal density
    % ribs
    % Muc = ...
    %     (rhoA1 * (tA1*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*(hexD^2-(hexD-2*tA1)^2)) + ...
    %     rhoM * (tM*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD-2*tA1)^2-(hexD-2*tA1-2*tM)^2)) + ...
    %     rhoA2 * (tA2*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD-2*tA1-2*tM)^2-(hexD-2*tA1-2*tM-2*tA2)^2))); % kg ... unit cell mass
    % trenches
    Muc = ...
        (rhoA1 * (tA1*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1).^2-hexD.^2)) + ...
        rhoM * (tM*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1+2*tM).^2-(hexD+2*tA1).^2)) + ...
        rhoA2 * (tA2*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1+2*tM+2*tA2).^2-(hexD+2*tA1+2*tM).^2))); % kg ... unit cell mass  
    Auc = sqrt(3)/2 * (hexD+hexW)^2; % unit cell area
    rhoAc = Muc/Auc *(1000/1); % g/m2 ... Corrugated areal density
    rhoA = rhoAc % g/m2
    G = 1;
    mechRob = G * (tA1*sigYldAl2O3 + tM*sigYldMoS2poly + tA2*sigYldAl2O3) % N/m 
end

% This study - optimized
% Uses the reflectivity actually measured by Caltech (this is The January 25, 2022 sample).
if dataSelect == 2
    lambda0 = 1.2; % um ... laser wavelength at zero sail velocity.  
    fNameRef = 'optimizedSpectrumFull_literatureNK.csv'; % First column is lambda (nm), second is reflectivity (0-1)
    refData = csvread(fNameRef,1,0); % read the CSV file into a matrix
    refData(:,1) = refData(:,1)*(1/1000); % nm --> um
    hexD = 70e-6; % m
    hexW = 2e-6; % m 
    hexH = 3e-6; % m 
    %tA1 = 8.9465e-9; % m 
    %tM = 74.6300e-9; % m 
    %tA2 = 8.9465e-9; % m 
    tA1 = 18.5485e-9; % m
    tM  = 62.6900e-9; % m
    tA2 = 18.5485e-9; % m
    rhoA1 = rhoAl2O3;
    rhoM = rhoMoS2;
    rhoA2 = rhoAl2O3; 
    rhoAnc = (rhoA1*tA1 + rhoM*tM + rhoA2*tA2)*(1000/1); % g/m2 ... Non-corrugated areal density
    % ribs
    % Muc = ...
    %     (rhoA1 * (tA1*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*(hexD^2-(hexD-2*tA1)^2)) + ...
    %     rhoM * (tM*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD-2*tA1)^2-(hexD-2*tA1-2*tM)^2)) + ...
    %     rhoA2 * (tA2*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD-2*tA1-2*tM)^2-(hexD-2*tA1-2*tM-2*tA2)^2))); % kg ... unit cell mass
    % trenches
    Muc = ...
        (rhoA1 * (tA1*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1).^2-hexD.^2)) + ...
        rhoM * (tM*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1+2*tM).^2-(hexD+2*tA1).^2)) + ...
        rhoA2 * (tA2*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1+2*tM+2*tA2).^2-(hexD+2*tA1+2*tM).^2))); % kg ... unit cell mass  
    Auc = sqrt(3)/2 * (hexD+hexW)^2; % unit cell area
    rhoAc = Muc/Auc *(1000/1); % g/m2 ... Corrugated areal density
    rhoA = rhoAc % g/m2
    G = 1;
    mechRob = G * (tA1*sigYldAl2O3 + tM*sigYldMoS2cryst + tA2*sigYldAl2O3) % N/m 
end

% Brewer 2022
if dataSelect == 3
    lambda0 = 1.2; % um ... laser wavelength at zero sail velocity.  
    fNameRef = 'Brewer2022reflectivity.csv'; % First column is lambda (nm), second is reflectivity (0-1)
    refData = csvread(fNameRef,1,0); % read the CSV file into a matrix
    refData(:,1) = refData(:,1)*(1/1000); % nm --> um
    periodY = 1.16e-6; % m
    periodX = periodY * sqrt(3); % m
    dHole = 0.9*periodY; % m
    tSN1 = 5e-9; % m ... Si3N4 thickness on top,  see Fig 2
    tM = 90e-9; % m ... MoS2 thickness
    tSN2 = 5e-9; % m ... Si3N4 thickness on bottom
    Muc = (periodY * periodX - 2*pi*(dHole/2)^2) * (rhoSi3N4 * (tSN1+tSN2) + rhoMoS2 * tM); % kg ... each unit cell contains 2 holes
    Auc = periodY * periodX; % m2
    rhoA = Muc/Auc*(1000/1) % g/m2
    G = (periodY-dHole)/periodY; 
    mechRob = G * (tSN1*sigYldSi3N4 + tM*sigYldMoS2poly + tSN2*sigYldSi3N4) % N/m
end

% Chang 2024... design from Fig 3
if dataSelect == 4
    lambda0 = 1.3; % um ... laser wavelength at zero sail velocity. 
    fNameRef = 'Chang2023exptRcalib.csv'; % First column is lambda (nm), second is reflectivity (0-1)
    refData = csvread(fNameRef,1,0); % read the CSV file into a matrix
    refData(:,1) = refData(:,1)*(1/1000); % nm --> um
    tSN = 400e-9; % m
    tSi = 321e-9; % m
    aPattern = 1440e-9; % m ... measured from SEM in Figure 3
    dHole = 2*500e-9; % m ... confirmed from Fig 3 SEM
    Muc = (rhoSi3N4 * tSN) * (aPattern^2 - pi*(dHole/2)^2) + rhoSi*tSi*aPattern^2; % kg
    Auc = aPattern^2; % m2
    rhoA = Muc/Auc*(1000/1) % g/m2
    G = (aPattern-dHole)/aPattern;
    mechRob = G * (tSN*sigYldSi3N4) + 1*tSi*sigYldSi % N/m ... G for the Si part = 1 since it's not patterned. 
end

% Ilic 2018... design A11 (see sup inf)
if dataSelect == 5
    lambda0 = 1.2; % um ... laser wavelength at zero sail velocity. 
    fNameRef = 'Ilic2018reflectivityDesignA11_usingRodriguez_de_Marcos.csv'; % First column is lambda (nm), second is reflectivity (0-1)
    refData = csvread(fNameRef,1,0); % read the CSV file into a matrix
    refData(:,1) = refData(:,1)*(1/1000); % nm --> um
    tSiO2total = sum([0.104, 0.124, 0.131, 0.131, 0.124, 0.104])*(1/1e6); % um --> m
    rhoA = rhoSiO2*tSiO2total*(1000/1) % kg/m2 ==> g/m2
    G = 1;
    mechRob = G * (tSiO2total*sigYldSiO2) % N/m
end

% Lien 2022
if dataSelect == 6
    lambda0 = 1.064; % um ... laser wavelength at zero sail velocity.  
    fNameRef = 'Lien2022reflectivity.csv'; % First column is lambda (nm), second is reflectivity (0-1)
    refData = csvread(fNameRef,1,0); % read the CSV file into a matrix
    refData(:,1) = refData(:,1)*(1/1000); % nm --> um
    tSN = 690e-9; % m
    aPattern = 1064e-9; % m
    dHole = 415e-9; % m
    Muc = (rhoSi3N4 * tSN) * (aPattern^2 - pi*(dHole/2)^2); % kg
    Auc = aPattern^2; % m2
    rhoA = Muc/Auc*(1000/1); % g/m2
    G = (aPattern-dHole)/aPattern;
    mechRob = G * (tSN*sigYldSi3N4) % N/m
end

% Salary 2019
if dataSelect == 7
    lambda0 = 1.3; % um ... laser wavelength at zero sail velocity.  
    fNameRef = 'Salary2019reflectivity.csv'; % First column is lambda (nm), second is reflectivity (0-1)
    refData = csvread(fNameRef,1,0); % read the CSV file into a matrix
    refData(:,1) = refData(:,1)*(1/1000); % nm --> um
    tSiO2 = 60e-9; % m
    tSi = 89e-9; % m
    %delta = 639e-9; % m
    %hDisk = 323e-9; % m
    dCircle = 300e-9; % m ... range is 180 to 450 nm, and I pulled optical properties at 300 nm.
    rhoA180 = 0.413; % g/m2 ... value given for 180 nm circles
    rhoA450 = 0.660; % g/m2 ... value given for 450 nm circles
    rhoA = (rhoA180+rhoA450)/2 % g/m2... just take an average
    G = 1;
    mechRob = G * (tSiO2*sigYldSiO2 + tSi*sigYldSi) % N/m
end

% Santi 2022 ... design with TiO2 + SiO2 + TiO2
if dataSelect == 8
    lambda0 = 1.064; % um ... laser wavelength at zero sail velocity. 
    fNameRef = 'Santi2022reflectivityTiO2SiO2TiO2_usingRodriguez_de_Marcos.csv'; % First column is lambda (nm), second is reflectivity (0-1)
    refData = csvread(fNameRef,1,0); % read the CSV file into a matrix
    refData(:,1) = refData(:,1)*(1/1000); % nm --> um
    tTiO2 = 121e-9; % m
    tSiO2 = 203e-9; % m
    rhoA = (rhoTiO2 * 2*tTiO2 + rhoSiO2 * tSiO2)*(1000/1); % g/m2
    G = 1;
    mechRob = G * (2*tTiO2*sigYldTiO2 + tSiO2*sigYldSiO2); % N/m ... G for the Si part = 1 since it's not patterned. 
end

% Taghavi 2022
if dataSelect == 9
    lambda0 = 1.3; % um ... laser wavelength at zero sail velocity. 
    fNameRef = 'Taghavi2022reflectivity.csv'; % First column is lambda (nm), second is reflectivity (0-1)
    refData = csvread(fNameRef,1,0); % read the CSV file into a matrix
    refData(:,1) = refData(:,1)*(1/1000); % nm --> um
    tSiO2 = 60e-9; % m
    tSi = 103e-9; % m
    dCircle = 300e-9; % m ... range is 180 to 450 nm, and I pulled optical properties at 300 nm.
    rhoA = 4.7 / (2*2) % g/m2... 2x2 m2 sail said to have mass of 4.7 g
    G = 1;
    mechRob = G * (tSiO2*sigYldSiO2 + tSi*sigYldSi) % N/m
end

% Sail diameter... assumes spherically curved sail with s=d (radius of curvature equal to diameter)
if constraintSelect == 1 % specify total mass
    constraintString = "mass";
    Mtot = 2; % g ... total mass of sail, tethers, and chip
    Msail = Mtot/2; % g ... mass of sail only
    Asail = Msail/rhoA; % m2
    sSail = sqrt(Asail/(pi*(2-sqrt(3)))); % m ... equation for radius of sail when s=d ... see Equation S44 in Campbell 2022 sup inf
    dSail = sSail; % m ... sail diameter 
    AsailPerp = pi*(dSail/2)^2; % m
end

if constraintSelect == 2 % specify sail diameter 
    constraintString = "diameter";
    dSail = 2; % m
    sSail = dSail; % guideline from Campbell 2022
    AsailPerp = pi*(dSail/2)^2; % m
    Asail = pi*(2-sqrt(3))*sSail^2; % m ... equation for spherically curved sail area when s=d ... see Equation S44 in Campbell 2022 sup inf
    Msail = Asail*rhoA; % g
    Mtot = 2*Msail;
end

% Light wavelengths at the sail
lambdaMax = refData(end,1); % um
digitsOld = digits(10); % reduce precision temporarily
syms betaFind; % declare variable
betaF = vpasolve(lambdaMax == lambda0*(1+betaFind)/sqrt(1-betaFind^2), betaFind, 0.2); % upper limit to the fraction of speed of light that we want to achieve... if sails could get going this fast. 
digits(digitsOld); % restore precision
betaF = double(betaF); % increase to double precision
if dataSelect == 4 % for Chang2023
    disp('Allowing extrapolation of reflectivity data!')
    betaF = 0.2; % force extrapolation... the original data aren't even enough to get to beta=0.2
end
betaVals = linspace(beta0, betaF, 10000); % fraction of speed of light
gammaVals = 1./sqrt(1-betaVals.^2); % relativistic Lorentz factor
lambdaInVals = lambda0 * sqrt((1+betaVals)./(1-betaVals)); % um ... wavelengths of light that hit the sail

% Iterpolate within reflectivity values at beta steps
interpFxnRef = griddedInterpolant(refData(:,1),refData(:,2),'linear','linear'); % create function to describe data and enable interpolations.  
RefSpectrum = interpFxnRef(lambdaInVals); % arbs ... find reflectivity at the wavelengths corresponding to the beta values we defined above

% Estmate maximum acceleration distance possible with focus on sail, given sail diameter and laser array diameter and wavelength
Lmax = dSail*dLaser0/(2*lambda0*(1/1e6)) % m

% Iteration parameters... want to find the power that takes the sail within "failureMargin" of breaking. 
failureMargin = 0.001; % 0.1-percent failure margin
maxPwr0nr = 1e15; % W ... maximum power for Newton Raphson algorithm... beyond this and it's probably an error.
minPwr0nr = 1e6; % W ... min power... lower than this and it's probably an error. 
maxIterPwr0 = 100; % max iterations
maxTimePwr0 = 3; % seconds 
tolerancePwr0 = 0.0005; % convergence tolerance
deltaPctPwr0 = 0.005; % percent difference for computing derivative

% Initial guesses
Pwr0guessA = 1e9; % W = kg-m2/s3 ... laser power
Pwr0guessB = Pwr0guessA + deltaPctPwr0*Pwr0guessA; % point B guess

% Iterative loop 
stopper=0;
count = 1;
tStartNR = tic;
while stopper==0
    % for point A
    [mechRatioMaxA, distValsA, timeValsA, accValsLaserA, accValsSailA, PvalsA, breakingForceValsA] = solveSail(...
            Pwr0guessA, cLight, Mtot, sSail, AsailPerp, Lmax, mechRob, betaVals, gammaVals, RefSpectrum);
    fRA = mechRatioMaxA;
    fA = fRA-(1-failureMargin);
    if abs(fA/(1-failureMargin))<=tolerancePwr0; break; end % exit if we're within tolerance
    % for point B.  
    [mechRatioMaxB, ~, ~, ~, ~, ~, ~] = solveSail(...
            Pwr0guessB, cLight, Mtot, sSail, AsailPerp, Lmax, mechRob, betaVals, gammaVals, RefSpectrum);
    fRB = mechRatioMaxB;
    fB = fRB-(1-failureMargin);
    % next guesses
    Pwr0guessAold = Pwr0guessA; % store previous A guess value
    Pwr0guessA = Pwr0guessA + (deltaPctPwr0*Pwr0guessA)*((fA)/(fA-fB));
    Pwr0guessB = Pwr0guessA + deltaPctPwr0*Pwr0guessA;
    % now check our next input... don't run away! 
    if Pwr0guessA > maxPwr0nr
        Pwr0guessA = (Pwr0guessAold+maxPwr0nr)/2; % take the average 
        Pwr0guessB = Pwr0guessA + deltaPctPwr0*Pwr0guessA;
    elseif Pwr0guessA < minPwr0nr
        Pwr0guessA = (Pwr0guessAold+minPwr0nr)/2;
        Pwr0guessB = Pwr0guessA + deltaPctPwr0*Pwr0guessA;
    end
    % prevent a runaway... too many iterations or too much time
    count = count+1;
    if (count > maxIterPwr0); fprintf('RUNAWAY: count > maxIter \n'); Pwr0guessA=NaN; break; end
    tElapsedNR = toc(tStartNR);
    if (tElapsedNR > maxTimePwr0); fprintf('RUNAWAY: tElapsed > maxTime \n'); Pwr0guessA=NaN; break; end
end
tElapsedNR = toc(tStartNR); % timer off
Pwr0found = Pwr0guessA; % W
actualPower = Pwr0found;

% Check if we went over the allowed power
usingMaxPower = 0;
if Pwr0found>maxPwr0 % if a greater power was allowable, than surely a lower power would be as well.  Just use the max allowable. 
    [mechRatioMaxA, distValsA, timeValsA, accValsLaserA, accValsSailA, PvalsA, breakingForceValsA] = solveSail(...
            maxPwr0, cLight, Mtot, sSail, AsailPerp, Lmax, mechRob, betaVals, gammaVals, RefSpectrum);
    usingMaxPower = 1;
    actualPower = maxPwr0;
else % if we used the max power or less than the max
    % check which Beta limited the power.  We want to make sure this beta is
    % less than or equal to betamax, otherwise the limit occurs outside of the
    % actual practical range and we could have had higher power. 
    betaWhereLimitHappened = betaVals(mechRatioMaxA==breakingForceValsA/mechRob)
end

powerAppliedGW = actualPower/1e9

tempVec = find(isnan(distValsA)); % vector of indicees corresponding to NaN, which are values where distance is too long 
indexLast = tempVec(1)-1; % last index for which the acceleration distance was less than the diffraction limit of distance
betaMax = betaVals(indexLast)
Lfinal = distValsA(indexLast)

% Final calculations
% Average reflectivity over beta=0.00-0.20 range... see Equation S100 in Campbell 2022 Sup Inf.  This is not just an average over wavelengths; it needs to be weighted by beta or by lambda (for fine enough gridding, these should be equal). 
% Turns out that averaging over beta gives a slightly different result from averaging over wavelength.
% In theory, doing this with enough points should make them converge, but I haven't seen that happen for some reason. 
% For this paper, I'm going to ****average over wavelength for everything***. 
% So I integrate Rdlambda then divide by integral dlambda = lambdaF-lambdaI.
betaFstarshot = 0.2;
[dummy,idxBetaFstarshot] = min(abs(betaVals-betaFstarshot));
betaFstarshotHave = betaVals(idxBetaFstarshot)
%meanRefStarshot = mean(RefSpectrum(1:idxBetaFstarshot)) % note that for linearly-spaced beta values this is equivalent to Equation S100: int(Ri*dBeta)/(betaF-betaI) = sum(Ri*dBeta)/(betaF-betaI) but dBeta = (betaF-betaI)/N --> sum(Ri*dBeta)/(betaF-betaI) = sum(Ri)*dBeta/(betaF-betaI) = sum(Ri)/N = average calculated by Matlab
meanRefStarshot = sum((1/2)*(RefSpectrum(1:idxBetaFstarshot-1)+RefSpectrum(2:idxBetaFstarshot)).*(lambdaInVals(2:idxBetaFstarshot)-lambdaInVals(1:idxBetaFstarshot-1)))/(lambdaInVals(idxBetaFstarshot)-lambdaInVals(1))

% check reflectivity
lambdaMaxStarshot = lambda0 * sqrt((1+betaFstarshot)./(1-betaFstarshot));
lambdaValsCheck = linspace(lambda0, lambdaMaxStarshot, 10000);
RefSpectrumCheck = interpFxnRef(lambdaValsCheck); % arbs ... find reflectivity at the wavelengths corresponding to the beta values we defined above
meanRefStarshotCheck = mean(RefSpectrumCheck)
% Looks like reflectivity values from beta and lambda don't quite match... but doing the integration this way and with a weighted average over the non-linearly-spaced lambda points is equivalent.  ***Should use lambda in this case then.***

% Output string
clipStr = sprintf('%10s, %10s, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e',...
    fNameRef,constraintString,maxPwr0*(1/1e9),dLaser0*(1/1000),lambda0,rhoA,mechRob,Mtot,Msail,Asail,AsailPerp,sSail,dSail,beta0,betaF,Lmax*(1/1e9),failureMargin,Pwr0found*(1/1e9),usingMaxPower,actualPower*(1/1e9),betaMax,Lfinal*(1/1e9),meanRefStarshot);
clipStr
clipboard('copy',clipStr);

% Plots
figure(1)
hold on;
plot(lambdaInVals,RefSpectrum,'b-')
plot(refData(:,1),refData(:,2),'k-')
xlabel('Wavelength [um]')
ylabel('Reflectivity [a.u.]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 1.5);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);

figure(2)
plot(betaVals,distValsA*(1/1e9))
xlabel('Beta')
ylabel('Distance [Gm]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 1.5);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);

figure(3)
plot(betaVals,PvalsA)
xlabel('Beta')
ylabel('Pressure [Pa]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 1.5);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);

figure(4)
hold on
plot(betaVals,breakingForceValsA)
plot([betaVals(1),betaVals(end)],[mechRob,mechRob],'k--')
hold off
xlabel('Beta')
ylabel('Breaking force [N/m]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 1.5);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);

% accelDistGm = accelDist*(1/1e9) % Gm
% dLaser0km = dLaser0*(1/1e3) % km
% accelTimeMinutes = accelTime*(1/60) % minutes
% laserOnTimeMinutes = laserOnTime*(1/60) % minutes
% 

%% Functions

function [mechRatioMax, distVals, timeVals, accValsLaser, accValsSail, Pvals, breakingForceVals] = solveSail(...
            Pwr0, cLight, Mtot, sSail, AsailPerp, Lmax, mechRob, betaVals, gammaVals, RefSpectrum)
    % This solves for the sail's acceleration parameters given a laser power value [W]
    % Integrate to find the acceleration distance, etc.
    % distance
    integrandVec1 = (1./RefSpectrum).*((gammaVals.*betaVals)./(1-betaVals).^2); % ## (no units)
    integralVec1 = zeros(size(integrandVec1)); % ## ... the first values of the integral are zero
    integralVec1(2:end) = cumsum((1/2)*(integrandVec1(1:end-1)+integrandVec1(2:end)).*diff(betaVals)); % ## ... this is doing cumulative trapezoidal integration 
    distVals = (cLight^3/(2*Pwr0))*(Mtot*(1/1000))*integralVec1; % m ... distance values throughout acceleration
    distVals(distVals>Lmax) = NaN; % remove distances outside bounds of focus
    % time
    integrandVec2 = (gammaVals.^3./RefSpectrum).*((1+betaVals)./(1-betaVals)); % ## (no units)
    integralVec2 = zeros(size(integrandVec2)); % ## ... the first values of the integral are zero
    integralVec2(2:end) = cumsum((1/2)*(integrandVec2(1:end-1)+integrandVec2(2:end)).*diff(betaVals)); % ## ... this is doing cumulative trapezoidal integration
    timeVals = (cLight^2/(2*Pwr0))*(Mtot*(1/1000))*integralVec2; % s ... time values throughout acceleration at beta values, in laser's reference frame
    % acceleration and pressure
    accValsLaser = ((2*Pwr0*RefSpectrum)./(Mtot*(1/1000)*cLight*gammaVals.^3)).*((1-betaVals)./(1+betaVals)); % m/s2 laser ref frame
    accValsSail = accValsLaser.*gammaVals.^3; % m/s2 sail ref frame
    Pvals = ((Mtot*(1/1000)).*accValsSail)/AsailPerp; % (kg*m/s2)/m2 = Pa ... photon pressure on sail 
    % mechanical robustness criterion
    breakingForceVals = Pvals*sSail/2; % N/m 
    mechRatioMax = max(breakingForceVals/mechRob);
end




