% Starshot composite corrugated spherically curved sail
% Matthew Campbell
% 2024-05-25
%
% Copyright (c) 2025 Matthew Campbell
% Permission is hereby granted, free of charge, to use, copy, modify, and distribute this software for any purpose, with or without modification, provided that the above copyright notice and this permission notice appear in all copies.
% THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
%
% Solve for constnat power that laser can tolerate, up to 100 GW. Accel to beta = 0.2

% Start fresh
clear all
close all
format short
format compact
clc

writeFiles = 0;

% Constants
cLight = 299792458; % m/s ... light speed
hPlanck = 6.62607015e-34; % J-s ... Planck's constant
kB = 1.380649e-23; % J/K ... Boltzmann constant
sigmaBB = (2*pi^5*kB^4)/(15*cLight^2*hPlanck^3); % W/m2-K4

% Inputs
Pwr0limit = 100e9; % W ... maximum power producible by light engine
%dLaser0 = 30e3; % m ... laser array diameter on Earth 
beta0 = 1e-10; % start from rest... Beta = v/cLight so 1e-10 is ~ 0.03 m/s
lambda0 = 1.2; % um ... laser wavelength at zero sail velocity.  
%theta = 0; % normal incidence
Toptical = 300; % K ... for evaluating optical properties.  It's the only T at which we have sufficient data. 
failureMargin = 0.0001; % 0.01-percent failure margin
kappaM0 = 1e-8 % Specify this for kappa of MoS2 at 1.2 um. 

% Material parameters
rhoAl2O3 = 3200; % kg/m3
rhoMoS2 = 5060; % kg/m3
sigYldAl2O3 = 2e9; % Pa ... was 5e9
sigYldMoS2 = 2.3e9; % Pa ... was 0.75e9
%TmeltMoS2 = 1050+273.15; % K ... actually this is more like a sublimination temperature
%TmeltMoS2 = 1458; % K
TmeltMoS2 = 1000; % K ... Cui 2018, this is an estimate at the sublimination temperature under vacuum (in space!)
Tmelt = TmeltMoS2; % K

% thickness values
tA = 18.5485*(1/1e9); % nm --> m ... Al2O3
tM = 62.6900*(1/1e9); % nm --> m ... MoS2
% hexagonal parameters
hexD = 70e-6; % m ... hexagonal corrugation parameters
hexW = 2e-6; % m 
hexH = 3e-6; % m 

% areal density and mechanical robustness
% ribs
%Muc = ...
%    rhoAl2O3 * (tA*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*(hexD^2-(hexD-2*tA).^2)) + ...
%    rhoMoS2 * (tM*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD-2*tA).^2-(hexD-2*tA-2*tM).^2)) + ...
%    rhoAl2O3 * (tA*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD-2*tA-2*tM).^2-(hexD-2*tA-2*tM-2*tA).^2)); % kg ... unit cell mass
% trenches
Muc = ...
    (rhoAl2O3 * (tA*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA).^2-hexD.^2)) + ...
    rhoMoS2 * (tM*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA+2*tM).^2-(hexD+2*tA).^2)) + ...
    rhoAl2O3 * (tA*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA+2*tM+2*tA).^2-(hexD+2*tA+2*tM).^2))); % kg ... unit cell mass  
Auc = sqrt(3)/2 * (hexD+hexW)^2; % unit cell area
rhoAc = Muc/Auc *(1000/1) % g/m2 ... Corrugated areal density
G = 1;
mechRob = G * (tA*sigYldAl2O3 + tM*sigYldMoS2 + tA*sigYldAl2O3) % N/m 

% Sail diameter... assumes spherically curved sail with s=d (radius of curvature equal to diameter)
% specify the total mass and derive the diameter. 
Mtot = 2; % g ... total mass of sail, tethers, and chip
Msail = Mtot/2; % g ... mass of sail only defined to be half of total
McMt = Msail; % g ... chip and tether mass combined
Asail = Msail./rhoAc; % m2
sSail = sqrt(Asail/(pi*(2-sqrt(3)))); % m ... equation for radius of sail when s=d ... see Equation S44 in Campbell 2022 sup inf
dSailP = sSail; % m ... sail diameter perpendicular to incoming laser photons
AsailP = pi*(dSailP/2).^2; % m

% Light wavelengths at the sail
betaF = 0.2; % Starshot
betaVals = linspace(beta0, betaF, 1001)'; % fraction of speed of light ... COLUMN vector!
gammaVals = 1./sqrt(1-betaVals.^2); % relativistic Lorentz factor ... COLUMN vector!
lambdaInVals = lambda0 * sqrt((1+betaVals)./(1-betaVals)); % um ... wavelengths of light that hit the sail ... COLUMN vector!

% thickness values for optical calculations
hVec = [NaN tA*(1e6/1) tM*(1e6/1) tA*(1e6/1) NaN]; % um... film layer thickness

% Read Lingart 1982 data for attenuation coefficient (alpha) for alumina
fNameL1982a = 'Lingart1981_T2_data.csv'; % data from table 2 in that publication.  rows show lambda in um, columns are T in K
tempMatrixL1982a = csvread(fNameL1982a,1,0); % read the CSV file into a matrix
lambdaValsL1982a = tempMatrixL1982a(2:end,1); % um ... first column is wavelength
TvalsL1982a = tempMatrixL1982a(1,2:end); % K ... first row is temperature 
spectraL1982a = tempMatrixL1982a(2:end,2:end); % cm-1, (alpha) Naperien attenuation coefficient data
spectraL1982k = (1/(4*pi))*spectraL1982a.*repmat(lambdaValsL1982a,1,size(spectraL1982a,2))*(100/1)*(1/1e6); % convert to extinction coefficient ... k = alpha*lambda/(4*pi) 
[lambdaMatL1982a, TmatL1982a] = ndgrid(lambdaValsL1982a, TvalsL1982a); % matricees saying the wavelength and temperature at all values in spectraAD
interpFxnL1982k = griddedInterpolant(lambdaMatL1982a,TmatL1982a,spectraL1982k,'linear','linear'); % create function to describe data and enable interpolations.  

% Read Lingart 1982 data for index of refraction (n) for alumina
fNameL1982n = 'Lingart1981_T4_data.csv'; % data from table 4 in that publication.  rows show lambda in um, columns are T in K
tempMatrixL1982n = csvread(fNameL1982n,1,0); % read the CSV file into a matrix
lambdaValsL1982n = tempMatrixL1982n(2:end,1); % um ... first column is wavelength
TvalsL1982n = tempMatrixL1982n(1,2:end); % K ... first row is temperature 
spectraL1982n = tempMatrixL1982n(2:end,2:end); % arbitrary units, (n) index of refraction data
[lambdaMatL1982n, TmatL1982n] = ndgrid(lambdaValsL1982n, TvalsL1982n); % matricees saying the wavelength and temperature at all values in spectraAD
interpFxnL1982n = griddedInterpolant(lambdaMatL1982n,TmatL1982n,spectraL1982n,'linear','linear'); % create function to describe data and enable interpolations.  

% Read Querry 1985 data for n and k at 300 K for alumina
fNameQ1985 = 'Querry1985-eo-complexRefIndData.csv'; % First column is lambda (um), second is n (arbs), third is k (arbs)
tempMatrixQ1985 = csvread(fNameQ1985,2,0); % read the CSV file into a matrix
lambdaValsQ1985 = tempMatrixQ1985(:,1); % um ... first column is wavelength
spectrumQ1985ne = tempMatrixQ1985(:,2); % arbitrary units, (n) index of refraction data @ 300 K... extraordinary
spectrumQ1985ke = tempMatrixQ1985(:,3); % arbitrary units, (k) extinction coefficient data @ 300 K
spectrumQ1985no = tempMatrixQ1985(:,4); % arbitrary units, (n) index of refraction data @ 300 K... ordinary
spectrumQ1985ko = tempMatrixQ1985(:,5); % arbitrary units, (k) extinction coefficient data @ 300 K
spectrumQ1985ke(spectrumQ1985ke<0)=0; % set all negative values to zero
spectrumQ1985ko(spectrumQ1985ko<0)=0; 
spectrumQ1985n = (spectrumQ1985ne+spectrumQ1985no)/2; % take the average
spectrumQ1985k = (spectrumQ1985ke+spectrumQ1985ko)/2; 
interpFxnQ1985n = griddedInterpolant(lambdaValsQ1985,spectrumQ1985n,'linear','linear'); % create function to describe data and enable interpolations.  
interpFxnQ1985k = griddedInterpolant(lambdaValsQ1985,spectrumQ1985k,'linear','linear'); % create function to describe data and enable interpolations.  

% Read Munkhbat 2022 MoS2 data
fNameMunk = 'Munkhbat2022MoS2nk.csv'; %fprintf('Using Single Layer MoS2 Data\n');
tempMatrixMunk = csvread(fNameMunk,2,0); % read the CSV file into a matrix
lambdaValsMunk = tempMatrixMunk(:,1); % um ... first column is wavelength
spectrumMunkN = tempMatrixMunk(:,2); % arbitrary units, (n) index of refraction data @ 300 K
spectrumMunkK = tempMatrixMunk(:,3); % arbitrary units, (k) extinction coefficient data @ 300 K

% Read Extension of Munkhbat from jason
fNameMunkExpans = 'MunkhbatExpansionLong.csv'; %fprintf('Using Single Layer MoS2 Data\n');
tempMatrixMunkExpans = csvread(fNameMunkExpans,2,0); % read the CSV file into a matrix
lambdaValsMunkExpans = tempMatrixMunkExpans(:,1); % um ... first column is wavelength
spectrumMunkExpansN = tempMatrixMunkExpans(:,2); % arbitrary units, (n) index of refraction data @ 300 K
spectrumMunkExpansK = tempMatrixMunkExpans(:,3); % arbitrary units, (k) extinction coefficient data @ 300 K

% concatenate extrap data to literature data and then create interpolaters
interpFxnMoS2n = griddedInterpolant([lambdaValsMunk; lambdaValsMunkExpans],[spectrumMunkN; spectrumMunkExpansN],'linear','linear'); % create function to describe data and enable interpolations.  
interpFxnMoS2k = griddedInterpolant([lambdaValsMunk; lambdaValsMunkExpans],[spectrumMunkK; spectrumMunkExpansK],'linear','linear'); % create function to describe data and enable interpolations.  

% Emission at longer wavelengths
% Range to explore.  20 um is about the highest we can reliably do for MoS2 before plasmons come into play and we don't have data on those. 
lambdaOutVals = [lambdaValsMunk(1):0.001:20]'; % um ... wavelengths to calculate emissivity.  Start at minimum of Ermolaev2020, go out to max of Querry. 

% Adjust kappa values of mos2
% Interpolate to find MoS2 Ermolaev 2020 n and k data at lambda and 300 K
nMe = interpFxnMoS2n(lambdaOutVals); % arbs
kMe = interpFxnMoS2k(lambdaOutVals); % arbs
% replace all values greater than laser wavelength with our assumed value
indOutLambda0 = find(lambdaOutVals==lambda0); % find index of 1.2 um
lambdaHave = lambdaOutVals(indOutLambda0) % um
kMe(indOutLambda0:end) = kappaM0; % assume this value for the laser wavelength out to long wavelengths.  This will underestimate emissivity slightly. 

% using a combination of Lingart and Querry: 
% Find Al2O3 Querry 1985 n and k data at lambda and 300 K
nQ1985e = interpFxnQ1985n(lambdaOutVals); % arbs
kQ1985e = interpFxnQ1985k(lambdaOutVals); % arbs
% Interpolate to find Lingart 1982 k data at lambda and 300 K
[lambdaValsCalcMat, TvalsCalcMat] = ndgrid(lambdaOutVals, Toptical); % generate input vectors for the griddedInterpolant interpolation
kL1982e = interpFxnL1982k(lambdaValsCalcMat, TvalsCalcMat); % arbs ...  Perform the interpolation in terms of wavelength and temperature. lambda increases down each column, and temperature increases across each row.
nL1982e = interpFxnL1982n(lambdaValsCalcMat, TvalsCalcMat);
% Fit power law to bridge between Lingart and Querry k-datasets
lambdaA = lambdaValsL1982a(end); % max wavelength given by Lingart
indLambdaA = find(lambdaOutVals==lambdaA);
lambdaHaveA = lambdaOutVals(indLambdaA) % um
kappaA = kL1982e(indLambdaA); % from Lingart
lambdaB = 10; % um
indLambdaB = find(lambdaOutVals==lambdaB);
lambdaHaveB = lambdaOutVals(indLambdaB) % um
kappaB = kQ1985e(indLambdaB) % from Querry
b = log(kappaA/kappaB)/log(lambdaA/lambdaB);
a = kappaA/lambdaA^b;
kFite = a*lambdaOutVals.^b; % evaluate the power law everywhere (we'll only use some of these values though)
% Assemble a unified k matrix from the pieces
kUnife = kL1982e; % start with Lingart 1982 data
kUnife(indLambdaA+1:indLambdaB-1) = kFite(indLambdaA+1:indLambdaB-1); % insert fit
kUnife(indLambdaB:end) = kQ1985e(indLambdaB:end); % insert Querry 1985 data that has been T-scaled.
nUnife = nQ1985e; % column vector ... just use n-data from Querry  

% Reflection and Absorption at laser wavelength (Doppler range)

% Interpolate to find MoS2 Ermolaev 2020 n and k data at lambda and 300 K
nMra = interpFxnMoS2n(lambdaInVals); % arbs
kAssumeMra = kappaM0*ones(size(lambdaInVals)); % we assume it's constant in laser band
indComplexInMra = nMra + 1i*kAssumeMra; % column vector

% Use Querry n-data for alumina
nAluminaRA = interpFxnQ1985n(lambdaInVals); % arbs
[lambdaValsCalcMat, TvalsCalcMat] = ndgrid(lambdaInVals, Toptical); % generate input vectors for the griddedInterpolant interpolation
kL1982ra = interpFxnL1982k(lambdaValsCalcMat, TvalsCalcMat); % arbs ...  Perform the interpolation in terms of wavelength and temperature. lambda increases down each column, and temperature increases across each row.
nL1982ra = interpFxnL1982n(lambdaValsCalcMat, TvalsCalcMat);
kAluminaRA = kL1982ra;
indComplexInAra = nAluminaRA + 1i*kAluminaRA; % column vector

% Assemble complex index of refraction
% BE CAREFUL WHEN ASSEMBLING COMPLEX NUMBERS!  MAKE SURE ALL ARE COLUMN VECTORS!!!
nMatra = [ones(length(lambdaInVals),1) indComplexInAra indComplexInMra indComplexInAra ones(length(lambdaInVals),1)];
% Assemble theta values over which to calculate spectra. They are dependent on the sail diameter and radius of curvature
thetaMaxra = asin((dSailP/2)/sSail); % rad ... maximum angle of reflection based on geometry is a function of sail diameter and curvature
nThetara = 99; % number of theta values to explore 
thetaVecra = linspace(0,thetaMaxra,nThetara); % rad ... angles to evaluate R and A

%plot it
figure(100)
hold on
plot(lambdaOutVals,nQ1985e,'b-')
plot(lambdaOutVals,nL1982e,'g-')
plot(lambdaOutVals,nUnife,'k-')
hold off
set(gca, 'XScale', 'log')
figure(101)
hold on
plot(lambdaOutVals,kQ1985e,'b-')
plot(lambdaOutVals,kL1982e,'g-')
plot(lambdaOutVals,kUnife,'k-')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off
figure(102)
hold on
plot(lambdaOutVals,nMe,'b-')
hold off
set(gca, 'XScale', 'log')
figure(103)
hold on
plot(lambdaOutVals,kMe,'b-')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off
figure(104)
hold on
plot(lambdaOutVals,nMe,'r:')
plot(lambdaOutVals,nUnife,'b--')
hold off
set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
figure(105)
hold on
plot(lambdaOutVals,kMe,'r:')
plot(lambdaOutVals,kUnife,'b--')
hold off
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')




% Assemble complex index of refraction
% BE CAREFUL WHEN ASSEMBLING COMPLEX NUMBERS!  MAKE SURE ALL ARE COLUMN VECTORS!!!
% ALSO: The prime (') operator works differently on complex vectors than on real vectors.  BE CAREFUL!!!
indComplexOutMe = nMe + 1i*kMe; % column vector
indComplexOutAe = nUnife + 1i*kUnife; % column vector
nMate = [ones(length(lambdaOutVals),1) indComplexOutAe indComplexOutMe indComplexOutAe ones(length(lambdaOutVals),1)];
% Theta values for integration
nThetae = 101; % if the maximum angle is 77.5 degrees, there are at least 0.775 degrees per step.
%thetaMaxe = 77.5*(pi/180); % degrees --> rad.  Larger angles than this produce strange results in the transform matrix.
thetaMaxe = 89*(pi/180); % degrees --> rad.  
thetaVece = linspace(0,thetaMaxe,nThetae); % rad ... angles to evaluate E
thetaMateEx = repmat([thetaVece 90*(pi/180)],[length(lambdaOutVals),1]); % tack on a 90-degree angle and pattern over all lambda values

% Gather optical properties
% Calculate R/T/A for laser wavelength range at angles corresponding to current sail geometry
pol = 1; % parallel polarization
[Rmatrap, ~, Amatrap, ~, ~] = FresnelMat(lambdaInVals, thetaVecra, hVec, nMatra, pol); 
pol = 0; % perpendicular polarization
[Rmatras, ~, Amatras, ~, ~] = FresnelMat(lambdaInVals, thetaVecra, hVec, nMatra, pol); 
Rmatra = mean(cat(3,Rmatrap,Rmatras),3); % average polarizations
Amatra = mean(cat(3,Amatrap,Amatras),3); % averaging here shouldn't matter because we're only looking normal to the sail for absorption, but in principle this should be done.
% Calculate R/T/A for larger wavelength range at all angles to find emissivity
pol = 1;
[~, ~, Amatep, ~, ~] = FresnelMat(lambdaOutVals, thetaVece, hVec, nMate, pol);
pol = 0;
[~, ~, Amates, ~, ~] = FresnelMat(lambdaOutVals, thetaVece, hVec, nMate, pol);
Amate = mean(cat(3,Amatep,Amates),3); % average polarizations
% Absorption spectrum
AbsSpectrum = Amatra(:,1); % Pull out only the 0-degree angle data to look at the center-of-sail temperature
% Angle-dependent reflectivity
thetaMatra = repmat(thetaVecra,[length(lambdaInVals),1]); % pattern over all lambda values
integrandMatR = Rmatra.*sin(thetaMatra).*(cos(thetaMatra)).^3;
RefSpectrum = 8*(sSail/dSailP)^2*sum((1/2)*(integrandMatR(:,1:end-1)+integrandMatR(:,2:end)).*diff(thetaMatra,1,2),2); % ## ... this is doing trapezoidal integration in the 2nd dimension (along rows, integrating over theta)
RefSpectrumRatio = RefSpectrum./Rmatra(:,1); % Ratio of net reflectivity to normal reflectivity.  
RefSpectrumRatioAvg = mean(RefSpectrumRatio); % average decrease in reflectivity due to curvature
% Integrate over theta to find emissivity
AmateEx = cat(2,Amate,zeros(length(lambdaOutVals),1)); % absorption is zero at 90 degrees
integrandMatE = AmateEx.*sin(thetaMateEx).*cos(thetaMateEx); 
EmisSpectrum = (2*pi)*(1/pi)*sum((1/2)*(integrandMatE(:,1:end-1)+integrandMatE(:,2:end)).*diff(thetaMateEx,1,2),2); % ## ... this is doing trapezoidal integration in the 2nd dimension (along rows, across angles)
% Effective emissivity for reference
% Black body spectrum
IbVec = (2*hPlanck*cLight^2)./((lambdaOutVals*(1/1e6)).^5.*(exp((hPlanck*cLight)./(kB*lambdaOutVals*(1/1e6)*Toptical))-1)); % W/m3-sr
integrandVecDen = IbVec; % W/m3-sr
integralDen = sum((1/2)*(integrandVecDen(1:end-1)+integrandVecDen(2:end)).*diff(lambdaOutVals)); % ## ... trapezoidal integration
% https://en.wikipedia.org/wiki/Kirchhoff%27s_law_of_thermal_radiation
integrandVecNum = IbVec.*EmisSpectrum; % W/m3-sr
integralNum = sum((1/2)*(integrandVecNum(1:end-1)+integrandVecNum(2:end)).*diff(lambdaOutVals)); % ## ...trapezoidal integration

Toptical 
disp('note: effective emissivity depends on Toptical!')
effEps = integralNum/integralDen % effective emissivity value design.
avgAbs = mean(Amatra(:,1)) % mean of all 0-degree (normal) absorption values in lambda space... but note that these are spaced evenly in BETA-space not in lambda. 

figure(200)
plot(lambdaOutVals,EmisSpectrum)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


% Find max constant power without failing
% Newton-Raphson loop to solve for maximum constant power that this design can handle for a given failure margin.
% Iteration parameters
maxPwr0 = 500*(1e9/1); % W ... maximum power... beyond this and it's probably an error.
minPwr0 = 1*(1e3/1); % W ... min power... lower than this and it's probably an error. 
maxIter = 100; % max iterations
maxTime = 10; % seconds 
tolerance = 0.0005; % convergence tolerance
deltaPctPwr0 = 0.005; % percent difference for computing derivative
% Starting guesses
Pwr0guessA = min(2*(effEps/avgAbs)*sigmaBB*AsailP*((1+betaVals)./(1-betaVals))*Tmelt^4);
Pwr0guessA = 10e9
Pwr0guessB = Pwr0guessA + deltaPctPwr0*Pwr0guessA % point B guess
% Guess for temperature
TsailGuess = 900; % K
% Iterative loop 
stopper=0;
count = 1;
tStartNR = tic;
while stopper==0
    % for point A
    [failureRatiosA, accValsLaser, accValsSail, Pvals, EdiffVals, EabsVals, TsailEqVals, photonInducedTensionVals, effEpsAtTsailVals] = solveSailConstPwr(...
        Pwr0guessA, TsailGuess, hPlanck, cLight, kB, Mtot, AsailP, sSail, mechRob, Tmelt, ...
        betaVals, gammaVals, RefSpectrum, RefSpectrumRatio, EmisSpectrum, lambdaOutVals, AbsSpectrum);
    [fRA, fRidA] = max(failureRatiosA);
    fA = fRA-(1-failureMargin);
    %if (fA/(1-failureMargin))<=0 && -(fA/(1-failureMargin))<=tolerance ; break; end % this expression forces us to be below the failure value, but it can't by itself guarantee convergence because the algorithm might never jump to the other side. 
    %if abs(fA/(1-failureMargin))<=tolerance; fprintf('(%1.0f-%1.0f) Success: Pwr0=%4.3e \n',m,n,Pwr0guessA); break; end % exit if we're within tolerance
    if abs(fA/(1-failureMargin))<=tolerance; break; end % exit if we're within tolerance
    % for point B.  Use the first T of the first iteration cycle as the guess now. 
    [failureRatiosB, ~, ~, ~, ~, ~, ~, ~, ~] = solveSailConstPwr(...
        Pwr0guessB, TsailEqVals(1), hPlanck, cLight, kB, Mtot, AsailP, sSail, mechRob, Tmelt, ...
        betaVals, gammaVals, RefSpectrum, RefSpectrumRatio, EmisSpectrum, lambdaOutVals, AbsSpectrum);
    [fRB, fRidB] = max(failureRatiosB);
    fB = fRB-(1-failureMargin);
    if(fRidA~=fRidB); fprintf('(%1.0f-%1.0f) Warning: Failure mode changed in differential \n', m,n); end
    % next guesses
    Pwr0guessAold = Pwr0guessA; % store previous A guess value
    Pwr0guessA = Pwr0guessA + (deltaPctPwr0*Pwr0guessA)*((fA)/(fA-fB));
    Pwr0guessB = Pwr0guessA + deltaPctPwr0*Pwr0guessA;
    % now check our next input... don't run away! 
    if Pwr0guessA > maxPwr0
        Pwr0guessA = (Pwr0guessAold+maxPwr0)/2; % take the average 
        Pwr0guessB = Pwr0guessA + deltaPctPwr0*Pwr0guessA;
    elseif Pwr0guessA < minPwr0
        Pwr0guessA = (Pwr0guessAold+minPwr0)/2;
        Pwr0guessB = Pwr0guessA + deltaPctPwr0*Pwr0guessA;
    end
    % prevent a runaway... too many iterations or too much time
    count = count+1;
    if (count > maxIter); fprintf('RUNAWAY: count > maxIter \n'); Pwr0guessA=NaN; break; end
    tElapsedNR = toc(tStartNR);
    if (tElapsedNR > maxTime); fprintf('RUNAWAY: tElapsed > maxTime \n'); Pwr0guessA=NaN; break; end
end
tElapsedNR = toc(tStartNR); % timer off
Pwr0 = Pwr0guessA % W
%%% Check power found against power allowable; if too high, reset power and rerun
if Pwr0 > Pwr0limit
    Pwr0 = Pwr0limit
    fprintf('Recalculating at power limit \n')
    [failureRatiosA, accValsLaser, accValsSail, Pvals, EdiffVals, EabsVals, TsailEqVals, photonInducedTensionVals, effEpsAtTsailVals] = solveSailConstPwr(...
        Pwr0, TsailGuess, hPlanck, cLight, kB, Mtot, AsailP, sSail, mechRob, Tmelt, ...
        betaVals, gammaVals, RefSpectrum, RefSpectrumRatio, EmisSpectrum, lambdaOutVals, AbsSpectrum);
end
%%% Finish Newton mechanics now that we know the laser output power
integrandVec1 = (1./RefSpectrum).*((gammaVals.*betaVals)./(1-betaVals).^2); % ## (no units)
integralVec1 = zeros(size(integrandVec1)); % ## ... the first values of the integral are zero
integralVec1(2:end) = cumsum((1/2)*(integrandVec1(1:end-1)+integrandVec1(2:end)).*diff(betaVals)); % ## ... this is doing cumulative trapezoidal integration 
distVals = (cLight^3/(2*Pwr0))*(Mtot*(1/1000))*integralVec1; % m ... distance values throughout acceleration
integrandVec2 = (gammaVals.^3./RefSpectrum).*((1+betaVals)./(1-betaVals)); % ## (no units)
integralVec2 = zeros(size(integrandVec2)); % ## ... the first values of the integral are zero
integralVec2(2:end) = cumsum((1/2)*(integrandVec2(1:end-1)+integrandVec2(2:end)).*diff(betaVals)); % ## ... this is doing cumulative trapezoidal integration
timeVals = (cLight^2/(2*Pwr0))*(Mtot*(1/1000))*integralVec2; % s ... time values throughout acceleration
timeValsLaser = timeVals-distVals/cLight; % s
PwrValsOnSail = Pwr0*(1-betaVals)./(1+betaVals);

Pwr0gw = Pwr0/1e9 % GW
Tmax = max(TsailEqVals) % K
sailTime = timeVals(end)*(1/60) % min
laserTime = timeValsLaser(end)*(1/60) % min
accelDist = distVals(end)*(1/1e9) % GM
laserDiam = (2*lambda0*(1/1e6)*distVals(end)/dSailP)*(1/1e3) % km
dSailP % m diameter... also spherical radius of curvature
AsailP % m2 perpendicular area


%%% Save ouput
if writeFiles
    fprintf('Writing output\n')
    fNameOut = sprintf('thermalSimulation.csv');
    h = fopen(fNameOut,'w');
    firstLineString = ['TimeSailMin, TimeLaserOnMin, BetaNum, VelocityMps, LambdaUm, LaserOutputPowerGW, PowerOnSailGW, ',...
                        'DistanceGm, AccelLaserMps2, AccelSailMps2, PressurePa, photonInducedTensionNpm,',...
                        'tTemperatureK, ReflectAvgNum, ReflectPerpNum, AbsPerpx1e9Num, ExitanceWpm2, effectiveEmissivityAtTsailNum\n'];
    fprintf(h, firstLineString);
    formatString = ['%6.6e, %6.6e, %6.6e, %6.6e, %6.6e, %6.6e, %6.6e, ',...
                    '%6.6e, %6.6e, %6.6e, %6.6e, %6.6e, ',...
                    '%6.6e, %6.6e, %6.6e, %6.6e, %6.6e, %6.6e \n'];
    for i=1:length(betaVals)
        fprintf(h, formatString, ...
                    timeVals(i)*(1/60), timeValsLaser(i)*(1/60), betaVals(i), betaVals(i)*cLight, lambdaInVals(i), Pwr0*(1/1e9), PwrValsOnSail(i)*(1/1e9), ...
                    distVals(i)*(1/1e9), accValsLaser(i), accValsSail(i), Pvals(i), photonInducedTensionVals(i), ...
                    TsailEqVals(i), RefSpectrum(i), Rmatra(i,1), AbsSpectrum(i)*1e9, EabsVals(i), effEpsAtTsailVals(i));
    end
    fclose('all');
else
    fprintf('NOT WRITING OUTPUT!\n')
end


%addpath('./DisplayTools/');
%TileFigures;



%% Functions


function [failureRatios, accValsLaser, accValsSail, Pvals, EdiffVals, EabsVals, TsailEqVals, photonInducedTensionVals, effEpsAtTsailVals] = solveSailConstPwr(...
            Pwr0, Tsail1Guess, hPlanck, cLight, kB, Mtot, AsailP, sSail, mechRob, Tmelt,  ...
            betaVals, gammaVals, RefSpectrum, RefSpectrumRatio, EmisSpectrum, lambdaOutVals, AbsSpectrum)
    % This solves for the sail's acceleration parameters given a laser power value [W]
    % 
    % Newtonian mechanics
        accValsLaser = ((2*Pwr0*RefSpectrum)./(Mtot*(1/1000)*cLight*gammaVals.^3)).*((1-betaVals)./(1+betaVals)); % m/s2 laser ref frame
        accValsSail = accValsLaser.*gammaVals.^3; % m/s2 sail ref frame
        Pvals = ((1./RefSpectrumRatio)*(Mtot*(1/1000)).*accValsSail)/AsailP; % (kg*m/s2)/m2 = Pa ... photon pressure on sail (occurring in the center).  The reflectivity ratio ensures that the force felt here is the local/central area value, which is higher than the average value on the entire sail
    % Energy balance
        % Energy in at all times through acceleration.  Looking at differential area at center of sail. 
        intensityInVec = ((1-betaVals)./(1+betaVals))*Pwr0/AsailP; % W/m2 ... the first term adjusts for the lower energy of red-shifted photons
        energyInVec = intensityInVec.*AbsSpectrum; % W/m2 ... absorbed power from laser as fxn of velocity (ie, beta). [both column vectors]
        % Energy output.  Need to loop over all trajectory times and iterate to solve for temperature at each one
        % Iteration parameters
        maxT = 10000; % K
        minT = 3; % K
        maxIter = 100; % max iterations
        maxTime = 3; % seconds 
        tolerance = 0.00005; % convergence tolerance. This needs to be very tight because it's the inner NR loop.  If it's too loose, it could throw off the outer NR loop. 
        deltaPctT = 0.005; % percent difference for computing derivative
        % Prefill
        EdiffVals = zeros(size(betaVals)); % W/m2. EdifVals<0 means Eout<Ein
        TsailEqVals = zeros(size(betaVals)); % K ... sail temperature as a function of velocity throughout trajectory
        effEpsAtTsailVals = zeros(size(betaVals)); % epsilon_eff at Tsail
        % Initial guess for k=1 case
        TsailGuessA = Tsail1Guess; 
        if isnan(TsailGuessA); TsailGuessA = sqrt(maxT*minT); end % geometric average
        % Loop over all beta (time) values and solve the energy equation at each one 
        for k=1:length(betaVals)
            % Starting guesses
            if k>1; TsailGuessA = TsailEqVals(k-1); end % from previous time step
            if isnan(TsailGuessA) || TsailGuessA<=minT || TsailGuessA>maxT; TsailGuessA = sqrt(maxT*minT); end
            TsailGuessB = TsailGuessA + deltaPctT*TsailGuessA; % point B guess
            % Iterative loop 
            stopper=0;
            count = 1;
            tStartNRT = tic;
            while stopper==0
                % for point A
                [EconsA, effEpsAtTsailA] = solveEnergyOneBetaStep(TsailGuessA, energyInVec(k), hPlanck, cLight, kB, EmisSpectrum, lambdaOutVals); % gives Econs=Eout-Ein value
                fA = EconsA/energyInVec(k); % W/m2. Econs<0 means Eout<Ein and Tguess was too high. 
                if abs(fA)<=tolerance; break; end % exit if we're within tolerance
                % for point B
                [EconsB, ~]= solveEnergyOneBetaStep(TsailGuessB, energyInVec(k), hPlanck, cLight, kB, EmisSpectrum, lambdaOutVals); % gives Econs=Eout-Ein value
                fB = EconsB/energyInVec(k); % W/m2.
                % next guesses
                TsailGuessAold = TsailGuessA; % store previous A guess value
                TsailGuessA = TsailGuessA + (deltaPctT*TsailGuessA)*((fA)/(fA-fB));
                TsailGuessB = TsailGuessA + deltaPctT*TsailGuessA;
                % now check our next input... don't run away! 
                if TsailGuessA > maxT
                    TsailGuessA = (TsailGuessAold+maxT)/2; % take the average 
                    TsailGuessB = TsailGuessA + deltaPctT*TsailGuessA;
                elseif TsailGuessA < minT
                    TsailGuessA = (TsailGuessAold+minT)/2;
                    TsailGuessB = TsailGuessA + deltaPctT*TsailGuessA;
                end
                % prevent a runaway... too many iterations or too much time
                count = count+1;
                if (count > maxIter); fprintf('RUNAWAY: count > maxIter \n'); TsailGuessA=NaN; break; end
                tElapsedNRT = toc(tStartNRT);
                if (tElapsedNRT > maxTime); fprintf('RUNAWAY: tElapsed > maxTime \n'); TsailGuessA=NaN; break; end
            end
            tElapsedNRT = toc(tStartNRT); % timer off
            TsailEqVals(k) = TsailGuessA; % K
            EdiffVals(k) = EconsA; % W/m2
            effEpsAtTsailVals(k) = effEpsAtTsailA; % no units
        end
        EabsVals = energyInVec; % W/m2 ... energy absorbed (= energy emitted) as a function of velocity throughout trajectory
        TsailRatioMax = max(TsailEqVals/Tmelt); 
    % mechanical robustness criterion
        photonInducedTensionVals = Pvals*sSail/2; % N/m 
        mechRatioMax = max(photonInducedTensionVals/mechRob);
    % Failure ratios
        failureRatios = [mechRatioMax TsailRatioMax];
end

function [Econs, effEpsAtTsail] = solveEnergyOneBetaStep(Tsail, energyIn, hPlanck, cLight, kB, EmisSpectrum, lambdaOutVals)
    % Solves the energy equation to find the sail temperature at a single point in time
    IbVec = (2*hPlanck*cLight^2)./((lambdaOutVals*(1/1e6)).^5.*(exp((hPlanck*cLight)./(kB*lambdaOutVals*(1/1e6).*Tsail))-1)); % W/m3-sr
    integrandVecEmis = pi*EmisSpectrum.*IbVec*(1/1e6); % W/m2-um 
    integralEmis = sum((1/2)*(integrandVecEmis(1:end-1)+integrandVecEmis(2:end)).*diff(lambdaOutVals)); % trapezoidal integration in lambda space
    energyOut = 2*integralEmis; % W/m2 ... emission can occur from both sides, hence the factor of 2. 
    Econs = energyOut-energyIn; % W/m2

    % Effective emissivity at this temperature (using constant optical properties)
    % Black body spectrum
    integrandVecDen = IbVec; % W/m3-sr
    integralDen = sum((1/2)*(integrandVecDen(1:end-1)+integrandVecDen(2:end)).*diff(lambdaOutVals)); % ## ... trapezoidal integration
    % https://en.wikipedia.org/wiki/Kirchhoff%27s_law_of_thermal_radiation
    integrandVecNum = IbVec.*EmisSpectrum; % W/m3-sr
    integralNum = sum((1/2)*(integrandVecNum(1:end-1)+integrandVecNum(2:end)).*diff(lambdaOutVals)); % ## ...trapezoidal integration
    effEpsAtTsail = integralNum/integralDen; % effective emissivity value design.
end







