% Find optimal Al2O3-MoS2-Al2O3 film thicknesses
% We're optimizing betaMax
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

% Calculate parameters for each sail
constraintSelect = 1;
saveContourMats = 0;
saveRefData = 0;

% Constants
cLight = 299792458; % m/s ... light speed
hPlanck = 6.62607015e-34; % J-s ... Planck's constant
kB = 1.380649e-23; % J/K ... Boltzmann constant
sigmaBB = (2*pi^5*kB^4)/(15*cLight^2*hPlanck^3); % W/m2-K4

% Inputs
maxPwr0 = 100e9; % W ... maximum power producible by light engine
dLaser0 = 30e3; % m ... laser array diameter on Earth 
beta0 = 1e-10; % start from rest... Beta = v/cLight so 1e-10 is ~ 0.03 m/s
lambda0 = 1.2; % um ... laser wavelength at zero sail velocity.  
theta = 0; % normal incidence
hexD = 70e-6; % m ... hexagonal corrugation parameters
hexW = 2e-6; % m 
hexH = 3e-6; % m 

% Material parameters
rhoAl2O3 = 3200; % kg/m3
rhoMoS2 = 5060; % kg/m3
sigYldAl2O3 = 2e9; % Pa ... was 5e9
sigYldMoS2 = 2.3e9; % Pa ... was 0.75e9

% thickness values
% Fine grid to test - takes a long time to run
tAvec = linspace(1,100,300)*(1/1e9); % nm --> m ... alumina face layer thickness (this value on both sides, so the total is double this)
tMvec = linspace(1,200,301)*(1/1e9); % nm --> m ... mos2 thickness
% Rough grid, faster. 
tAvec = linspace(1,100,30)*(1/1e9); % nm --> m ... alumina face layer thickness (this value on both sides, so the total is double this)
tMvec = linspace(1,200,31)*(1/1e9); % nm --> m ... mos2 thickness
[tAmat, tMmat] = ndgrid(tAvec, tMvec); 
numCases = numel(tAmat)

% Read Kischkat 2012 Al2O3 data
fNameK2012 = 'Kischkat2012_Al2O3nk.csv'; 
tempMatrixK2012 = csvread(fNameK2012,1,0); % read the CSV file into a matrix
lambdaValsK2012 = tempMatrixK2012(1:end,1); % um ... first column is wavelength
spectrumK2012n = tempMatrixK2012(1:end,2); % n
spectrumK2012k = tempMatrixK2012(1:end,3); % k
interpFxnKischkatN = griddedInterpolant(lambdaValsK2012,spectrumK2012n,'linear','linear'); % create function to describe data and enable interpolations.  
interpFxnKischkatK = griddedInterpolant(lambdaValsK2012,spectrumK2012k,'linear','linear'); % create function to describe data and enable interpolations.  
pKischkat2012smallExtrap = polyfit(lambdaValsK2012(1:791),spectrumK2012k(1:791),3);

% Read Munkhbat 2022 MoS2 data
fNameMunk = 'Munkhbat2022MoS2nk.csv'; %fprintf('Using Single Layer MoS2 Data\n');
tempMatrixMunk = csvread(fNameMunk,2,0); % read the CSV file into a matrix
lambdaValsMunk = tempMatrixMunk(:,1); % um ... first column is wavelength
spectrumMunkN = tempMatrixMunk(:,2); % arbitrary units, (n) index of refraction data @ 300 K
spectrumMunkK = tempMatrixMunk(:,3); % arbitrary units, (k) extinction coefficient data @ 300 K
interpFxnMoS2n = griddedInterpolant(lambdaValsMunk,spectrumMunkN,'linear','linear'); % create function to describe data and enable interpolations.  
interpFxnMoS2k = griddedInterpolant(lambdaValsMunk,spectrumMunkK,'linear','linear'); % create function to describe data and enable interpolations.  




% areal density and mechanical robustness
% trenches
MucMat = ...
    rhoAl2O3 * (tAmat*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD+2*tAmat).^2-hexD^2)) + ...
    rhoMoS2 * (tMmat*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD+2*tAmat+2*tMmat).^2-(hexD+2*tAmat).^2)) + ...
    rhoAl2O3 * (tAmat*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD+2*tAmat+2*tMmat+2*tAmat).^2-(hexD+2*tAmat+2*tMmat).^2)); % kg ... unit cell mass
% ribs
% MucMat = ...
%     rhoAl2O3 * (tAmat*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*(hexD^2-(hexD-2*tAmat).^2)) + ...
%     rhoMoS2 * (tMmat*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD-2*tAmat).^2-(hexD-2*tAmat-2*tMmat).^2)) + ...
%     rhoAl2O3 * (tAmat*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD-2*tAmat-2*tMmat).^2-(hexD-2*tAmat-2*tMmat-2*tAmat).^2)); % kg ... unit cell mass
Auc = sqrt(3)/2 * (hexD+hexW)^2; % unit cell area
rhoAcMat = MucMat/Auc *(1000/1); % g/m2 ... Corrugated areal density
G = 1;
mechRobMat = G * (tAmat*sigYldAl2O3 + tMmat*sigYldMoS2 + tAmat*sigYldAl2O3); % N/m 

% now calculate beta max

% Sail diameter... assumes spherically curved sail with s=d (radius of curvature equal to diameter)
if constraintSelect == 1 % specify total mass
    constraintString = "mass";
    MtotMat = 2*ones(size(tAmat)); % g ... total mass of sail, tethers, and chip
    MsailMat = MtotMat/2; % g ... mass of sail only
    AsailMat = MsailMat./rhoAcMat; % m2
    sSailMat = sqrt(AsailMat/(pi*(2-sqrt(3)))); % m ... equation for radius of sail when s=d ... see Equation S44 in Campbell 2022 sup inf
    dSailMat = sSailMat; % m ... sail diameter 
    AsailPerpMat = pi*(dSailMat/2).^2; % m
end

if constraintSelect == 2 % specify sail diameter 
    constraintString = "diameter";
    dSailMat = 2*ones(size(tAmat)); % m
    sSailMat = dSailMat; % guideline from Campbell 2022
    AsailPerpMat = pi*(dSailMat/2).^2; % m
    AsailMat = pi*(2-sqrt(3))*sSail.^2; % m ... equation for spherically curved sail area when s=d ... see Equation S44 in Campbell 2022 sup inf
    MsailMat = AsailMat*rhoAcMat; % g
    MtotMat = 2*MsailMat;
end

% Estmate maximum acceleration distance possible with focus on sail, given sail diameter and laser array diameter and wavelength
LmaxMat = dSailMat*dLaser0/(2*lambda0*(1/1e6)); % m

% Light wavelengths at the sail
lambdaMax = 2; % um
digitsOld = digits(10); % reduce precision temporarily
syms betaFind; % declare variable
betaF = vpasolve(lambdaMax == lambda0*(1+betaFind)/sqrt(1-betaFind^2), betaFind, 0.2); % upper limit to the fraction of speed of light that we want to achieve... if sails could get going this fast. 
digits(digitsOld); % restore precision
betaF = double(betaF); % increase to double precision
betaVals = linspace(beta0, betaF, 6000)'; % fraction of speed of light
gammaVals = 1./sqrt(1-betaVals.^2); % relativistic Lorentz factor
lambdaInVals = lambda0 * sqrt((1+betaVals)./(1-betaVals)); % um ... wavelengths of light that hit the sail

% find index cooresponding to beta=0.2
betaWantStarshot = 0.2; 
[dummy,idxBetaStarshot] = min(abs(betaVals-betaWantStarshot));
betaHaveStarshot = betaVals(idxBetaStarshot);

% optical properties
AnData = interpFxnKischkatN(lambdaInVals);
AkData = interpFxnKischkatK(lambdaInVals);
MnData = interpFxnMoS2n(lambdaInVals); 
MkData = interpFxnMoS2k(lambdaInVals);


figure(100)
plot(lambdaInVals,AnData)

figure(101)
hold on
plot(lambdaInVals,AkData)
plot(lambdaValsK2012,spectrumK2012k)
hold off


% Assemble complex index of refraction
% BE CAREFUL WHEN ASSEMBLING COMPLEX NUMBERS!  MAKE SURE ALL ARE COLUMN VECTORS!!!
if ~iscolumn(lambdaInVals); warning('check complex vectors'); end
indComplexInA = AnData + 1i*AkData; % face layers ... here lambda increases down columns 
indComplexInM = MnData + 1i*MkData; % middle later ... here lambda increases down columns
nMatIn = [ones(length(lambdaInVals),1)  indComplexInA  indComplexInM  indComplexInA  ones(length(lambdaInVals),1)]; % the index of refraction values depend on temperature! In nMat, lambda increases down columns and layer number increases across rows.

% Iteration parameters... want to find the power that takes the sail within "failureMargin" of breaking. 
failureMargin = 0.0001; % 0.1-percent failure margin
maxPwr0nr = 1e15; % W ... maximum power for Newton Raphson algorithm... beyond this and it's probably an error.
minPwr0nr = 1e6; % W ... min power... lower than this and it's probably an error. 
maxIterPwr0 = 100; % max iterations
maxTimePwr0 = 3; % seconds 
tolerancePwr0 = 0.0005; % convergence tolerance
deltaPctPwr0 = 0.005; % percent difference for computing derivative

% Loop over thickness matricees and iteratively solve for betaMax for each.
betaMaxMat = zeros(size(tAmat)); % prefill
Pwr0mat = betaMaxMat;
Pwr0foundMat = betaMaxMat;
LfinalMat = betaMaxMat;
RavgMat = betaMaxMat;
AavgMat = betaMaxMat;
TavgMat = betaMaxMat;
betaLimitMat = -1*ones(size(betaMaxMat)); % keep track of the beta where the power was limited; need this to be less than betaMax
betaLimitIndexMat = betaMaxMat; % full of zeros
RefSpectrumMat = zeros(size(tAmat,1),size(tAmat,2),length(lambdaInVals));
for m=1:size(tAmat,1)
    fprintf('m=%4.0f\n',m)
    for n=1:size(tAmat,2)
        % where are we?
        %fprintf('m = %4.0f n=%4.0f\n',m,n)
        keepCalculating = true;
        indexQuit = length(lambdaInVals);
        while keepCalculating == true
            keepCalculating = false;
            % thicknesses and sail parameters
            tA = tAmat(m,n); % m
            tM = tMmat(m,n); % m
            Mtot = MtotMat(m,n); % g
            sSail = sSailMat(m,n); % m
            AsailPerp = AsailPerpMat(m,n); % m2
            mechRob = mechRobMat(m,n); % N/m
            Lmax = LmaxMat(m,n); % m
            % Layer thicknesses
            hVec = [NaN, tA, tM, tA, NaN]*(1e6/1); % m --> um... film layer thicknesses
            % Run Transfer Matrix
            pol = 1;    
            [Rvecp, Tvecp, Avecp, Mdummy, Vdummy] = FresnelMat(lambdaInVals(1:indexQuit), theta, hVec, nMatIn(1:indexQuit,:), pol); 
            pol = 0;
            [Rvecs, Tvecs, Avecs, Mdummy, Vdummy] = FresnelMat(lambdaInVals(1:indexQuit), theta, hVec, nMatIn(1:indexQuit,:), pol); 
            Rvec = mean(cat(2,Rvecp,Rvecs),2); % average polarizations
            Avec = mean(cat(2,Avecp,Avecs),2); % 
            Tvec = mean(cat(2,Tvecp,Tvecs),2); % 
            RefSpectrum = Rvec;
            % Initial power guesses
            if m==1 && n==1
                Pwr0guessA = 1e9; % W = kg-m2/s3 ... laser power
            elseif m==1 && n>1
                Pwr0guessA = Pwr0foundMat(m,n-1); % W
            elseif m>1 && n>1
                Pwr0guessA = mean([Pwr0foundMat(m-1,n) Pwr0foundMat(m,n-1)]); % W
            end
            Pwr0guessB = Pwr0guessA + deltaPctPwr0*Pwr0guessA; % point B guess
            % Iterative loop 
            stopper=0;
            count = 1;
            tStartNR = tic;
            while stopper==0
                % for point A
                [mechRatioMaxA, distValsA, timeValsA, accValsLaserA, accValsSailA, PvalsA, breakingForceValsA] = solveSail(...
                        Pwr0guessA, cLight, Mtot, sSail, AsailPerp, Lmax, mechRob, betaVals(1:indexQuit), gammaVals(1:indexQuit), RefSpectrum);
                fRA = mechRatioMaxA;
                fA = fRA-(1-failureMargin);
                if abs(fA/(1-failureMargin))<=tolerancePwr0; break; end % exit if we're within tolerance
                % for point B.  
                [mechRatioMaxB, ~, ~, ~, ~, ~, ~] = solveSail(...
                        Pwr0guessB, cLight, Mtot, sSail, AsailPerp, Lmax, mechRob, betaVals(1:indexQuit), gammaVals(1:indexQuit), RefSpectrum);
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
            actualPower = Pwr0found; % W
            % Check if we went over the allowed power
            usingMaxPower = 0;
            if Pwr0found>maxPwr0 % if a greater power was allowable, than surely a lower power would be as well.  Just use the max allowable. 
                [mechRatioMaxA, distValsA, timeValsA, accValsLaserA, accValsSailA, PvalsA, breakingForceValsA] = solveSail(...
                        maxPwr0, cLight, Mtot, sSail, AsailPerp, Lmax, mechRob, betaVals(1:indexQuit), gammaVals(1:indexQuit), RefSpectrum);
                usingMaxPower = 1;
                actualPower = maxPwr0; % W
                betaLimitMat(m,n) = -1;
                betaLimitIndexMat(m,n) = 0;
            else % if we used the max power or less than the max
                % check which Beta limited the power.  We want to make sure this beta is
                % less than or equal to betamax, otherwise the limit occurs outside of the
                % actual practical range and we could have had higher power. 
                ratioVec = breakingForceValsA/mechRob;
                indexLimit = find(ratioVec == max(ratioVec),1);
                betaLimitMat(m,n) = betaVals(indexLimit);
                betaLimitIndexMat(m,n) = indexLimit;
            end
            % find betaMax
            tempVec = find(isnan(distValsA)); % vector of indicees corresponding to NaN, which are values where distance is too long 
            if(isempty(tempVec)) % if there's no NaN values, than the last value is the limit Lmax. 
                indexLast = length(distValsA);
            else
                indexLast = tempVec(1)-1; % last index for which the acceleration distance was less than the diffraction limit of distance
            end
            betaMax = betaVals(indexLast);
            Lfinal = distValsA(indexLast); % m
            % check that betaMax is greater than power-limiting beta
            if ~(betaLimitMat(m,n)<=betaMax) % if the limiting beta is not less than or equal to betamax.  Note that the matrix was prefilled with -1. 
                %fprintf('Beta limit out of range: m=%1.0f n=%1.0f betaLim=%5.4f betaMax=%5.4f indexLimit=%2.0f indexQuit=%2.0f indexLast=%2.0f\n',m,n,betaLimitMat(m,n),betaMax,indexLimit,indexQuit,indexLast);
                indexQuit = indexQuit - 10; % take off a beta value and do it again... keep doing this until the limit is within range. 
                keepCalculating = true;
            end
            if indexQuit == length(lambdaInVals) % doing this here stores the full vectors without being cut off.  
                RefSpectrumMat(m,n,:) = Rvec; % save for later
                % The code here averages over beta space. 
                %RavgMat(m,n) = mean(Rvec(1:idxBetaStarshot)); % average the spectrum over the starshot range
                %AavgMat(m,n) = mean(Avec(1:idxBetaStarshot)); 
                %TavgMat(m,n) = mean(Tvec(1:idxBetaStarshot));
                % A better way is to average over lambda space via trapezoidal integration, which is how I calculate average optical properties everywhere else. 
                RavgMat(m,n) = sum((1/2)*(Rvec(1:idxBetaStarshot-1)+Rvec(2:idxBetaStarshot)).*(lambdaInVals(2:idxBetaStarshot)-lambdaInVals(1:idxBetaStarshot-1)))/(lambdaInVals(idxBetaStarshot)-lambdaInVals(1));
                AavgMat(m,n) = sum((1/2)*(Avec(1:idxBetaStarshot-1)+Avec(2:idxBetaStarshot)).*(lambdaInVals(2:idxBetaStarshot)-lambdaInVals(1:idxBetaStarshot-1)))/(lambdaInVals(idxBetaStarshot)-lambdaInVals(1));
                TavgMat(m,n) = sum((1/2)*(Tvec(1:idxBetaStarshot-1)+Tvec(2:idxBetaStarshot)).*(lambdaInVals(2:idxBetaStarshot)-lambdaInVals(1:idxBetaStarshot-1)))/(lambdaInVals(idxBetaStarshot)-lambdaInVals(1));
            end
        end
        % store data
        betaMaxMat(m,n) = betaMax;
        LfinalMat(m,n) = Lfinal; % m
        Pwr0mat(m,n) = actualPower; % W
        Pwr0foundMat(m,n) = Pwr0found; % W
    end
end

% save output
diaryName = 'commandWindowOutput.txt';
diary(diaryName)
disp(datestr(now))
diary on

idxBetaStarshot
betaHaveStarshot

% find thickness with highest betaMax value
optimalBetaMaxLinearIndexVals = find(betaMaxMat==max(max(betaMaxMat))) % linear index; increases down columns then across rows.
optimalBetaMaxLinearIndex = optimalBetaMaxLinearIndexVals(ceil(length(optimalBetaMaxLinearIndexVals)/2)) % if there's more than one, take the middle value in the vector. Note use of "ceil" so that if there's only one value, you still get an answer! 
betaMaxOptimal = betaMaxMat(optimalBetaMaxLinearIndex) % fractoin of speed of light
tAoptimal = tAmat(optimalBetaMaxLinearIndex)*(1e9/1) % nm
tMoptimal = tMmat(optimalBetaMaxLinearIndex)*(1e9/1) % nm
Pwr0optimal = Pwr0mat(optimalBetaMaxLinearIndex)*(1/1e9) % GW
FracMaxDistOptimal = LfinalMat(optimalBetaMaxLinearIndex)/LmaxMat(optimalBetaMaxLinearIndex) % fraction
rhoAcOptimal = rhoAcMat(optimalBetaMaxLinearIndex) % g/m2
mechRobOptimal = mechRobMat(optimalBetaMaxLinearIndex) % N/m
[mOpt, nOpt] = ind2sub(size(betaMaxMat),optimalBetaMaxLinearIndex);

% Calculate average optical properties for best sail design
% wavelengths for average optical props
betaF = 0.2; % 
betaVals = linspace(beta0, betaF, 5000)'; % fraction of speed of light... column vector! 
gammaVals = 1./sqrt(1-betaVals.^2); % relativistic Lorentz factor
lambdaVals = lambda0 * sqrt((1+betaVals)./(1-betaVals)); % um ... wavelengths of light that hit the sail

% For numerical consistency, use evenly-spaced lambda points rather than
% evenly-spaced beta points
lambdaMaxStarshot = lambda0 * sqrt((1+betaF)./(1-betaF)); % um
lambdaVals = linspace(lambda0, lambdaMaxStarshot, 5000)'; % um

% optical properties within beta = 0-0.2 range
AnData = interpFxnKischkatN(lambdaVals);
AkData = interpFxnKischkatK(lambdaVals);
MnData = interpFxnMoS2n(lambdaVals); 
MkData = interpFxnMoS2k(lambdaVals);

% Assemble complex index of refraction
% BE CAREFUL WHEN ASSEMBLING COMPLEX NUMBERS!  MAKE SURE ALL ARE COLUMN VECTORS!!!
if ~iscolumn(lambdaVals); warning('check complex vectors'); end
indComplexInA = AnData + 1i*AkData; % face layers ... here lambda increases down columns 
indComplexInM = MnData + 1i*MkData; % middle later ... here lambda increases down columns
nMatIn = [ones(length(lambdaVals),1)  indComplexInA  indComplexInM  indComplexInA  ones(length(lambdaVals),1)]; % the index of refraction values depend on temperature! In nMat, lambda increases down columns and layer number increases across rows.

% Layer thicknesses
hVec = [NaN, tAoptimal*(1/1e9), tMoptimal*(1/1e9), tAoptimal*(1/1e9), NaN]*(1e6/1); % m --> um... film layer thicknesses
% Run Transfer Matrix
pol = 1;    
[Rvecp, Tvecp, Avecp, Mdummy, Vdummy] = FresnelMat(lambdaVals, theta, hVec, nMatIn, pol); 
pol = 0;
[Rvecs, Tvecs, Avecs, Mdummy, Vdummy] = FresnelMat(lambdaVals, theta, hVec, nMatIn, pol); 
Rvec = mean(cat(2,Rvecp,Rvecs),2); % average polarizations
Avec = mean(cat(2,Avecp,Avecs),2); % 
Tvec = mean(cat(2,Tvecp,Tvecs),2); % 
RavgOptimal = mean(Rvec) % average the spectrum over beta = [0, 0.2]
AavgOptimal = mean(Avec) 
TavgOptimal = mean(Tvec)

diary off 



figure(9)
hold on;
plot(lambdaVals*(1000/1),Rvec,'g-')
plot(lambdaVals*(1000/1),Avec,'r-')
plot(lambdaVals*(1000/1),Tvec,'b-')
hold off;
xlabel('Wavelength [nm]')
ylabel('R, A, T')
saveas(gcf,'09.png')

% output contour plot data
% Top row is TOTAL alumina thickness (sum of both sides)
% First column is MoS2 thickness
if saveContourMats % these are formatted for Veusz
    fNameOut = sprintf('optimalFilmBetaMax.csv'); 
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',2*tAvec*(1e9/1)); fprintf(h,'\n'); % header line.  Note the initial blank (leading comma) 
    for n=length(tMvec):-1:1 % MoS2 thickness
        fprintf(h,'%7.5e',tMvec(n)*(1e9/1)); fprintf(h,' , %7.5e',betaMaxMat(:,n)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = sprintf('optimalFilmPwr0.csv');  % GW
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',2*tAvec*(1e9/1)); fprintf(h,'\n'); % header line.  Note the initial blank (leading comma) 
    for n=length(tMvec):-1:1 % MoS2 thickness
        fprintf(h,'%7.5e',tMvec(n)*(1e9/1)); fprintf(h,' , %7.5e',Pwr0mat(:,n)*(1/1e9)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = sprintf('optimalFracMaxDist.csv');  % GW
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',2*tAvec*(1e9/1)); fprintf(h,'\n'); % header line.  Note the initial blank (leading comma) 
    for n=length(tMvec):-1:1 % MoS2 thickness
        fprintf(h,'%7.5e',tMvec(n)*(1e9/1)); fprintf(h,' , %7.5e',LfinalMat(:,n)./LmaxMat(:,n)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = sprintf('optimalAccelDist.csv');  % GW
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',2*tAvec*(1e9/1)); fprintf(h,'\n'); % header line.  Note the initial blank (leading comma) 
    for n=length(tMvec):-1:1 % MoS2 thickness
        fprintf(h,'%7.5e',tMvec(n)*(1e9/1)); fprintf(h,' , %7.5e',LfinalMat(:,n)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = sprintf('optimalRavg.csv');  % GW
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',2*tAvec*(1e9/1)); fprintf(h,'\n'); % header line.  Note the initial blank (leading comma) 
    for n=length(tMvec):-1:1 % MoS2 thickness
        fprintf(h,'%7.5e',tMvec(n)*(1e9/1)); fprintf(h,' , %7.5e',RavgMat(:,n)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = sprintf('optimalAavg.csv');  % GW
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',2*tAvec*(1e9/1)); fprintf(h,'\n'); % header line.  Note the initial blank (leading comma) 
    for n=length(tMvec):-1:1 % MoS2 thickness
        fprintf(h,'%7.5e',tMvec(n)*(1e9/1)); fprintf(h,' , %7.5e',AavgMat(:,n)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = sprintf('optimalRhoAc.csv');  % GW
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',2*tAvec*(1e9/1)); fprintf(h,'\n'); % header line.  Note the initial blank (leading comma) 
    for n=length(tMvec):-1:1 % MoS2 thickness
        fprintf(h,'%7.5e',tMvec(n)*(1e9/1)); fprintf(h,' , %7.5e',rhoAcMat(:,n)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = sprintf('optimalMechRob.csv');  % GW
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',2*tAvec*(1e9/1)); fprintf(h,'\n'); % header line.  Note the initial blank (leading comma) 
    for n=length(tMvec):-1:1 % MoS2 thickness
        fprintf(h,'%7.5e',tMvec(n)*(1e9/1)); fprintf(h,' , %7.5e',mechRobMat(:,n)); fprintf(h,'\n'); % next series of lines
    end
end

% output reflectivity data
if saveRefData
    fNameOut = 'ThisWorkOptimalDesignRAT.csv';
    headerString = "WLnm, R, A, T\n";
    h = fopen(fNameOut,'w');
    fprintf(h,headerString);
    for m=1:length(lambdaVals)
        fprintf(h,'%10.6e, %10.6e, %10.6e, %10.6e\n',lambdaVals(m)*(1000/1),Rvec(m),Avec(m),Tvec(m));
    end
end
fclose('all');



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


