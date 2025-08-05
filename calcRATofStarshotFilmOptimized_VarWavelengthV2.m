% Calculate R-A-T of laminate film
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

% film thicknesses
%TfilmA = 0.0089465; % um ... Alumina thickness on both sides of MoS2
%TfilmM = 0.0746300 ; % um ... MoS2

TfilmA = 18.5485e-3; % um
TfilmM = 62.6900e-3;

% wavelength range for starshot with beta=0.2
lambda0 = 1.2; % um
betaFstarshot = 0.2;
lambdaMaxStarshot = lambda0 * sqrt((1+betaFstarshot)./(1-betaFstarshot));
lambdaValsAvgRef = linspace(lambda0,lambdaMaxStarshot,5000)'; % um

% wavelength range
lambdaVals = linspace(1.0,1.687,10000)'; % um
thetaIn = 0; % radians... normal to surface

% Optical data
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

% calculate optical constants
AnData = interpFxnKischkatN(lambdaVals); % arbs ... find n at the wavelengths corresponding to the beta values we defined above
AkData = interpFxnKischkatK(lambdaVals);
MnData = interpFxnMoS2n(lambdaVals); 
MkData = interpFxnMoS2k(lambdaVals);

% calculate optical constants in lambda range for starshot
AnDataAvgRef = interpFxnKischkatN(lambdaValsAvgRef); % arbs ... find n at the wavelengths corresponding to the beta values we defined above
AkDataAvgRef = interpFxnKischkatK(lambdaValsAvgRef);
MnDataAvgRef = interpFxnMoS2n(lambdaValsAvgRef); 
MkDataAvgRef = interpFxnMoS2k(lambdaValsAvgRef);

% Assemble complex index of refraction
% BE CAREFUL WHEN ASSEMBLING COMPLEX NUMBERS!  MAKE SURE ALL ARE COLUMN VECTORS!!!
indComplexInA = AnData + 1i*AkData; % face layers ... here lambda increases down columns 
indComplexInM = MnData + 1i*MkData; % middle later ... here lambda increases down columns
nMatIn = [ones(length(lambdaVals),1)  indComplexInA  indComplexInM  indComplexInA  ones(length(lambdaVals),1)]; % the index of refraction values depend on temperature! In nMat, lambda increases down columns and layer number increases across rows.

indComplexInAavgRef = AnDataAvgRef + 1i*AkDataAvgRef; % face layers ... here lambda increases down columns 
indComplexInMavgRef = MnDataAvgRef + 1i*MkDataAvgRef; % middle later ... here lambda increases down columns
nMatInAvgRef = [ones(length(lambdaValsAvgRef),1)  indComplexInAavgRef  indComplexInMavgRef  indComplexInAavgRef  ones(length(lambdaValsAvgRef),1)]; % the index of refraction values depend on temperature! In nMat, lambda increases down columns and layer number increases across rows.

% Layer thicknesses
hVec = [NaN, TfilmA, TfilmM, TfilmA, NaN]; % um... film layer thicknesses

% Run Transfer Matrix
pol = 1;    
[Rvecp, Tvecp, Avecp, Mdummy, Vdummy] = FresnelMat(lambdaVals, thetaIn, hVec, nMatIn, pol); 
pol = 0;
[Rvecs, Tvecs, Avecs, Mdummy, Vdummy] = FresnelMat(lambdaVals, thetaIn, hVec, nMatIn, pol); 
Rvec = mean(cat(2,Rvecp,Rvecs),2); % average
Avec = mean(cat(2,Avecp,Avecs),2); % average
Tvec = mean(cat(2,Tvecp,Tvecs),2); % average

pol = 1;    
[RvecpAvg, TvecpAvg, AvecpAvg, Mdummy, Vdummy] = FresnelMat(lambdaValsAvgRef, thetaIn, hVec, nMatInAvgRef, pol); 
pol = 0;
[RvecsAvg, TvecsAvg, AvecsAvg, Mdummy, Vdummy] = FresnelMat(lambdaValsAvgRef, thetaIn, hVec, nMatInAvgRef, pol); 
RvecAvg = mean(cat(2,RvecpAvg,RvecsAvg),2); % average
AvecAvg = mean(cat(2,AvecpAvg,AvecsAvg),2); % average
TvecAvg = mean(cat(2,TvecpAvg,TvecsAvg),2); % average
Ravg = mean(RvecAvg)
Aavg = mean(AvecAvg)
Tavg = mean(TvecAvg)

% Save new spectrum over larger range
writeFiles = 0;
if writeFiles == 1
    if opticalData==1 || opticalData==4 || opticalData==5
        refIndSrc = 'literatureNK';
    elseif opticalData==2
        refIndSrc = 'measuredNKjason';
    elseif opticalData==3
        refIndSrc = 'measuredNKpawan';
    end
    fNameOut = ['optimizedSpectrumFull_',refIndSrc,'.csv'];
    headerString = "simWLnm, simR, simA, simT\n";
    h = fopen(fNameOut,'w');
    fprintf(h,headerString);
    for m=1:length(lambdaVals)
        fprintf(h,'%10.6e, %10.6e, %10.6e, %10.6e\n',lambdaVals(m)*(1000/1),Rvec(m),Avec(m),Tvec(m));
    end
    fclose('all');
end

% plot
figure(1)
hold on
plot(lambdaVals*(1000/1),Rvec,'g-','LineWidth', 2)
plot(lambdaVals*(1000/1),Avec,'r-','LineWidth', 2)
plot(lambdaVals*(1000/1),Tvec,'b-','LineWidth', 2)
hold off
set(gca, 'YScale', 'log')
xlabel('Wavelength [nm]')
ylabel('R, A, T [#]')
title('R, A, T');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 1.5);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);


% calculate extrapolation range from Kischkat 2012
lambdaValsExtrap = linspace(1.2,lambdaValsK2012(1),100)'; % um
AnDataExtrap = interpFxnKischkatN(lambdaValsExtrap); % 
AkDataExtrap = interpFxnKischkatK(lambdaValsExtrap);
AkDataPoly = polyval(pKischkat2012smallExtrap,lambdaValsExtrap);

lambdaValsRight = linspace(lambdaValsK2012(1),5,100)'; % um
AnDataRight = interpFxnKischkatN(lambdaValsRight); % 
AkDataRight = interpFxnKischkatK(lambdaValsRight);

figure(10)
hold on
plot(lambdaValsExtrap,AnDataExtrap)
plot(lambdaValsRight,AnDataRight )
hold off

figure(11)
hold on
plot(lambdaValsExtrap,AkDataExtrap)
plot(lambdaValsRight,AkDataRight)
plot(lambdaValsExtrap,AkDataPoly)
hold off
set(gca, 'YScale', 'log')
