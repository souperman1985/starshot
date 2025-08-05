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


writeFiles = 0;

% Constants
cLight = 299792458; % m/s ... light speed
hPlanck = 6.62607015e-34; % J-s ... Planck's constant
hPlanckBar = 1.054571817e-34; % J-s ... Planck's constant with a bar... h/(2*pi)
kB = 1.380649e-23; % J/K ... Boltzmann constant
sigmaBB = (2*pi^5*kB^4)/(15*cLight^2*hPlanck^3); % W/m2-K4

% film thicknesses
TfilmAtop = 0.051; % um ... Alumina thickness on both sides of MoS2
TfilmM = 0.053; % um ... MoS2
TfilmAbot = 0.021; % um ... Alumina thickness on both sides of MoS2

% wavelength range for starshot with beta=0.2
lambda0 = 1.2; % um
betaFstarshot = 0.2;
lambdaMaxStarshot = lambda0 * sqrt((1+betaFstarshot)./(1-betaFstarshot));
lambdaValsAvgRef = linspace(lambda0,lambdaMaxStarshot,1000)'; % um

% wavelength range
lambdaVals = linspace(1.0,1.6,1000)'; % um

thetaIn = 0; % radians... normal to surface

% Optical data
fNameAnk = 'Al2O3nkFromJasonJan2024.csv'; % First column is lambda (nm), second is n, third is kappa
AnkData = csvread(fNameAnk,1,0); % read the CSV file into a matrix
AnkData(:,1) = AnkData(:,1)*(1/1000); % nm --> um
interpFxnAl2O3n = griddedInterpolant(AnkData(:,1),AnkData(:,2),'linear','linear'); % create function to describe data and enable interpolations.  
interpFxnAl2O3k = griddedInterpolant(AnkData(:,1),AnkData(:,3),'linear','linear');

fNameMnk = 'MoS2nkFromJasonApril2024.csv'; % First column is lambda (nm), second is n, third is kappa
MnkData = csvread(fNameMnk,1,0); % read the CSV file into a matrix
MnkData(:,1) = MnkData(:,1)*(1/1000); % nm --> um
interpFxnMoS2n = griddedInterpolant(MnkData(:,1),MnkData(:,2),'linear','linear'); % create function to describe data and enable interpolations.  
interpFxnMoS2k = griddedInterpolant(MnkData(:,1),MnkData(:,3),'linear','linear');

% calculate optical constants
AnData = interpFxnAl2O3n(lambdaVals); % arbs ... find n at the wavelengths corresponding to the beta values we defined above
AkData = interpFxnAl2O3k(lambdaVals);
MnData = interpFxnMoS2n(lambdaVals); 
MkData = interpFxnMoS2k(lambdaVals);

% calculate optical constants in lambda range for starshot
AnDataAvgRef = interpFxnAl2O3n(lambdaValsAvgRef); % arbs ... find n at the wavelengths corresponding to the beta values we defined above
AkDataAvgRef = interpFxnAl2O3k(lambdaValsAvgRef);
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
hVec = [NaN, TfilmAtop, TfilmM, TfilmAbot, NaN]; % um... film layer thicknesses

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

% spectrum from Caltech
caltechMeasOctAvgFname = 'octoberAvgPlusMinusWithWL.csv'; % First column is lambda (nm), 2 = R, 5=T, 8=A. 3/4=+/-R, 6/7=+/-T, 9/10=+/-A
caltechData = csvread(caltechMeasOctAvgFname,1,0); % read the CSV file into a matrix 

% Save new spectrum over larger range
if writeFiles == 1
    if opticalData==1
        refIndSrc = 'literatureNK';
    elseif opticalData==2
        refIndSrc = 'measuredNKjason';
    elseif opticalData==3
        refIndSrc = 'measuredNKpawan';
    elseif opticalData==4
        refIndSrc = 'measuredNKcorrugatedCaltech';
    elseif opticalData==5
        refIndSrc = 'measuredNKjasonApril2024';
    end
    fNameOut = ['prototypeSimulatedSpectrumFull_',refIndSrc,'.csv'];
    headerString = "WLnm, Rdata, Adata, Tdata\n";
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
plot(caltechData(:,1),caltechData(:,2),'g:') % R
plot(caltechData(:,1),caltechData(:,2)+caltechData(:,3),'g:') % R
plot(caltechData(:,1),caltechData(:,2)+caltechData(:,4),'g:') % R
plot(caltechData(:,1),caltechData(:,8),'r:') % A
plot(caltechData(:,1),caltechData(:,8)+caltechData(:,9),'r:') % A
plot(caltechData(:,1),caltechData(:,8)+caltechData(:,10),'r:') % A
plot(caltechData(:,1),caltechData(:,5),'b:') % T
plot(caltechData(:,1),caltechData(:,5)+caltechData(:,6),'b:') % T
plot(caltechData(:,1),caltechData(:,5)+caltechData(:,7),'b:') % T
hold off
%set(gca, 'YScale', 'log')
xlabel('Wavelength [nm]')
ylabel('R, A, T [#]')
title('R, A, T');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 1.5);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);


