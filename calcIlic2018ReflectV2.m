% Calculate reflectance Ilic 2018
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

% wavelength range for starshot with beta=0.2
lambda0 = 1.2; % um
betaFstarshot = 0.2;
lambdaMaxStarshot = lambda0 * sqrt((1+betaFstarshot)./(1-betaFstarshot));
lambdaValsAvgRef = linspace(lambda0,lambdaMaxStarshot,1000)'; % um

thetaIn = 0; % radians... normal to surface

% larger wavelength range for evaluating max speed in separate code
lambdaVals = linspace(lambda0,3.0,3000)'; % um

% Optical stack Layer thicknesses in microns
hVec = [NaN, 0.104, 0.514, 0.124, 0.492, 0.131, 0.486, 0.131, 0.492, 0.124, 0.514, 0.104, NaN]; % um... film layer thicknesses
%hVec = [NaN, 0.206, NaN]

% Evaluate n and kappa
fNameSiO2 = 'Rodriguez_de_Marcos_2016_SiO2_nk.csv'; % First column is lambda (um), second is n, third is kappa
SiO2nkData = csvread(fNameSiO2,1,0); % read the CSV file into a matrix
%SiO2nkData(:,1) = SiO2nkData(:,1)*(1/1000); % nm --> um
interpFxnSiO2n = griddedInterpolant(SiO2nkData(:,1),SiO2nkData(:,2),'linear','linear'); % create function to describe data and enable interpolations.  
interpFxnSiO2k = griddedInterpolant(SiO2nkData(:,1),SiO2nkData(:,3),'linear','linear');
SiO2nVals = interpFxnSiO2n(lambdaVals); 
SiO2kVals = interpFxnSiO2k(lambdaVals); 
SiO2nValsAvg = interpFxnSiO2n(lambdaValsAvgRef); 
SiO2kValsAvg = interpFxnSiO2k(lambdaValsAvgRef); 

% assemble index of refraction matrix
% BE CAREFUL WHEN ASSEMBLING COMPLEX NUMBERS!  MAKE SURE ALL ARE COLUMN VECTORS!!!
indComplexSiO2 = SiO2nVals + 1i*SiO2kVals; % ... lambda increases down columns 
indComplexGap = ones(size(lambdaVals)) + 1i*zeros(size(lambdaVals)); % gap laters 
nMatIn = [ones(size(lambdaVals)) indComplexSiO2 indComplexGap indComplexSiO2 indComplexGap indComplexSiO2 indComplexGap indComplexSiO2 indComplexGap indComplexSiO2 indComplexGap indComplexSiO2 ones(length(lambdaVals),1)]; % the index of refraction values depend on temperature! In nMat, lambda increases down columns and layer number increases across rows.
%nMatIn = [ones(size(lambdaVals)) indComplexSiO2 ones(size(lambdaVals))];
indComplexSiO2avg = SiO2nValsAvg + 1i*SiO2kValsAvg; % ... lambda increases down columns 
indComplexGapAvg = ones(size(lambdaValsAvgRef)) + 1i*zeros(size(lambdaValsAvgRef)); % gap laters 
nMatInAvg = [ones(size(lambdaValsAvgRef)) indComplexSiO2avg indComplexGapAvg indComplexSiO2avg indComplexGapAvg indComplexSiO2avg indComplexGapAvg indComplexSiO2avg indComplexGapAvg indComplexSiO2avg indComplexGapAvg indComplexSiO2avg ones(length(lambdaValsAvgRef),1)]; % the index of refraction values depend on temperature! In nMat, lambda increases down columns and layer number increases across rows.
%nMatInAvg = [ones(size(lambdaValsAvgRef)) indComplexSiO2avg ones(size(lambdaValsAvgRef))];

% transfer matrix for reflectivity over wide lambda range
pol = 1;    
[Rvecp, Tvecp, Avecp, Mdummy, Vdummy] = FresnelMat(lambdaVals, thetaIn, hVec, nMatIn, pol); 
pol = 0;
[Rvecs, Tvecs, Avecs, Mdummy, Vdummy] = FresnelMat(lambdaVals, thetaIn, hVec, nMatIn, pol); 
Rvec = mean(cat(2,Rvecp,Rvecs),2); % average
Avec = mean(cat(2,Avecp,Avecs),2); % average
Tvec = mean(cat(2,Tvecp,Tvecs),2); % average

% transfer matrix for average reflectivity in beta = 0 - 0.2 range
pol = 1;    
[RvecpAvg, TvecpAvg, AvecpAvg, Mdummy, Vdummy] = FresnelMat(lambdaValsAvgRef, thetaIn, hVec, nMatInAvg, pol); 
pol = 0;
[RvecsAvg, TvecsAvg, AvecsAvg, Mdummy, Vdummy] = FresnelMat(lambdaValsAvgRef, thetaIn, hVec, nMatInAvg, pol); 
RvecAvg = mean(cat(2,RvecpAvg,RvecsAvg),2); % average
AvecAvg = mean(cat(2,AvecpAvg,AvecsAvg),2); % average
TvecAvg = mean(cat(2,TvecpAvg,TvecsAvg),2); % average
Ravg = mean(RvecAvg)
Aavg = mean(AvecAvg)
Tavg = mean(TvecAvg)

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

% plot
figure(2)
hold on
plot(SiO2nkData(:,1)*(1000/1),SiO2nkData(:,2),'c-','LineWidth', 2)
plot(SiO2nkData(:,1)*(1000/1),SiO2nkData(:,3),'m-','LineWidth', 2)
hold off
set(gca, 'YScale', 'log')
xlabel('Wavelength [nm]')
ylabel('n, k')
title('n, k for SiO2');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 1.5);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);

% Save output file for calculation later
% Save new spectrum over larger range
% fNameOut = 'Ilic2018reflectivityDesignA11_usingRodriguez_de_Marcos.csv';
% headerString = "WLnm, Rdata\n";
% h = fopen(fNameOut,'w');
% fprintf(h,headerString);
% for m=1:length(lambdaVals)
%     fprintf(h,'%10.6e, %10.6e\n',lambdaVals(m)*(1000/1),Rvec(m));
% end
% fclose('all');

