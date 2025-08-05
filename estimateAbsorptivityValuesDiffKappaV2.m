% Estimate the absorptivity of different designs in the literature
%
% Copyright (c) 2025 Matthew Campbell
% Permission is hereby granted, free of charge, to use, copy, modify, and distribute this software for any purpose, with or without modification, provided that the above copyright notice and this permission notice appear in all copies.
% THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
%
% absorptivity = (4*pi*kappa/lambda)*thickness*fillFactor

% Start fresh
clear all
close all
format long e
format compact
clc

% Constants
cLight = 299792458; % m/s ... light speed
hPlanck = 6.62607015e-34; % J-s ... Planck's constant
hPlanckBar = 1.054571817e-34; % J-s ... Planck's constant with a bar... h/(2*pi)
kB = 1.380649e-23; % J/K ... Boltzmann constant
eV = 1.602176634e-19; % J per eV
sigmaBB = (2*pi^5*kB^4)/(15*cLight^2*hPlanck^3); % W/m2-K4

% Parameters
betaFinal = 0.2; % for starshot

% Si kappa data
fNameSi = 'Poruba 2000 fig5b data.csv'; % energy [eV] and alpha [cm-1]
absDataSi = csvread(fNameSi,1,0); % read the CSV file into a matrix
lambdaSi = (hPlanck*cLight./(eV*absDataSi(:,1)))*(1e6/1); % (J-s * m/s)/(J/eV --> um
kappaSi = ((absDataSi(:,2).*lambdaSi)/(4*pi))*(100/1)*(1/1e6); % cm-1 * um * cm/m * m/um --> no units
interpFxnKappaSi = griddedInterpolant(flip(lambdaSi),flip(kappaSi),'linear','linear'); % create function to describe data and enable interpolations.  

% SiO2 kappa data
fNameSiO2 = 'Rodriguez_de_Marcos_2016_SiO2_nk.csv'; % First column is lambda (um), second is n, third is kappa
SiO2nkData = csvread(fNameSiO2,1,0); % read the CSV file into a matrix
interpFxnSiO2k = griddedInterpolant(SiO2nkData(:,1),SiO2nkData(:,3),'linear','linear');

% Si3N4 
fNameSN = 'Kischkat2012SiNnk.csv'; % First column is lambda (um), second is n, third is kappa 
SNnkData = csvread(fNameSN,1,0); % read the CSV file into a matrix
SNlambdaData = SNnkData(:,1); % um
SNkappaData = SNnkData(:,3); 
interpFxnSNk = griddedInterpolant(SNlambdaData,SNkappaData,'linear','linear');
% fit power law extrapolation using simple two-point fit
indLambdaA = 1;
lambdaA = SNlambdaData(indLambdaA); % um
% kappaA = SNkappaData(indLambdaA)
lambdaBwant = 2.5; % um
indLambdaB = find(SNlambdaData-lambdaBwant>=0,1);
lambdaB = SNlambdaData(indLambdaB); % um
% kappaB = SNkappaData(indLambdaB)
% b = log(kappaA/kappaB)/log(lambdaA/lambdaB)
% a = kappaA/lambdaA^b
lambdaFitSN = linspace(1,lambdaB,100); % um
% kappaFitSNsimple = a*lambdaFitSN.^b;
% fit power law extrapolation using better fit
pSN = polyfit(log(SNlambdaData(indLambdaA:indLambdaB)),log(SNkappaData(indLambdaA:indLambdaB)),1);
bSN = pSN(1);
aSN = exp(pSN(2));
kappaFitSN = aSN*lambdaFitSN.^bSN;

% TiO2 
fNameTiO2 = 'Kischkat2012TiO2nk.csv'; % First column is lambda (um), second is n, third is kappa 
TiO2nkData = csvread(fNameTiO2,1,0); % read the CSV file into a matrix
TOlambdaData = TiO2nkData(:,1);
TOkappaData = TiO2nkData(:,3);
interpFxnTiO2k = griddedInterpolant(TOlambdaData,TOkappaData,'linear','linear');
% fit power law extrapolation using simple two-point fit
indLambdaA = 1;
lambdaA = TOlambdaData(indLambdaA); % um
lambdaBwant = 2.5; % um
indLambdaB = find(TOlambdaData-lambdaBwant>=0,1);
lambdaB = TOlambdaData(indLambdaB); % um
lambdaFitTO = linspace(1,lambdaB,100); % um
% fit power law extrapolation using better fit
pTO = polyfit(log(TOlambdaData(indLambdaA:indLambdaB)),log(TOkappaData(indLambdaA:indLambdaB)),1);
bTO = pTO(1);
aTO = exp(pTO(2));
kappaFitTO = aTO*lambdaFitTO.^bTO;

% Al2O3
fNameK2012 = 'Kischkat2012_Al2O3nk.csv'; 
tempMatrixK2012 = csvread(fNameK2012,1,0); % read the CSV file into a matrix
lambdaValsK2012 = tempMatrixK2012(1:end,1); % um ... first column is wavelength
spectrumK2012k = tempMatrixK2012(1:end,3); % k
interpFxnAl2O3k = griddedInterpolant(lambdaValsK2012,spectrumK2012k,'linear','linear'); % create function to describe data and enable interpolations.  
% fit power law extrapolation using simple two-point fit
indLambdaA = 1;
lambdaA = lambdaValsK2012(indLambdaA); % um
lambdaBwant = 2.5; % um
indLambdaB = find(lambdaValsK2012-lambdaBwant>=0,1);
lambdaB = lambdaValsK2012(indLambdaB); % um
lambdaFitAO = linspace(1,lambdaB,100); % um
% fit power law extrapolation using better fit
pAO = polyfit(log(lambdaValsK2012(indLambdaA:indLambdaB)),log(spectrumK2012k(indLambdaA:indLambdaB)),1);
bAO = pAO(1);
aAO = exp(pAO(2));
kappaFitAO = aAO*lambdaFitAO.^bAO;

% MoS2 kappa data
fNameMunk = 'Munkhbat2022MoS2nk.csv'; %fprintf('Using Single Layer MoS2 Data\n');
tempMatrixMunk = csvread(fNameMunk,2,0); % read the CSV file into a matrix
lambdaValsMunk = tempMatrixMunk(:,1); % um ... first column is wavelength
spectrumMunkK = tempMatrixMunk(:,3); % arbitrary units, (k) extinction coefficient data @ 300 K
interpFxnMoS2k = griddedInterpolant(lambdaValsMunk,spectrumMunkK,'linear','linear'); % create function to describe data and enable interpolations.  

% 

figure(100)
hold on
plot(lambdaSi,kappaSi,'k-')
plot(SiO2nkData(:,1),SiO2nkData(:,3),'g-')
plot(SNlambdaData,SNkappaData,'c-')
plot(TOlambdaData,TOkappaData,'m-')
plot(lambdaValsK2012,spectrumK2012k,'b-')
plot(lambdaValsMunk,spectrumMunkK,'r-')
%plot(lambdaFitSN,kappaFitSNsimple,'k:')
plot(lambdaFitSN,kappaFitSN,'k.')
plot(lambdaFitTO,kappaFitTO,'k.')
plot(lambdaFitAO,kappaFitAO,'k.')
hold off
ax = axis;
axis([0, 3, ax(3), ax(4)]);
set(gca,'yscale','log')
xlabel('Wavelength [um]')
ylabel('Kappa [arbs]')
title('Extinction coefficient data')
legend('Si','SiO2','Si3N4','TiO2','Al2O3','MoS2')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);

% Justificaiton for using de Marcos for SiO2... from Jason
% The Rodriguez-de Marcos one is probably best (https://doi.org/10.1364/OME.6.003622). 
% It was deposited using E-beam evaporation, and the films were 54 and 125 nm thick. 
% Kischkat (https://doi.org/10.1364/OME.6.003622) is 510 nm thick, used reactive ion sputtering, and didn't directly measure the laser wavelengths. 
% Franta (https://doi.org/10.1117/12.2227580) had 878 nm thick films, used E-beam evaporation, and the SiO2 was polycrystalline with water/organic contaminants adsorbed in the pores.


% Fabricated prototype
disp('---------------')
% wavelengths
lambda0 = 1.2; % um
lambdaFinal = lambda0 * sqrt((1+betaFinal)/(1-betaFinal)) % um
lambdaValsFabProto = linspace(lambda0,lambdaFinal,1000);
% kappa values
kappaValsAl2O3FabProto = aAO*lambdaValsFabProto.^bAO;
kappaValsMoS2FabProto = interpFxnMoS2k(lambdaValsFabProto);
% thicknesses
tA1 = 21; % nm ... ellipsometry, bottom alumina thickness
tM = 53; % nm ... fit to optical data, MoS2 thickness
tA2 = 51; % nm ... ellipsometry, top alumina thickness
% fill factors
ff1 = 1;
ff2 = 1;
ff3 = 1;
% absorptivity
alphaVecFabProto = (4*pi./lambdaValsFabProto).*(kappaValsAl2O3FabProto*tA1*ff1 + kappaValsMoS2FabProto*tM*ff2 + kappaValsAl2O3FabProto*tA2*ff3)*(1/1000);
alphaFabProto = mean(alphaVecFabProto)
%alphaVecFabProto2 = 1 - ((1-(4*pi./lambdaValsFabProto).*kappaValsAl2O3FabProto*tA1*ff1*(1/1000)).*(1-(4*pi./lambdaValsFabProto).*kappaValsMoS2FabProto*tM*ff2*(1/1000)).*(1-(4*pi./lambdaValsFabProto).*kappaValsAl2O3FabProto*tA2*ff3*(1/1000)));
%alphaFabProto2 = mean(alphaVecFabProto2)

%keyboard

% plot
figure(1)
hold on;
plot(lambdaValsFabProto,kappaValsAl2O3FabProto);
plot(lambdaValsFabProto,kappaValsMoS2FabProto)
hold off;
set(gca,'yscale','log')
xlabel('Wavelength [um]')
ylabel('Kappa [arbs]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);


% Optimized proposed design
disp('---------------')
% wavelengths
lambda0 = 1.2; % um
lambdaFinal = lambda0 * sqrt((1+betaFinal)/(1-betaFinal)); % um
lambdaValsOptDesign = linspace(lambda0,lambdaFinal,1000);
% kappa values
kappaValsAl2O3OptDesign = aAO*lambdaValsOptDesign.^bAO;
kappaValsMoS2OptDesign = interpFxnMoS2k(lambdaValsOptDesign);
% thicknesses
%tA1 = 8.9465; % nm ... Dec 2024 optimization
%tM = 74.6300; % nm 
%tA2 = 8.9465; % nm 
tA1 = 18.5485; % nm ... March 2025 optimization
tM = 62.6900;
tA1 = 18.5485;
% fill factors
ff1 = 1;
ff2 = 1;
ff3 = 1;
% absorptivity
alphaVecOptDesign = (4*pi./lambdaValsOptDesign).*(kappaValsAl2O3OptDesign*tA1*ff1 + kappaValsMoS2OptDesign*tM*ff2 + kappaValsAl2O3OptDesign*tA2*ff3)*(1/1000);
alphaOptDesign = mean(alphaVecOptDesign)
% plot
figure(2)
hold on;
plot(lambdaValsOptDesign,kappaValsAl2O3OptDesign);
plot(lambdaValsOptDesign,kappaValsMoS2OptDesign)
hold off;
set(gca,'yscale','log')
xlabel('Wavelength [um]')
ylabel('Kappa [arbs]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);


% Brewer 2022
disp('---------------')
% wavelengths
lambda0 = 1.2; % um
lambdaFinal = lambda0 * sqrt((1+betaFinal)/(1-betaFinal)); % um
lambdaValsBrewer = linspace(lambda0,lambdaFinal,1000);
% kappa
kappaValsSiN3Brewer = aSN*lambdaValsBrewer.^bSN;
kappaValsMoS2Brewer = interpFxnMoS2k(lambdaValsBrewer);
% dimensions
periodY = 1.16; % um
periodX = periodY * sqrt(3); % um
dHole = 0.9*periodY; % um
tSN1 = 5; % nm ... Si3N4 thickness on top,  see Fig 2
tM = 90; % nm ... MoS2 thickness
tSN2 = 5; % nm ... Si3N4 thickness on bottom
% fill factor
Amaterial = (periodY * periodX - 2*pi*(dHole/2)^2) % um^2 ... each unit cell contains 2 holes
Auc = periodY * periodX; % um2
ff1 = Amaterial/Auc
ff2 = ff1;
ff3 = ff1;
% absorptivity
alphaVecBrewer2022 = (4*pi./lambdaValsBrewer).*(kappaValsSiN3Brewer*tSN1*ff1 + kappaValsMoS2Brewer*tM*ff2 + kappaValsSiN3Brewer*tSN2*ff3)*(1/1000);
alphaBrewer2022 = mean(alphaVecBrewer2022)
% plot
figure(3)
hold on;
plot(lambdaValsBrewer,kappaValsSiN3Brewer);
plot(lambdaValsBrewer,kappaValsMoS2Brewer)
hold off;
set(gca,'yscale','log')
xlabel('Wavelength [um]')
ylabel('Kappa [arbs]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);


% Chang 2024... design from Fig 3
disp('---------------')
% wavelengths
lambda0 = 1.3; % um ... laser wavelength at zero sail velocity. 
lambdaFinal = lambda0 * sqrt((1+betaFinal)/(1-betaFinal)); % um
lambdaValsChang = linspace(lambda0,lambdaFinal,1000);
% kappa
kappaValsSiN3Chang = aSN*lambdaValsChang.^bSN;
kappaValsSiChang = interpFxnKappaSi(lambdaValsChang);
% dimensions
tSN = 400; % nm
tSi = 321; % nm
aPattern = 1440e-9; % m ... measured from SEM in Figure 3
dHole = 2*500e-9; % m ... confirmed from Fig 3 SEM
Amaterial = aPattern^2 - pi*(dHole/2)^2;
Auc = aPattern^2; % m2
ff1 = Amaterial/Auc % holes in Si3N4 only
ff2 = 1; % Si layer is unpatterned 
% absorptivity
alphaVecChang2024 = (4*pi./lambdaValsChang).*(kappaValsSiN3Chang*tSN*ff1 + kappaValsSiChang*tSi*ff2)*(1/1000);
alphaChang2024 = mean(alphaVecChang2024)
% plot
figure(4)
hold on;
plot(lambdaValsChang,kappaValsSiN3Chang);
plot(lambdaValsChang,kappaValsSiChang)
hold off;
set(gca,'yscale','log')
xlabel('Wavelength [um]')
ylabel('Kappa [arbs]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);


% Ilic 2018... design A11 (see sup inf)
disp('---------------')
% wavelengths
lambda0 = 1.2; % um ... laser wavelength at zero sail velocity. 
lambdaFinal = lambda0 * sqrt((1+betaFinal)/(1-betaFinal)); % um
lambdaValsIlic = linspace(lambda0,lambdaFinal,1000);
% kappa
kappaValsSiO2Ilic = interpFxnSiO2k(lambdaValsIlic);
% dimensions
tSiO2total = sum([0.104, 0.124, 0.131, 0.131, 0.124, 0.104])*(1000/1); % um --> nm
ff = 1;
% absorptivity
alphaMaxVecIlic2018 = (4*pi./lambdaValsIlic).*(kappaValsSiO2Ilic*tSiO2total*ff)*(1/1000);
alphaMaxIlic2018 = mean(alphaMaxVecIlic2018)
% plot
figure(5)
hold on;
plot(lambdaValsIlic,kappaValsSiO2Ilic);
hold off;
set(gca,'yscale','log')
xlabel('Wavelength [um]')
ylabel('Kappa [arbs]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);

% Lien 2022
disp('---------------')
% wavelengths
lambda0 = 1.064; % um ... laser wavelength at zero sail velocity. 
lambdaFinal = lambda0 * sqrt((1+betaFinal)/(1-betaFinal)); % um
lambdaValsLien = linspace(lambda0,lambdaFinal,1000);
% kappa
kappaValsSi3N4Lien = aSN*lambdaValsLien.^bSN;
% dimensions
tSN = 690; % nm
aPattern = 1064e-9; % m
dHole = 415e-9; % m
Amaterial = aPattern^2 - pi*(dHole/2)^2; % m2
Auc = aPattern^2; % m2
ff1 = Amaterial/Auc % 
% absorptivity
alphaVecLien2022 = (4*pi./lambdaValsLien).*(kappaValsSi3N4Lien*tSN*ff1)*(1/1000);
alphaLien2022 = mean(alphaVecLien2022)
% plot
figure(6)
hold on;
plot(lambdaValsLien,kappaValsSi3N4Lien);
hold off;
set(gca,'yscale','log')
xlabel('Wavelength [um]')
ylabel('Kappa [arbs]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);


% Salary 2019
disp('---------------')
% wavelengths
lambda0 = 1.3; % um 
lambdaFinal = lambda0 * sqrt((1+betaFinal)/(1-betaFinal)); % um
lambdaValsSalary = linspace(lambda0,lambdaFinal,1000);
% kappa
kappaValsSiO2Salary = interpFxnSiO2k(lambdaValsSalary);
kappaValsSiSalary = interpFxnKappaSi(lambdaValsSalary);
% dimensions
tSiO2 = 60; % nm
tSi = 89; % nm
hDisk = 323; % nm
delta = 639e-9; % m
dCircle = 300e-9; % m ... range is 180 to 450 nm, and I pulled optical properties at 300 nm.
Amaterial = pi*(dCircle/2)^2; % m2
Auc = delta^2; % m2
ff1 = Amaterial/Auc;
ff2 = 1;
ff3 = 1;
% absorptivity
alphaVecSalary2019 = (4*pi./lambdaValsSalary).*(kappaValsSiSalary*hDisk*ff1 + kappaValsSiSalary*tSi*ff2 + kappaValsSiO2Salary*tSiO2*ff3)*(1/1000);
alphaSalary2019 = mean(alphaVecSalary2019)
% plot
figure(7)
hold on;
plot(lambdaValsSalary,kappaValsSiO2Salary);
plot(lambdaValsSalary,kappaValsSiSalary);
hold off;
set(gca,'yscale','log')
xlabel('Wavelength [um]')
ylabel('Kappa [arbs]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);


% Santi 2022 ... design with TiO2 + SiO2 + TiO2
disp('---------------')
% wavelengths
lambda0 = 1.064; % um 
lambdaFinal = lambda0 * sqrt((1+betaFinal)/(1-betaFinal)); % um
lambdaValsSanti = linspace(lambda0,lambdaFinal,1000);
% kappa
kappaValsTOSanti = aTO*lambdaValsSanti.^bTO;
kappaValsSiO2Santi = interpFxnSiO2k(lambdaValsSanti);
% dimensions
tTiO2 = 121; % nm
tSiO2 = 203; % nm
ff1 = 1;
ff2 = 1;
ff3 = 1;
% absorptivity
alphaVecSanti2022 = (4*pi./lambdaValsSanti).*(kappaValsTOSanti*tTiO2*ff1 + kappaValsSiO2Santi*tSiO2*ff2 + kappaValsTOSanti*tTiO2*ff3)*(1/1000);
alphaSanti2022 = mean(alphaVecSanti2022)
% plot
figure(8)
hold on;
plot(lambdaValsSanti,kappaValsTOSanti);
plot(lambdaValsSanti,kappaValsSiO2Santi);
hold off;
set(gca,'yscale','log')
xlabel('Wavelength [um]')
ylabel('Kappa [arbs]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);


% Taghavi 2022
disp('---------------')
% wavelengths
lambda0 = 1.3; % um 
lambdaFinal = lambda0 * sqrt((1+betaFinal)/(1-betaFinal)); % um
lambdaValsTaghavi = linspace(lambda0,lambdaFinal,1000);
% kappa
kappaValsSiO2Taghavi = interpFxnSiO2k(lambdaValsTaghavi);
kappaValsSiTaghavi = interpFxnKappaSi(lambdaValsTaghavi);
% dimensions
tSiO2 = 60; % m
tSi = 103; % m
hDisk = 306; % nm ...hR
delta = 653e-9; % m ... deltaR
dCircle = 300e-9; % m ... range is 180 to 450 nm, and I pulled optical properties at 300 nm.
Amaterial = pi*(dCircle/2)^2; % m2
Auc = delta^2; % m2
ff1 = Amaterial/Auc;
ff2 = 1;
ff3 = 1;
% absorptivity
alphaVecTaghavi2022 = (4*pi./lambdaValsTaghavi).*(kappaValsSiTaghavi*hDisk*ff1 + kappaValsSiTaghavi*tSi*ff2 + kappaValsSiO2Taghavi*tSiO2*ff3)*(1/1000);
alphaTaghavi2022 = mean(alphaVecTaghavi2022)
% plot
figure(9)
hold on;
plot(lambdaValsTaghavi,kappaValsSiO2Taghavi);
plot(lambdaValsTaghavi,kappaValsSiTaghavi);
hold off;
set(gca,'yscale','log')
xlabel('Wavelength [um]')
ylabel('Kappa [arbs]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);


    
    





