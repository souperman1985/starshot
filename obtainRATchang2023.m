% Compare R/A/T for Chang 2023 paper
% See Figure 4b

% Start fresh
clear all
close all
format short
format compact
clc

% wavelengths
betaStarshot = 0.2;
lambda0 = 1300; % nm
lambdaMax = lambda0*sqrt((1+betaStarshot)/(1-betaStarshot)) % nm
%lambda0 = 1309 % ... even if we only look at their experimental range without extrapolation, their absorption is still higher than ours. 
%lambdaMax = 1500
lambdaVals = linspace(lambda0,lambdaMax,1000);

% Input R and T data
fNameRef = 'Chang2023exptRmeas.csv'; % First column is lambda (nm), second is reflectivity (0-1)
refData = csvread(fNameRef,1,0); % read the CSV file into a matrix
fNameRefCalib = 'Chang2023exptRcalib.csv'; % First column is lambda (nm), second is reflectivity (0-1)
refDataCalib = csvread(fNameRefCalib,1,0); % read the CSV file into a matrix
fNameTrans = 'Chang2023exptT.csv'; % First column is lambda (nm), second is transmission (0-1)
transData = csvread(fNameTrans,1,0); % read the CSV file into a matrix
fNameLoss = 'Chang2023exptLoss.csv'; % First column is lambda (nm), second is reflectivity (0-1)
lossData = csvread(fNameLoss,1,0); % read the CSV file into a matrix ... system loss determined using Au mirror in place of film
fNameAuLoss = 'Chang2023exptAuLoss.csv'; % First column is lambda (nm), second is transmission (0-1)
AuLossData = csvread(fNameAuLoss,1,0); % read the CSV file into a matrix ... loss for Au mirror alone (from manufacturer, ripped from paper)

% Interpolate/extrapolate upon common lambda values
interpFxnRef = griddedInterpolant(refData(:,1),refData(:,2),'linear','linear'); % create function to describe data and enable interpolations.  
RefSpectrum = interpFxnRef(lambdaVals); % arbs ... find reflectivity at the wavelengths corresponding to the beta values we defined above
interpFxnRefCalib = griddedInterpolant(refDataCalib(:,1),refDataCalib(:,2),'linear','linear'); % create function to describe data and enable interpolations.  
RefSpectrumCalib = interpFxnRefCalib(lambdaVals); % arbs ... find reflectivity at the wavelengths corresponding to the beta values we defined above
interpFxnTrans = griddedInterpolant(transData(:,1),transData(:,2),'linear','linear'); 
TransSpectrum = interpFxnTrans(lambdaVals); % arbs 
interpFxnLoss = griddedInterpolant(lossData(:,1),lossData(:,2),'linear','nearest');  
LossSpectrumWithAu = interpFxnLoss(lambdaVals); % arbs ...
interpFxnAuLoss = griddedInterpolant(AuLossData(:,1),AuLossData(:,2),'linear','nearest');  
AuLossSpectrum = interpFxnAuLoss(lambdaVals); % arbs 

% seems like the Au loss etc are constant with wavelength.  In this case, just use an average value for all.
%LossSpectrumWithAu = mean(lossData(:,2))*ones(size(lambdaVals)); LossSpectrumWithAu(1)
%AuLossSpectrum = mean(AuLossData(:,2))*ones(size(lambdaVals)); AuLossSpectrum(1)

% Calc system loss including sail film
AbsSpectrumOld = 1-RefSpectrum-TransSpectrum;
AbsSpectrumAlt = 1-RefSpectrumCalib-TransSpectrum;
LossSpectrumWithFilm = 1-RefSpectrum-TransSpectrum; % use the non-calibrated reference spectrum for this, otherwise we're double-calibrating. 

% Account for losses
% 1-TotalLossWithAu = (1-AuLoss)*(1-SystemLoss) --> SystemLoss = 1-[(1-TotalLossWithAu)/(1-AuLoss)]
SystemLossSpect = 1-((1-LossSpectrumWithAu)./(1-AuLossSpectrum)); SystemLossSpect(1)
% 1-LossSpectrumWithFilm = RefSpectrum+TransSpectrum = (1-FilmAbsorptivity)*(1-SystemLoss) --> FilmAbsorptivity = 1-[(1-RefSpectrum-TransSpectrum)/(1-SystemLoss)]
AbsSpectrum = 1-((1-LossSpectrumWithFilm)./(1-SystemLossSpect));

% Calculate average R and T
% The calibrated R spectrum was corrected for system loss by hte authors so we'll use that for the refectivity value
Ravg = mean(RefSpectrum)
RavgCheck = trapz(lambdaVals,RefSpectrum)/(lambdaVals(end)-lambdaVals(1));
RavgCalib = mean(RefSpectrumCalib)
RavgCheckCalib = trapz(lambdaVals,RefSpectrumCalib)/(lambdaVals(end)-lambdaVals(1));
Tavg = mean(TransSpectrum)
TavgCheck = trapz(lambdaVals,TransSpectrum)/(lambdaVals(end)-lambdaVals(1));
AavgOld = mean(AbsSpectrumOld)
AavgAlt = mean(AbsSpectrumAlt)
Aavg = mean(AbsSpectrum)
AavgCheck = trapz(lambdaVals,AbsSpectrum)/(lambdaVals(end)-lambdaVals(1));

% Plot
figure(1)
hold on;
plot(lambdaVals,RefSpectrum,'g-')
plot(lambdaVals,RefSpectrumCalib,'c-')
plot(lambdaVals,TransSpectrum,'b-')
plot(lambdaVals,AbsSpectrum,'r-')
plot(lambdaVals,LossSpectrumWithAu,'k-')
plot(lambdaVals,AuLossSpectrum,'-','color',"#EDB120")
plot(refData(:,1),refData(:,2),'go')
plot(refDataCalib(:,1),refDataCalib(:,2),'co')
plot(transData(:,1),transData(:,2),'bo')
plot(lossData(:,1),lossData(:,2),'ko')
plot(AuLossData(:,1),AuLossData(:,2),'o','color',"#EDB120")
plot(lambdaVals,AbsSpectrumOld,'c:')
plot(lambdaVals,SystemLossSpect,'m:')
hold off;
xlabel('Wavelength [um]')
ylabel('Optical Data Chang 2023 [a.u.]')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 1.5);
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);

RavgCalib+Tavg+Aavg
RavgCalib+Tavg+AavgAlt




% After everything, it appears that the most effective way to calculate this is to use the corrected reflection spectrum and the one provided transmission spectrum.  
% THis is "AavgAlt" in this script. 