% Calculate R-A-T of laminate film

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
    fNameOut = 'prototypeSimulatedSpectrumFull_measuredNKjasonApril2024.csv';
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








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate reflectivity as a function of angle through the acceleration trajectory (for the data range of Caltech)
lambdaValsContour = linspace(1.0,1.6,300)'; % um
thetaValsContour = linspace(0,80,101)*(pi/180); % deg-->rad

% calculate optical constants
AnDataContour = interpFxnAl2O3n(lambdaValsContour); % arbs ... find n at the wavelengths corresponding to the beta values we defined above
AkDataContour = interpFxnAl2O3k(lambdaValsContour);
MnDataContour = interpFxnMoS2n(lambdaValsContour); 
MkDataContour = interpFxnMoS2k(lambdaValsContour);

% Assemble complex index of refraction
% BE CAREFUL WHEN ASSEMBLING COMPLEX NUMBERS!  MAKE SURE ALL ARE COLUMN VECTORS!!!
indComplexInAcontour = AnDataContour + 1i*AkDataContour; % face layers ... here lambda increases down columns 
indComplexInMcontour = MnDataContour + 1i*MkDataContour; % middle later ... here lambda increases down columns
nMatInContour = [ones(length(lambdaValsContour),1)  indComplexInAcontour  indComplexInMcontour  indComplexInAcontour  ones(length(lambdaValsContour),1)]; % the index of refraction values depend on temperature! In nMat, lambda increases down columns and layer number increases across rows.

% Calculate R/T/A for larger wavelength range at all angles to find emissivity
pol = 1;
[Rmatp, ~, Amatp, ~, ~] = FresnelMat(lambdaValsContour, thetaValsContour, hVec, nMatInContour, pol);
pol = 0;
[Rmats, ~, Amats, ~, ~] = FresnelMat(lambdaValsContour, thetaValsContour, hVec, nMatInContour, pol);
Rmat = mean(cat(3,Rmatp,Rmats),3); % average polarizations
Amat = mean(cat(3,Amatp,Amats),3); % average polarizations



% find the maximum deviation from 0-degree reflectivity up to 23.6 degrees 
angleNA04 = 23.6;
indWantAngle = find(thetaValsContour-angleNA04*(pi/180)>=0,1);
angleHave = thetaValsContour(indWantAngle)*(180/pi)
indWantLambdaLow = find(lambdaValsContour-lambda0>=0,1);
lambdaLowHave = lambdaValsContour(indWantLambdaLow)
indWantLambdaHigh = find(lambdaValsContour-lambdaMaxStarshot>=0,1);
lambdaHighHave = lambdaValsContour(indWantLambdaHigh)

RnormalEverywhereMat = repmat(Rmat(:,1),1,length(thetaValsContour)); 
AnormalEverywhereMat = repmat(Amat(:,1),1,length(thetaValsContour)); 
RpctDiffFromNormalMat = abs(Rmat-RnormalEverywhereMat)./Rmat;
ApctDiffFromNormalMat = abs(Amat-AnormalEverywhereMat)./Amat;

RmaxPctDiff = max(RpctDiffFromNormalMat(indWantLambdaLow:indWantLambdaHigh,1:indWantAngle),[],2);
max(RmaxPctDiff)
AmaxPctDiff = max(ApctDiffFromNormalMat(indWantLambdaLow:indWantLambdaHigh,1:indWantAngle),[],2);
max(AmaxPctDiff)

figure(100); 
hold on;
plot(lambdaValsContour(indWantLambdaLow:indWantLambdaHigh),RmaxPctDiff*100)
plot(lambdaValsContour(indWantLambdaLow:indWantLambdaHigh),AmaxPctDiff*100)
hold off;
xlabel('Wavelength [um]')
ylabel('PctDiff')

keyboard



% plot it
figure(20)
levels = 1000;
levelVals = [1.1, 1.25, 1.5, 2,0, 3.0];
%levelVals = 10.^(-2:1:10);
hold on;
contourf(lambdaValsContour,thetaValsContour*(180/pi),Amat',levels,'LineStyle','none')
[C,h] = contour(lambdaValsContour,thetaValsContour*(180/pi),Amat',levelVals,'LineStyle','-','LineColor','black'); 
%plot(DWratioPrototype, HTratioPrototype,'wo','MarkerSize',10)
hold off;
xlabel('Wavelength [um]')
ylabel('Incident angle [deg]')
% xlabel('D [\mum]')
% ylabel('H [\mum]')
title('Absorptivity of Fabricated Prototype')
colorbar;
clabel(C,h,'FontSize',15)
h.LineWidth = 2;
h.LineStyle = ':';
figHand = gca;
figHand.XAxis.FontSize = 18;
figHand.YAxis.FontSize = 18;
figHand.Title.FontSize = 20;
figHand.FontWeight = 'bold';
figHand.FontName = 'Calibri';
figHand.LineWidth = 1.5; % tick and border thickness

% plot it
figure(21)
levels = 1000;
levelVals = [1.1, 1.25, 1.5, 2,0, 3.0];
%levelVals = 10.^(-2:1:10);
hold on;
contourf(lambdaValsContour,thetaValsContour*(180/pi),Rmat',levels,'LineStyle','none')
[C,h] = contour(lambdaValsContour,thetaValsContour*(180/pi),Rmat',levelVals,'LineStyle','-','LineColor','black'); 
%plot(DWratioPrototype, HTratioPrototype,'wo','MarkerSize',10)
hold off;
xlabel('Wavelength [um]')
ylabel('Incident angle [deg]')
% xlabel('D [\mum]')
% ylabel('H [\mum]')
title('Reflectivity of Fabricated Prototype')
colorbar;
clabel(C,h,'FontSize',15)
h.LineWidth = 2;
h.LineStyle = ':';
figHand = gca;
figHand.XAxis.FontSize = 18;
figHand.YAxis.FontSize = 18;
figHand.Title.FontSize = 20;
figHand.FontWeight = 'bold';
figHand.FontName = 'Calibri';
figHand.LineWidth = 1.5; % tick and border thickness



outputContour = 1;
if outputContour
    fNameOut = append('fabricatedPrototypeReflectivityContour.csv'); % g/m2
    RmatT = Rmat';
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',lambdaValsContour); fprintf(h,'\n'); % header line = x-axis = wavelength [um].  Note the initial blank (leading comma) 
    for n=length(thetaValsContour):-1:1 % y-axis AND data.  y-axis = angle from normal axis [rad-->deg]
        fprintf(h,'%7.5e',thetaValsContour(n)*(180/pi)); fprintf(h,' , %7.5e',RmatT(n,:)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = append('fabricatedPrototypeAbsorptivityContour.csv'); % g/m2
    h = fopen(fNameOut,'w');
    AmatT = Amat';
    fprintf(h, ',%7.5e',lambdaValsContour); fprintf(h,'\n'); % header line = x-axis = wavelength [um].  Note the initial blank (leading comma) 
    for n=length(thetaValsContour):-1:1 % y-axis AND data.  y-axis = angle from normal axis [rad-->deg]
        fprintf(h,'%7.5e',thetaValsContour(n)*(180/pi)); fprintf(h,' , %7.5e',AmatT(n,:)); fprintf(h,'\n'); % next series of lines
    end
end