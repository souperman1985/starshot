% Script to show bending properties of hexagons 
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

trenches = 1; % 1 --> trenches... 0 --> ribs.  Matters for areal density! 

saveContourMats = 0;

% Constants
rhoA = 3200; % kg/m3 ... Al2O3
youngA = 170e9; % Pa 
poissonA = 0.24; 
rhoM = 5060; % kg/m3 ... value for bulk MoS2
youngM = 15e9; % Pa ... Graczykowski 2017, for polycrystalline MoS2.  
poissonM = 0.25; % ... poisson (in the middle of Peng 2013, Cooper 2013, Woo 2016, Graczykowski 2017)

model = 2;

% Parameters for prototype
if model == 1
    hexD = 77; % um 
    hexW = 15; % um 
    hexH = 10; % um 
    tA1 = 21e-3; % um 
    tM = 53e-3; % um 
    tA2 = 51e-3; % um 
    nameOut = "Prototype";
end

% Parameters for optimal film
if model == 2
    hexD = 70; % um
    hexW = 2; % um
    hexH = 3; % um
    tA1 = 18.5485e-3; % um 
    tM = 62.6900e-3;
    tA2 = 18.5485e-3;
    nameOut = "Improved";
end

% Calculate properties of composite: Young and Poisson
% See Liu 2009 Eqn 29 and Raju 2018 Eqn 6
hexT = tA1 + tM + tA2; % um
fracM = tM/hexT; % fractional thickness of MoS2
fracA = 1-fracM;
poissonComp = fracA*poissonA + fracM*poissonM % composite Poisson
youngComp =(fracA*youngA + fracM*youngM) + (fracA*fracM*youngA*youngM*(poissonA-poissonM)^2)/(fracA*youngA*(1-poissonM^2)+fracM*youngM*(1-poissonA^2))

DWratioPrototype = hexD/hexW
HTratioPrototype = hexH/hexT

% vectors for plotting
hexDbyhexWvec = linspace(1,50,201);
hexHbyhexTvec = linspace(1,100,202);
[hexDbyhexWmat, hexHbyhexTmat] = meshgrid(hexDbyhexWvec,hexHbyhexTvec);
hexDmat = hexDbyhexWmat*hexW; % um
hexHmat = hexHbyhexTmat*hexT; % um
hexAmat = hexDmat/sqrt(3); % um

% bending stiffness 
AunitCellmat = (sqrt(3)/2)*(hexDmat+hexW).^2; % um2...  Includes base hexagon and half of rib width
Acorn = (sqrt(3)/2)*hexW^2; % um2...  six diamond-shaped areas in corners of hexagon on tops of ribs
AmainMat = AunitCellmat-Acorn; % um2 ... the rest of the cell
Kcorn = (youngComp*hexT^3)/(12*(1-poissonComp^2)); % Pa*um3
KmainMat = youngComp/(12*(1-poissonComp^2))*(hexT^3+3*hexT*hexHmat.^2+(2*hexT*hexHmat.^3)./(hexAmat+hexW)); % Pa*um3
KeffMat = AunitCellmat./((AmainMat./KmainMat)+(Acorn/Kcorn)); % Pa*um3
BsefMat = KeffMat/Kcorn; % bending stiffness enhancement factor
simpleBsefMat = (hexDmat/hexW+1).^2; % for large hexH

% single point
hexA = hexD/sqrt(3); % um
AunitCell = (sqrt(3)/2)*(hexD+hexW).^2; % um2...  Includes base hexagon and half of rib width
Acorn = (sqrt(3)/2)*hexW^2; % um2...  six diamond-shaped areas in corners of hexagon on tops of ribs
Amain = AunitCell-Acorn; % um2 ... the rest of the cell
Kcorn = (youngComp*hexT^3)/(12*(1-poissonComp^2)); % Pa*um3
Kmain = youngComp/(12*(1-poissonComp^2))*(hexT^3+3*hexT*hexH.^2+(2*hexT*hexH.^3)./(hexA+hexW)); % Pa*um3
Keff = AunitCell./((Amain./Kmain)+(Acorn/Kcorn)); % Pa*um3
Bsef = Keff/Kcorn % bending stiffness enhancement factor
simpleBsef = (hexD/hexW+1).^2; % for large hexH


% tensile stiffness 
alpha = 0.234; 
beta = 0.062; 
%!!!AchannelMat = hexT*(2*(hexHmat-hexT)+hexW); % um2 ... cross sectional area of U-section
AchannelMat = hexT*(2*hexHmat+hexW); % um2 ... cross sectional area of U-section
IchannelMat = (hexT*hexW^3/12)+(hexHmat*hexT^3/6)+(hexHmat*hexT*hexW^2/2);% um4 ... moment of inertia for channels... I'm not sure how this was derived though. 
%!!!IchannelMat = (hexHmat*hexW^3)/12 - ((hexHmat-hexT)*(hexW-2*hexT)^3)/12;% um4 ... moment of inertia for channels through centroid and line parallel to rib sides.  U = rectangle - void 
LambdaMmat = AchannelMat.*(hexDmat.^2 + 3*(2.4+1.5*poissonComp)*hexW^2); % um4
poissonMeffMat = (LambdaMmat-36*IchannelMat)./(LambdaMmat+108*IchannelMat); % no units ... poisson for longitudinal tension
PsiAmat = hexT^2*hexHmat.^4*hexW^4; % um^10
%PsiAmat = ((hexT*hexHmat.^2*hexW^2/4)+(IchannelMat/2).*((2*hexHmat.^2+hexT*hexW)./(2*hexHmat+hexW))).^2; % um10
youngMeffMat = (youngComp./(1-poissonMeffMat.^2)).*(((alpha*IchannelMat)./(hexT*hexDmat.^3))+((beta*PsiAmat)./(hexT^4*hexDmat.^5*(2*hexT+hexW)))); % Pa ... Young's modulus for longitudinal tension
TsrfMat = (hexHmat/hexT).*(youngMeffMat/youngComp); % tensile stiffnes reduction factor
SplanM = (youngComp*hexT)/(1-poissonComp^2); % Pa*um ... bending stiffness for planar film (non-corrugated)
ScorrMmat = TsrfMat*SplanM; % Pa*um ... bending stiffness for corrugated film

% single
%!!!Achannel = hexT*(2*(hexH-hexT)+hexW); % um2 ... cross sectional area of U-section
Achannel = hexT*(2*hexH+hexW); % um2 ... cross sectional area of U-section
Ichannel = (hexT*hexW^3/12)+(hexH*hexT^3/6)+(hexH*hexT*hexW^2/2);% um4 ... moment of inertia for channels... I'm not sure how this was derived though. 
%!!!Ichannel = (hexH*hexW^3)/12 - ((hexH-hexT)*(hexW-2*hexT)^3)/12;% um4 ... moment of inertia for channels through centroid and line parallel to rib sides.  U = rectangle - void 
LambdaM = Achannel.*(hexD.^2 + 3*(2.4+1.5*poissonComp)*hexW^2); % um4
poissonMeff = (LambdaM-36*Ichannel)./(LambdaM+108*Ichannel); % no units ... poisson for longitudinal tension
PsiA = hexT^2*hexH.^4*hexW^4; % um^10
%PsiAmat = ((hexT*hexHmat.^2*hexW^2/4)+(IchannelMat/2).*((2*hexHmat.^2+hexT*hexW)./(2*hexHmat+hexW))).^2; % um10
%youngMeff = (youngM./(1-poissonMeff.^2)).*(((alpha*Ichannel)./(hexT*hexD.^3))+((beta*PsiA)./(hexT^4*hexD.^5*(2*hexT+hexW)))); % Pa ... Young's modulus for longitudinal tension
% found an old bug: looks like I had 2*t+w in last part before, rather than 2*h+w. 
youngMeff = (youngComp./(1-poissonMeff.^2)).*(((alpha*Ichannel)./(hexT*hexD.^3))+((beta*PsiA)./(hexT^4*hexD.^5*(2*hexH+hexW)))); % Pa ... Young's modulus for longitudinal tension
Tsrf = (hexH/hexT).*(youngMeff/youngComp) % tensile stiffnes reduction factor
SplanM = (youngComp*hexT)/(1-poissonComp^2); % Pa*um ... bending stiffness for planar film (non-corrugated)
ScorrM = Tsrf*SplanM; % Pa*um ... bending stiffness for corrugated film

% Equivalent stiffness and thickness for a flat plate
% K = (E*t^3)/(12*(1-nu^2)) ... bending stiffness of flat plate
% S = (E*t)/(1-nu^2) ... tensile stiffness of flat plate
% assume nu is the flat plate value rather than the equivalent value for the corrugated plate
youngCompEq = ((ScorrM.^3*(1-poissonComp^2)^2)./(12*Keff)).^(1/2) % Pa
hexTeq = (12*Keff./ScorrM).^(1/2) % um



% Non-corrugated areal density
rhoAnc = (rhoA*tA1 + rhoM*tM + rhoA*tA2)*(1000/1)*(1/1e6) % g/m2

% Corrugated areal density
if trenches == 1
    massMat = ...
        (rhoA * (tA1*sqrt(3)/2*(hexDmat+hexW).^2 + hexHmat*sqrt(3)/2.*((hexDmat+2*tA1).^2-hexDmat.^2)) + ...
        rhoM * (tM*sqrt(3)/2*(hexDmat+hexW).^2 + hexHmat*sqrt(3)/2.*((hexDmat+2*tA1+2*tM).^2-(hexDmat+2*tA1).^2)) + ...
        rhoA * (tA2*sqrt(3)/2*(hexDmat+hexW).^2 + hexHmat*sqrt(3)/2.*((hexDmat+2*tA1+2*tM+2*tA2).^2-(hexDmat+2*tA1+2*tM).^2)))*(1000/1)*(1/1e6)^3; % g
else % ribs
    massMat = ...
        (rhoA * (tA1*sqrt(3)/2*(hexDmat+hexW).^2 + hexHmat*sqrt(3)/2.*(hexDmat.^2-(hexDmat-2*tA1).^2)) + ...
        rhoM * (tM*sqrt(3)/2*(hexDmat+hexW).^2 + hexHmat*sqrt(3)/2.*((hexDmat-2*tA1).^2-(hexDmat-2*tA1-2*tM).^2)) + ...
        rhoA * (tA2*sqrt(3)/2*(hexDmat+hexW).^2 + hexHmat*sqrt(3)/2.*((hexDmat-2*tA1-2*tM).^2-(hexDmat-2*tA1-2*tM-2*tA2).^2)))*(1000/1)*(1/1e6)^3; % g
end
AplaneMat = sqrt(3)/2 * (hexDmat+hexW).^2 *(1/1e6)^2; % m2
rhoAcMat = massMat./AplaneMat; % g/m2
ADratio = rhoAcMat/rhoAnc;

% single
if trenches == 1
    mass = ...
        (rhoA * (tA1*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1).^2-hexD.^2)) + ...
        rhoM * (tM*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1+2*tM).^2-(hexD+2*tA1).^2)) + ...
        rhoA * (tA2*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1+2*tM+2*tA2).^2-(hexD+2*tA1+2*tM).^2)))*(1000/1)*(1/1e6)^3; % g    
else % ribs
    mass = ...
        (rhoA * (tA1*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*(hexD.^2-(hexD-2*tA1).^2)) + ...
        rhoM * (tM*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD-2*tA1).^2-(hexD-2*tA1-2*tM).^2)) + ...
        rhoA * (tA2*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD-2*tA1-2*tM).^2-(hexD-2*tA1-2*tM-2*tA2).^2)))*(1000/1)*(1/1e6)^3; % g
end
Aplane = sqrt(3)/2 * (hexD+hexW).^2 *(1/1e6)^2; % m2
rhoAc = mass./Aplane % g/m2
ADratio1 = rhoAc/rhoAnc


% Plot

figure(1)
numLevels = 1000;
levelVals = [100, 300, 1000, 3000, 5000, 8000];
%levelVals = 10.^(-2:1:7);
hold on;
contourf(hexDbyhexWvec,hexHbyhexTvec,BsefMat,numLevels,'LineStyle','none')
[C,h] = contour(hexDbyhexWvec,hexHbyhexTvec,BsefMat,levelVals,'LineStyle','-','LineColor','black'); 
plot(DWratioPrototype, HTratioPrototype,'wo','MarkerSize',10)
hold off;
xlabel('D/W')
ylabel('H/T')
% xlabel('D [\mum]')
% ylabel('H [\mum]')
title('Bending Stiffness Enhancement Factor')
clabel(C,h,'FontSize',15)
h.LineWidth = 2;
h.LineStyle = ':';
figHand = gca;
figHand.XAxis.FontSize = 18;
figHand.YAxis.FontSize = 18;
figHand.Title.FontSize = 20;
figHand.FontWeight = 'bold';
figHand.FontName = 'Calibri';
figHand.LineWidth = 1.5;

figure(2)
levels = 1000;
levelVals = [0.01 1 10 100 1000 10000];
%levelVals = 10.^(-2:1:10);
hold on;
%contourf(hexDbyhexWmat,hexHbyhexTmat,log10(TsrfAmat),levels,'LineStyle','none')
%[C,h] = contour(hexDbyhexWmat,hexHbyhexTmat,log10(TsrfAmat),'LineStyle','-','LineColor','black'); 
contourf(hexDbyhexWvec,hexHbyhexTvec,TsrfMat,levels,'LineStyle','none')
[C,h] = contour(hexDbyhexWvec,hexHbyhexTvec,TsrfMat,levelVals,'LineStyle','-','LineColor','black'); 
plot(DWratioPrototype, HTratioPrototype,'wo','MarkerSize',10)
hold off;
xlabel('D/W')
ylabel('H/T')
% xlabel('D [\mum]')
% ylabel('H [\mum]')
title('Tensile Stiffness Reduction Factor')
%title('Tensile Stiffness Reduction Factor')
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



figure(5)
levels = 1000;
levelVals = [0.48, 0.5, 0.53, 0.6, 0.7, 1, 1.5];
%levelVals = 10.^(-2:1:10);
hold on;
%contourf(hexDbyhexWmat,hexHbyhexTmat,log10(TsrfAmat),levels,'LineStyle','none')
%[C,h] = contour(hexDbyhexWmat,hexHbyhexTmat,log10(TsrfAmat),'LineStyle','-','LineColor','black'); 
contourf(hexDbyhexWvec,hexHbyhexTvec,rhoAcMat,levels,'LineStyle','none')
[C,h] = contour(hexDbyhexWvec,hexHbyhexTvec,rhoAcMat,levelVals,'LineStyle','-','LineColor','black'); 
plot(DWratioPrototype, HTratioPrototype,'wo','MarkerSize',10)
hold off;
xlabel('D/W')
ylabel('H/T')
% xlabel('D [\mum]')
% ylabel('H [\mum]')
title('Corrugated Areal Density [g/m2]')
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


figure(6)
levels = 1000;
levelVals = [1.1, 1.25, 1.5, 2,0, 3.0];
%levelVals = 10.^(-2:1:10);
hold on;
%contourf(hexDbyhexWmat,hexHbyhexTmat,log10(TsrfAmat),levels,'LineStyle','none')
%[C,h] = contour(hexDbyhexWmat,hexHbyhexTmat,log10(TsrfAmat),'LineStyle','-','LineColor','black'); 
contourf(hexDbyhexWvec,hexHbyhexTvec,ADratio,levels,'LineStyle','none')
[C,h] = contour(hexDbyhexWvec,hexHbyhexTvec,ADratio,levelVals,'LineStyle','-','LineColor','black'); 
plot(DWratioPrototype, HTratioPrototype,'wo','MarkerSize',10)
hold off;
xlabel('D/W')
ylabel('H/T')
% xlabel('D [\mum]')
% ylabel('H [\mum]')
title('Corrugated-To-Flat Areal Density Ratio')
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



if saveContourMats % these are formatted for Veusz
    fNameOut = append('bsefContour_',nameOut,'Design.csv'); % no units
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',hexDbyhexWvec); fprintf(h,'\n'); % header line = x-axis = D/W ratio.  Note the initial blank (leading comma) 
    for n=length(hexHbyhexTvec):-1:1 % y-axis AND data.  y-axis = h/t ratio
        fprintf(h,'%7.5e',hexHbyhexTvec(n)); fprintf(h,' , %7.5e',BsefMat(n,:)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = append('tsrfContour_',nameOut,'Design.csv'); % no units
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',hexDbyhexWvec); fprintf(h,'\n'); % header line = x-axis = D/W ratio.  Note the initial blank (leading comma) 
    for n=length(hexHbyhexTvec):-1:1 % y-axis AND data.  y-axis = h/t ratio
        fprintf(h,'%7.5e',hexHbyhexTvec(n)); fprintf(h,' , %7.5e',TsrfMat(n,:)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = append('arealdensityContour_',nameOut,'Design.csv'); % g/m2
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',hexDbyhexWvec); fprintf(h,'\n'); % header line = x-axis = D/W ratio.  Note the initial blank (leading comma) 
    for n=length(hexHbyhexTvec):-1:1 % y-axis AND data.  y-axis = h/t ratio
        fprintf(h,'%7.5e',hexHbyhexTvec(n)); fprintf(h,' , %7.5e',rhoAcMat(n,:)); fprintf(h,'\n'); % next series of lines
    end

    fNameOut = append('arealdensityratioContour_',nameOut,'Design.csv'); % g/m2
    h = fopen(fNameOut,'w');
    fprintf(h, ',%7.5e',hexDbyhexWvec); fprintf(h,'\n'); % header line = x-axis = D/W ratio.  Note the initial blank (leading comma) 
    for n=length(hexHbyhexTvec):-1:1 % y-axis AND data.  y-axis = h/t ratio
        fprintf(h,'%7.5e',hexHbyhexTvec(n)); fprintf(h,' , %7.5e',ADratio(n,:)); fprintf(h,'\n'); % next series of lines
    end
end


