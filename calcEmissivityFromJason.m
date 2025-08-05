% Calculate emissivity
%
% Copyright (c) 2025 Matthew Campbell
% Permission is hereby granted, free of charge, to use, copy, modify, and distribute this software for any purpose, with or without modification, provided that the above copyright notice and this permission notice appear in all copies.
% THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
%
% Matthew Campbell
% 2025-02-21

% Start fresh
clear all
close all
format long
format compact
clc

% Constants
cLight = 299792458; % m/s ... light speed
hPlanck = 6.62607015e-34; % J-s ... Planck's constant
kB = 1.380649e-23; % J/K ... Boltzmann constant
sigmaBB = (2*pi^5*kB^4)/(15*cLight^2*hPlanck^3); % W/m2-K4

Toptical = 1000 % K ... for evaluating optical properties. 

% Data selection
dataSelect = 3;

fNameJason = '';
if dataSelect == 1;     fNameJason = 'BrewerAngleDependentAbs_withMunkhbat';
elseif dataSelect == 2; fNameJason = 'Chang'; 
elseif dataSelect == 3; fNameJason = 'Experimental'; 
elseif dataSelect == 4; fNameJason = 'Ilic'; 
elseif dataSelect == 5; fNameJason = 'Lien';
elseif dataSelect == 6; fNameJason = 'Optimized_newDesign'; 
elseif dataSelect == 7; fNameJason = 'Salary'; 
elseif dataSelect == 8; fNameJason = 'Santi'; 
elseif dataSelect == 9; fNameJason = 'Taghavi'; 
end

%% Calculate 

% Read file
M = readmatrix([fNameJason,'.xlsx']);

% parse data
thetaVece = M(2:end,1)'*(pi/180); % first column, deg-->rad, shift from column to row
lambdaOutVals = M(1,2:end)'*(1/1000); % first row, nm-->um, shift from row to column
AmateJason = M(2:end,2:end)'; % now angle increases across rows, wavelength increases down columns 

% Add 90-degree angle with 0 absorption
thetaMateEx = repmat([thetaVece 90*(pi/180)],[length(lambdaOutVals),1]); % tack on a 90-degree angle and pattern over all lambda values
AmateExJason = cat(2,AmateJason,zeros(length(lambdaOutVals),1)); % absorption is zero at 90 degrees

% emissivity calc
integrandMatE = AmateExJason.*sin(thetaMateEx).*cos(thetaMateEx); 
EmisSpectrumJason = (2*pi)*(1/pi)*sum((1/2)*(integrandMatE(:,1:end-1)+integrandMatE(:,2:end)).*diff(thetaMateEx,1,2),2); % ## ... this is doing trapezoidal integration in the 2nd dimension (along rows, across angles)
% Effective emissivity for reference
% Black body spectrum
IbVec = (2*hPlanck*cLight^2)./((lambdaOutVals*(1/1e6)).^5.*(exp((hPlanck*cLight)./(kB*lambdaOutVals*(1/1e6)*Toptical))-1)); % W/m3-sr
integrandVecDen = IbVec; % W/m3-sr
integralDen = sum((1/2)*(integrandVecDen(1:end-1)+integrandVecDen(2:end)).*diff(lambdaOutVals)); % ## ... trapezoidal integration
% https://en.wikipedia.org/wiki/Kirchhoff%27s_law_of_thermal_radiation
integrandVecNum = IbVec.*EmisSpectrumJason; % W/m3-sr
integralNum = sum((1/2)*(integrandVecNum(1:end-1)+integrandVecNum(2:end)).*diff(lambdaOutVals)); % ## ...trapezoidal integration
effEpsJason = integralNum/integralDen; % effective emissivity value design.
effEpsJason
