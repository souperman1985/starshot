% Calculate Accel Dist Starshot Sail Prototypes
%
% Copyright (c) 2025 Matthew Campbell
% Permission is hereby granted, free of charge, to use, copy, modify, and distribute this software for any purpose, with or without modification, provided that the above copyright notice and this permission notice appear in all copies.
% THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
%
% This script calculates the acceleration distance for the sail film that we actually fabricated
% It uses the reflectivity actually measured by Caltech (this is The January 25, 2022 sample).
% We adopt a flat circular sail model (rather than spherically curved) because we only have perpendicular reflectivity data. 
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

% Inputs
lambda0 = 1.2; % um ... laser wavelength at zero sail velocity.  
beta0 = 1e-10; % start from rest... Beta = v/cLight so 1e-10 is ~ 0.03 m/s
betaF = 0.01; % fraction of speed of light that we want to achieve
Pwr0 = 1*(1e9/1); % GW-->W = kg-m2/s3 ... laser power
Mtot = 10; % g ... total mass of sail, tethers, and chip
Msail = 5; % g ... mass of sail only

% Parameters
hexD = 77e-6; % m ... optical microscope view. Technically the diameter of the hexagons on the Si mold
hexW = 15e-6; % m ... optical microscope view. Technically the width of the ribs on the Si mold
hexH = 10e-6; % m ... 95 DRIE smooth loops
tA1 = 21e-9; % m ... ellipsometry, bottom alumina thickness
tM = 53e-9; % m ... fit to optical data, MoS2 thickness
tA2 = 51e-9; % m ... ellipsometry, top alumina thickness

rhoA1 = 3200; % kg/m3 ... Al2O3 density 
rhoM = 5060; % kg/m3 ... value for bulk MoS2
rhoA2 = 3200; % kg/m3 ... Al2O3 density 

% Non-corrugated areal density
rhoAnc = (rhoA1*tA1 + rhoM*tM + rhoA2*tA2)*(1000/1) % g/m2

% Corrugated areal density
% ribs
% M = ...
%     (rhoA1 * (tA1*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*(hexD^2-(hexD-2*tA1)^2)) + ...
%     rhoM * (tM*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD-2*tA1)^2-(hexD-2*tA1-2*tM)^2)) + ...
%     rhoA2 * (tA2*sqrt(3)/2*(hexD+hexW)^2 + hexH*sqrt(3)/2*((hexD-2*tA1-2*tM)^2-(hexD-2*tA1-2*tM-2*tA2)^2))); % kg ... unit cell mass
% trenches
M = ...
    (rhoA1 * (tA1*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1).^2-hexD.^2)) + ...
    rhoM * (tM*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1+2*tM).^2-(hexD+2*tA1).^2)) + ...
    rhoA2 * (tA2*sqrt(3)/2*(hexD+hexW).^2 + hexH*sqrt(3)/2.*((hexD+2*tA1+2*tM+2*tA2).^2-(hexD+2*tA1+2*tM).^2))); % kg ... unit cell mass  
Aplane = sqrt(3)/2 * (hexD+hexW)^2;
rhoAc = (M/Aplane)*(1000/1) % g/m2

% Sail diameter... assumes a flat (non-curved) circular sail
dSail = 2*(Msail/(pi*rhoAc))^(1/2) % m
Asail = Msail/rhoAc % m2

% Light wavelengths at the sail
betaVals = linspace(beta0, betaF, 200); % fraction of speed of light
gammaVals = 1./sqrt(1-betaVals.^2); % relativistic Lorentz factor
lambdaInVals = lambda0 * sqrt((1+betaVals)./(1-betaVals)); % um ... wavelengths of light that hit the sail

% load experimental reflectivity data
% this is from October 2010, and is the average of three measurements.
fNameRefCaltech = 'octoberAvgPlusMinusWithWL.csv'; % First column is lambda (nm), second is reflectivity (0-1)
caltechRefData = csvread(fNameRefCaltech,1,0); % read the CSV file into a matrix
interpFxnCaltechRef = griddedInterpolant(caltechRefData(:,1)*(1/1000),caltechRefData(:,2),'linear','linear'); % create function to describe data and enable interpolations.  
RefSpectrum = interpFxnCaltechRef(lambdaInVals); % arbs ... find reflectivity at the wavelengths corresponding to the beta values we defined above

% Integrate to find the acceleration distance, etc.
% distance
integrandVec1 = (1./RefSpectrum).*((gammaVals.*betaVals)./(1-betaVals).^2); % ## (no units)
integralVec1 = zeros(size(integrandVec1)); % ## ... the first values of the integral are zero
integralVec1(2:end) = cumsum((1/2)*(integrandVec1(1:end-1)+integrandVec1(2:end)).*diff(betaVals)); % ## ... this is doing cumulative trapezoidal integration 
distVals = (cLight^3/(2*Pwr0))*(Mtot*(1/1000))*integralVec1; % m ... distance values throughout acceleration
accelDist = distVals(end); % m
dLaser0 = 2*accelDist*lambda0*(1/1e6)/dSail; % m ... laser array diameter on earth
% time
integrandVec2 = (gammaVals.^3./RefSpectrum).*((1+betaVals)./(1-betaVals)); % ## (no units)
integralVec2 = zeros(size(integrandVec2)); % ## ... the first values of the integral are zero
integralVec2(2:end) = cumsum((1/2)*(integrandVec2(1:end-1)+integrandVec2(2:end)).*diff(betaVals)); % ## ... this is doing cumulative trapezoidal integration
timeVals = (cLight^2/(2*Pwr0))*(Mtot*(1/1000))*integralVec2; % s ... time values throughout acceleration at beta values, in laser's reference frame
timeValsPhoton = timeVals - distVals/cLight; % s ... time that the laser is on at beta values. 
accelTime = timeVals(end); % s
laserOnTime = timeValsPhoton(end); % s
% acceleration and pressure
accValsLaser = ((2*Pwr0*RefSpectrum)./(Mtot*(1/1000)*cLight*gammaVals.^3)).*((1-betaVals)./(1+betaVals)); % m/s2 laser ref frame
accValsSail = accValsLaser.*gammaVals.^3; % m/s2 sail ref frame
Pvals = ((Mtot*(1/1000)).*accValsSail)/Asail; % (kg*m/s2)/m2 = Pa ... photon pressure on sail 

accelDistGm = accelDist*(1/1e9) % Gm
dLaser0km = dLaser0*(1/1e3) % km
accelTimeMinutes = accelTime*(1/60) % minutes
laserOnTimeMinutes = laserOnTime*(1/60) % minutes

% approximate time to Jupiter
% assume we just travel the entire time at the final speed, because the acceleration time is so short relative to the journey
distEarthJupiter = 715e9; % m
timeEarthJupiterApprox = distEarthJupiter/(cLight*betaF) * (1/3600)*(1/24) % days
