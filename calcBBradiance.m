
% Script to calculate black body radiation @ temperatures
%
% Copyright (c) 2025 Matthew Campbell
% Permission is hereby granted, free of charge, to use, copy, modify, and distribute this software for any purpose, with or without modification, provided that the above copyright notice and this permission notice appear in all copies.
% THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
%
% start fresh
clear
close all
format short
format compact
clc

% Constants
cLight = 299792458; % m/s ... light speed
hPlanck = 6.62607015e-34; % J-s ... Planck's constant
kB = 1.380649e-23; % J/K ... Boltzmann constant
sigmaBB = (2*pi^5*kB^4)/(15*cLight^2*hPlanck^3); % W/m2-K4

Tvec = [300 500 750 1000 1250 1500]; % K
lambdaVec = logspace(-1,2,1000)'; % um

[lambdaMat, Tmat] = ndgrid(lambdaVec,Tvec);

IbMat = (2*hPlanck*cLight^2)./((lambdaMat*(1/1e6)).^5.*(exp((hPlanck*cLight)./(kB*lambdaMat*(1/1e6).*Tmat))-1)); % W/m3-sr

figure(1)
hold on;
plot(lambdaVec,IbMat(:,1))
plot(lambdaVec,IbMat(:,2))
plot(lambdaVec,IbMat(:,3))
plot(lambdaVec,IbMat(:,4))
plot(lambdaVec,IbMat(:,5))
plot(lambdaVec,IbMat(:,6))
hold off;
%set(gca,'yscale','log');
%set(gca,'xscale','log');

norm1 = IbMat(:,1)/max(IbMat(:,1));
norm2 = IbMat(:,2)/max(IbMat(:,2));
norm3 = IbMat(:,3)/max(IbMat(:,3));
norm4 = IbMat(:,4)/max(IbMat(:,4));
norm5 = IbMat(:,5)/max(IbMat(:,5));
norm6 = IbMat(:,6)/max(IbMat(:,6));

figure(2)
hold on;
plot(lambdaVec,norm1)
%plot(lambdaVec,norm2)
plot(lambdaVec,norm3)
plot(lambdaVec,norm4)
%plot(lambdaVec,norm5)
plot(lambdaVec,norm6)
hold off;
%set(gca,'yscale','log');
set(gca,'xscale','log');


b = 2.897771955185172661; % mm*K
lambdaPeakVec = (b./Tvec)*(1/1000)*(1e6/1) % mm --> m --> um