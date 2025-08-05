function [FRmat,FTmat,FAmat,transMatMat,angleMat]=FresnelMat(lambdaVec,thetaIvec,hVec,nMat,pol)
%
% Copyright (c) 2025 Matthew Campbell
% Permission is hereby granted, free of charge, to use, copy, modify, and distribute this software for any purpose, with or without modification, provided that the above copyright notice and this permission notice appear in all copies.
% THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
%
% Transform matrix calculations
% INPUT
%   * lambdaVec = [distance] vector of wavelengths over which to run calculations
%   * thetaIvec = [radians] vector of incident angles, from normal=0, over which to run calculations
%   * hVec = [same units as lambda] vector of layer heights.  For zero thickness, use NaN.
%   * nMat = Matrix with numCol=length(hVec) and numRow=length(lambdaVec). Contains complex index of refraction data (n+i*k) as a function of wavelength.  If vacuum/air, use 1+0i. 
%   * pol = polarization, 0 for s (senkrecht=perpendicular) and 1 for p (parallel)
%
% OUTPUT
%   * FRmat = Fresnel Reflectance matrix, with lambda increasing down columns and theta increasing across rows
%   * FTmat = Fresnel Transmittance matrix
%   * FAmat = Absorption matrix
%   * transMatMat = 4-D matrix of transfer matricees, with each transfer matrix located in the 3rd and 4th dimensions
%   * angleMat = 3-D matrix of angles, with angles occupying 3rd dimension
%

% Check input
if length(hVec)==1 
    error('No interface\n')
end

% Form grid of values to check
[lambdaMat, thetaImat, hMat] = ndgrid(lambdaVec, thetaIvec, hVec); % lambda (wavelength) changes down columns, theta (incident angle) changes across rows, and layer thickness (h) changes up pages (3rd dimension).  Snell angles also change up 3rd dimension because they coincide with layers.
nMatL = repmat(reshape(nMat,[size(nMat,1),1,size(nMat,2)]),1,length(thetaIvec),1); % reshape flips nMat into 3rd dimension and repmat patterns it across rows to cover all theta values.

% Snell's law
angleMat = zeros(length(lambdaVec), length(thetaIvec), length(hVec));
angleMat(:,:,1) = thetaImat(:,:,1); 
for m=2:size(angleMat,3) % need to loop because next value of theta depends on previous
    angleMat(:,:,m) = real(asin((nMatL(:,:,m-1)./nMatL(:,:,m)).*sin(angleMat(:,:,m-1)))) - 1i*abs(imag(asin((nMatL(:,:,m-1)./nMatL(:,:,m)).*sin(angleMat(:,:,m-1))))); % rad
end

%Fresnel coefficients
if pol==0 % s/perpendicular polarization
    FrMat = (nMatL(:,:,1:end-1).*cos(angleMat(:,:,1:end-1))-nMatL(:,:,2:end).*cos(angleMat(:,:,2:end)))./(nMatL(:,:,1:end-1).*cos(angleMat(:,:,1:end-1))+nMatL(:,:,2:end).*cos(angleMat(:,:,2:end)));
    FtMat = 2*nMatL(:,:,1:end-1).*cos(angleMat(:,:,1:end-1))./(nMatL(:,:,1:end-1).*cos(angleMat(:,:,1:end-1))+nMatL(:,:,2:end).*cos(angleMat(:,:,2:end)));
elseif pol==1 % p/parallel polarization
    FrMat = (nMatL(:,:,1:end-1).*cos(angleMat(:,:,2:end))-nMatL(:,:,2:end).*cos(angleMat(:,:,1:end-1)))./(nMatL(:,:,1:end-1).*cos(angleMat(:,:,2:end))+nMatL(:,:,2:end).*cos(angleMat(:,:,1:end-1)));
    FtMat = 2*nMatL(:,:,1:end-1).*cos(angleMat(:,:,1:end-1))./(nMatL(:,:,1:end-1).*cos(angleMat(:,:,2:end))+nMatL(:,:,2:end).*cos(angleMat(:,:,1:end-1)));
else
    error('Wrong Polarization\n')
end

% What comes next depends on how many layers there are. 
if length(hVec)==2 % Special case of single interface
    FrTotMat = FrMat; % Total Fresnel coefficients. Note that in this case, FrMat is only one page long
    FtTotMat = FtMat;
    transMatMat = zeros(size(FrMat,1),size(FrMat,2),2,2);
    transMatMat(:,:,1,1) = 1./FtMat;     % upper left corner of 2x2
    transMatMat(:,:,1,2) = FrMat./FtMat; % upper right corner
    transMatMat(:,:,2,1) = FrMat./FtMat; % lower left corner
    transMatMat(:,:,2,2) = 1./FtMat;     % lower right corner
else 
    % Phase shift factors.  These ultimately have numRow=length(lambdaVec), numCol=length(thetaIvec), and numPage3=(length(hVec)-1)
    deltaMat = 2*pi*(hMat(:,:,2:end-1)./lambdaMat(:,:,2:end-1)).*nMatL(:,:,2:end-1).*cos(angleMat(:,:,2:end-1));
    deltaMat = cat(3,deltaMat,zeros(size(deltaMat,1),size(deltaMat,2))); % tack on a page of zeros on top.  Note that exp(0)=1.     
    % Form transfer matrix
    transMatMatP11 = (1./FtMat).*exp(-1i*deltaMat); % this has numRow=length(lambdaVec), numCol=length(thetaIvec), and numPage3=(length(hVec)-1)
    transMatMatP12 = (FrMat./FtMat).*exp(1i*deltaMat);
    transMatMatP21 = (FrMat./FtMat).*exp(-1i*deltaMat);
    transMatMatP22 = (1./FtMat).*exp(1i*deltaMat);
    transMatMatP = zeros(length(lambdaVec),length(thetaIvec),2,2,size(deltaMat,3));
    transMatMatP(:,:,1,1,:) = reshape(transMatMatP11,[size(transMatMatP11,1),size(transMatMatP11,2),1,1,size(transMatMatP11,3)]); % keep rows and columns the same, but flip the 3rd dimension into the 5th dimension, then stick it into the matrix of transfer matricees.
    transMatMatP(:,:,1,2,:) = reshape(transMatMatP12,[size(transMatMatP12,1),size(transMatMatP12,2),1,1,size(transMatMatP12,3)]); % 
    transMatMatP(:,:,2,1,:) = reshape(transMatMatP21,[size(transMatMatP21,1),size(transMatMatP21,2),1,1,size(transMatMatP21,3)]); % 
    transMatMatP(:,:,2,2,:) = reshape(transMatMatP22,[size(transMatMatP22,1),size(transMatMatP22,2),1,1,size(transMatMatP22,3)]); % 
    % Multiply all 5th dimensional pages of transMatMatP together. This involves 2x2 matrix multiplication. The result is a 4th dimensional matrix (transMatMat).  
    % In transMatMatP: lambda increases down columns, theta increases across columns, and at each lambda-theta combo there are several 2x2 matricees.
    % In transMatMat: lambda increases down columns, theta increases across columns, and at each lambda-theta combo there is a single 2x2 matrix.
    % [a b]   [e f]   [ae+bg af+bh]
    % [c d] x [g h] = [ce+dg cf+dh]
    % Note that we're using the dot multipler within these statements, because we're multiplying a zillion 2x2 elements in one go. 
    transMatMat = transMatMatP(:,:,:,:,1);    
    for m=2:size(deltaMat,3)
        transMatMat11 = transMatMat(:,:,1,1).*transMatMatP(:,:,1,1,m) + transMatMat(:,:,1,2).*transMatMatP(:,:,2,1,m); 
        transMatMat12 = transMatMat(:,:,1,1).*transMatMatP(:,:,1,2,m) + transMatMat(:,:,1,2).*transMatMatP(:,:,2,2,m); 
        transMatMat21 = transMatMat(:,:,2,1).*transMatMatP(:,:,1,1,m) + transMatMat(:,:,2,2).*transMatMatP(:,:,2,1,m); 
        transMatMat22 = transMatMat(:,:,2,1).*transMatMatP(:,:,1,2,m) + transMatMat(:,:,2,2).*transMatMatP(:,:,2,2,m); 
        transMatMat(:,:,1,1) = transMatMat11; % we need to do this in two steps; otherwise we're over-writing values in the middle of our calculation!
        transMatMat(:,:,1,2) = transMatMat12;
        transMatMat(:,:,2,1) = transMatMat21;
        transMatMat(:,:,2,2) = transMatMat22;
    end
    % Total Fresnel coefficients... now the resulting 2-D matricees have numRow=length(lambdaVec) and numCol=length(thetaIvec)
    FrTotMat = transMatMat(:,:,2,1)./transMatMat(:,:,1,1);
    FtTotMat = 1./transMatMat(:,:,1,1);
end

% Total Fresnel coefficients in intensity
FRmat = (abs(FrTotMat)).^2;
FTmat = (abs(FtTotMat)).^2.*(real(nMatL(:,:,end).*cos(angleMat(:,:,end)))./real(nMatL(:,:,1).*cos(angleMat(:,:,1))));
FAmat = 1 - FRmat - FTmat;

end