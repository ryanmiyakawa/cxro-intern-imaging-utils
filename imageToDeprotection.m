function [thresholded, latentImage] = imageToDeprotection(imgIntens, sz_um, dose_mJpcm2, ...
                                blur_nm, rate, threshold)

mpm addpath

[sr, sc] = size(imgIntens);

if (size(sz_um) == 1)
    sx = sz_um;
    sy = sz_um;
else
    sx = sz_um(2);
    sy = sz_um(1); 
end

pxArea = sx * sy/(sr*sc);
pxLen = sx/sr;

% Get dose per pixel

totalImageSize_um2 = sx * sy;
totalImageSize_cm2 = totalImageSize_um2 / 10^8;

% Energy in entire image
energyInImage_mJ = dose_mJpcm2 * totalImageSize_cm2;
energyInImage_J = energyInImage_mJ/1000;

% Photons in entire image
Ep_J = 6.626E-34*300000000/0.0000000135;
Np = energyInImage_J/Ep_J;


% Photons per pixel:
ppp = Np / (sr * sc);

% Max img 
mxPixel = max(abs(imgIntens(:)));
imgIntens = imgIntens/mxPixel * ppp;

% Add photon noise noise
imgIntens = imgIntens + sqrt(imgIntens).*randn(size(imgIntens));
imgIntens = abs(real(imgIntens));

% Normalize to area
cellSize = 2/1000;
pxRat = pxLen/cellSize;
Nr = round(pxRat * sr);
Nc = round(pxRat * sc);

scaledImg = bin2(imgIntens, Nr, Nc);

% log events
events = rand(size(scaledImg)).*scaledImg;
events(events > rate) = 1;
events(events ~= 1) = 0;


% create grid:
xidx = linspace(-sx/2, sx/2, sc) * 1000;
yidx = linspace(-sy/2, sy/2, sr) * 1000;

[X, Y] = meshgrid(xidx, yidx);

% create gaussian:
gaus = (exp(-((X).^2 + (Y).^2)/(2*blur_nm^2)));

% resize and conv
events = bin2(events, sr, sc);
latentImage = cconv2(gaus, events);

mxVal = max(latentImage(:));

thresholded = latentImage;
thresholded(thresholded < mxVal * threshold) = 0;
thresholded(thresholded ~= 0) = 1;








