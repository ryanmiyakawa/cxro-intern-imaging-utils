function out = imageToDeprotection(img, sz, dose, blur, rate)
[sr, sc] = size(img);

if (size(sz) == 1)
    sx = sz;
    sy = sz;
else
    sx = sz(2);
    sy = sz(1); 
end

% Get dose per pixel

totalImageSize_um2 = sx * sy;

totalImageSize_cm2 = totalImageSize_um2 / 10^8;

% Energy in entire image
energyInImage_mJ = dose * totalImageSize_cm2;
energyInImage_J = energyInImage_mJ/1000;

% Photons in entire image
Ep_J = 6.626E-34*300000000/0.0000000135;
Np = energyInImage_J/Ep_J;


% Photons per pixel:
ppp = Np / (sr * sc);




out = ppp;

