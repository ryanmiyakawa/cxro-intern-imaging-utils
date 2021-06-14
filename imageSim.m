function filtered = imageSim(object, size_um, lambda_nm, NA)


% Check if 1D or 2D:
[sr, sc] = size(object);


if (sr == 1) % then 1D
    
    objSpectrum = (fft(object));
    F0 = 1/(size_um); % Smallest spatial frequency
    fCutoff = 1/(lambda_nm/NA); % Largest frequency supported by lens
    pxCutoff = ceil(fCutoff/F0);
    
    % Filter frequencies that miss the lens
    objSpectrum(pxCutoff + 1 : end - (pxCutoff - 1)) = 0;
  
    filtered = ifft(objSpectrum);
    filtered(abs(filtered) < 1e-8) = 0;
    
else
    if length(size_um) ~= 2
        error('Please enter an array of sizes [size_x, size_y]');
    end
    
    F0x = 1/(size_um(2)); % Smallest spatial frequency
    F0y = 1/(size_um(1)); % Smallest spatial frequency
    
    fCutoff = 1/(lambda_nm/NA); % Largest frequency supported by lens
    
    xidx = 1:sc;
    yidx = 1:sr;
    
    [Y, X] = meshgrid(yidx, xidx);
    
    filtMask = sqrt((Y*F0y).^2 + (X*F0x).^2) < fCutoff;
    
    filtMask = filtMask | flipud(filtMask);
    filtMask = filtMask | fliplr(filtMask);
    
    objSpectrum = fft2(object);
    filtered = ifft2(objSpectrum.*filtMask);

    
end





