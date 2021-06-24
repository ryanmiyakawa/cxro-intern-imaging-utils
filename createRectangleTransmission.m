function [output, idealMask] = createRectangleTransmission(xCoords, yCoords, center, zElevation, ...
    sz, idx)

lambda_nm = 13.5; % EUV wavelength relevant for semiconductors
aoi = 6; % degrees

% Lookup common idx:
if ~isnumeric(idx)
    switch(idx)
        case 'Mo'
            delta = 0.07620647;
            beta = 0.00643542456;
        case 'Pt'
            delta = 0.109337106;
            beta = 0.0600363575;
        case 'Ru'
            delta = 0.113639966;
            beta = 0.0170648936;
        case 'Rh'
            delta = 0.124951169;
            beta = 0.0311778523;
        case 'TaN'
            delta = 0.0518;
            beta = 0.0727;
    end
    idx = 1 - delta + 1i*beta;
end

% Initialize empty matrix with resolution equal to coordinates:
output = ones(size(xCoords));

% Deconstruct inputs:
xCenter = center(1);
yCenter = center(2);

xSize   = sz(1);
ySize   = sz(2);
h       = sz(3);

n       = real(idx);
k       = imag(idx);

n_vac   = 1;
k_vac   = 0;


% compute rectangle bounds:
leftEdgeX           = xCenter - xSize/2;
rightEdgeX          = xCenter + xSize/2;
botEdgeY            = yCenter - ySize/2;
topEdgeY            = yCenter + ySize/2;

leftLeadingEdgeX    = xCenter - xSize/2 + zElevation*tand(aoi);
leftFlatTopEdgeX    = xCenter - xSize/2 + (h + zElevation)*tand(aoi);
rightFlatTopEdgeX   = xCenter + xSize/2 + zElevation*tand(aoi);
rightFallingEdgeX   = xCenter + xSize/2 + (h + zElevation)*tand(aoi);

% Create ideal mask:
region_ideal = leftEdgeX <= xCoords & xCoords <= rightEdgeX & ...
    botEdgeY <= yCoords & yCoords <= topEdgeY;
idealMask = ones(size(xCoords));
idealMask(region_ideal) = 0;

flatTopRayLength = h/cosd(aoi);

% Split into 4 regions:  Vac, flat_top, leading_edge, falling_edge
region_flat_top = (botEdgeY <= yCoords & yCoords <= topEdgeY) & ... % Y-coord constraints
    leftFlatTopEdgeX <= xCoords & xCoords <= rightFlatTopEdgeX; % X-coord constraints

region_leadingEdge = (botEdgeY <= yCoords & yCoords <= topEdgeY) & ... % Y-coord constraints
    leftLeadingEdgeX <= xCoords & xCoords <= leftFlatTopEdgeX; % X-coord constraints

region_fallingEdge = (botEdgeY <= yCoords & yCoords <= topEdgeY) & ... % Y-coord constraints
    rightFlatTopEdgeX <= xCoords & xCoords <= rightFallingEdgeX; % X-coord constraints

region_vac = ~(region_flat_top | region_leadingEdge | region_fallingEdge);


% Compute total effective ray path thickness for absorber and vacuum
% regions:
absorber_thickness = zeros(size(xCoords));
absorber_thickness(region_flat_top) = flatTopRayLength;

absorber_thickness(region_leadingEdge) = ...
    (xCoords(region_leadingEdge) - leftLeadingEdgeX) / (leftFlatTopEdgeX - leftLeadingEdgeX) ...
    * flatTopRayLength;

absorber_thickness(region_fallingEdge) = ...
    (rightFallingEdgeX - xCoords(region_fallingEdge)) / (rightFallingEdgeX - rightFlatTopEdgeX) ...
    * flatTopRayLength;

vacuum_thickness = flatTopRayLength - absorber_thickness;


% Compute amplitude and phase in each region:
output = exp(-2*pi*k * absorber_thickness/lambda_nm) .* exp(2i*pi*n*absorber_thickness/lambda_nm) ...
    .*  exp(2i*pi*n_vac*(vacuum_thickness/lambda_nm));





