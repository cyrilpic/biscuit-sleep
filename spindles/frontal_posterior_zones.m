function [fZones,pZones] = frontal_posterior_zones(Electrodes)
%FRONTAL_POSTERIOR_ZONES Return the channels involved in the frontal and
%posterior zones

% adjust the electrode coordinates
Th = pi/180*[Electrodes.theta];        % Calculate theta values from x,y,z e_loc
Rd = [Electrodes.radius];              % Calculate radian values from x,y,z e_loc

x = Rd.*cos(Th);                            % Calculate 2D projected X
y = Rd.*sin(Th);                            % Calculate 2D projected Y

% Squeeze the coordinates into a -0.5 to 0.5 box
intrad = min(1.0,max(abs(Rd))); intrad = max(intrad,0.5);
squeezefac = 0.5/intrad;
x = x * squeezefac; y = y * squeezefac;

% fZones
distance_from_center = 0.10; % x and y distance from the center in each direction
circle_radius = 0.08;

centers_x = -0.16;
%centers_y = centers_x';
centers_y = (-1:1) * distance_from_center;

% pre-allocate
fZones = false(1, length(x));

% loop for each region
for n = 1 : 3
    % calculate distance to center for each electrode
    distances = ((x + centers_x).^ 2 ...
        + (y + centers_y(n)).^ 2).^ 0.5;
    % take electrodes within maximum distance
    fZones(n, :) = distances < circle_radius;
end

% pZones

distance_from_center = 0.10; % x and y distance from the center in each direction
circle_radius = 0.08;

centers_x = [0.08, 0.10, 0.08];
%centers_y = centers_x';
centers_y = (-1:1) * distance_from_center;

% pre-allocate
pZones = false(1, length(x));

% loop for each region
for n = 1 : 3
    % calculate distance to center for each electrode
    distances = ((x + centers_x(n)).^ 2 ...
        + (y + centers_y(n)).^ 2).^ 0.5;
    % take electrodes within maximum distance
    pZones(n, :) = distances < circle_radius;
end

end

