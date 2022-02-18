% Script to plot the different detection regions for frontal and posterior
% spindles

% Load Info structure

%% Original method

frontal_channels = {'E36', 'E29', 'E22', 'E14', 'E5', 'E224',...
    'E215', 'E207', 'E24', 'E30', 'E23', 'E16', 'E6', 'E7', 'E15'};
parietal_channels = {'E142', 'E141', 'E140', 'E127', 'E118',...
    'E119', 'E109', 'E99', 'E88', 'E129', 'E128', 'E100', 'E101', 'E110'};

% convert to indices
frontal_channels = any(cell2mat(cellfun(@(x) strcmp({Info.Electrodes.labels}', x), frontal_channels, ...
    'uniform', false)), 2);
parietal_channels = any(cell2mat(cellfun(@(x) strcmp({Info.Electrodes.labels}', x), parietal_channels, ...
    'uniform', false)), 2);


H = csc_Topoplot(nan(size(Info.Electrodes)), Info.Electrodes);

set(H.Channels(frontal_channels), ...
    'visible', 'on', ...
    'fontSize', 15, ...
    'string', num2str('f'));

set(H.Channels(parietal_channels), ...
    'visible', 'on', ...
    'fontSize', 15, ...
    'string', num2str('p'));


%% New method

H = csc_Topoplot(nan(size(Info.Electrodes)), Info.Electrodes);


% adjust the electrode coordinates
Th = pi/180*[Info.Electrodes.theta];        % Calculate theta values from x,y,z e_loc
Rd = [Info.Electrodes.radius];              % Calculate radian values from x,y,z e_loc

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

set(H.Channels(any(fZones,1)), ...
    'visible', 'on', ...
    'fontSize', 15, ...
    'string', num2str('f'));

plot(-centers_y, -centers_x, 'rx')

for n=1:3
    pos = [[-centers_y(n), -centers_x]-circle_radius, 2*circle_radius 2*circle_radius];
    rectangle('Position', pos, 'Curvature', [1, 1], 'FaceColor', [0. 0.5 0.5 0.3]);
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

set(H.Channels(any(pZones,1)), ...
    'visible', 'on', ...
    'fontSize', 15, ...
    'string', num2str('p'));

plot(-centers_y, -centers_x, 'rx')

for n=1:3
    pos = [[-centers_y(n), -centers_x(n)]-circle_radius, 2*circle_radius 2*circle_radius];
    rectangle('Position', pos, 'Curvature', [1, 1], 'FaceColor', [0.5 0.5 0.0 0.3]);
end