clear all;
close all;
clc;
%% Create a satellite scenario
startTime = datetime(2020, 05, 04, 18,45,50);
stopTime = datetime(2020, 05, 04, 19,02,20);
sampleTime=10;
satscene = satelliteScenario(startTime,stopTime,sampleTime);

% Add satellites from TLE file.
tleFile = "leoSatelliteConstellation.tle";
constellation = satellite(satscene, tleFile);

% Define ue to be ground station Lat Long Alt;
ueStationLLA = [40.786648, 29.449502, 182];
ueStationECEF = lla2ecef(ueStationLLA);

gsUE = groundStation(satscene, ...
                     "Latitude",  ueStationLLA(1), ...
                     "Longitude", ueStationLLA(2), ...
                     "Altitude",  ueStationLLA(3));


%% Find the access intervals
ac = access(constellation,gsUE);
intvls = accessIntervals(ac);
% play(satscene)

%% Play the sat scene defined on tle with ue 
% play(satscene);

%% Any satellite can be selected from the viewer 
% Be just sure it it can access to the groun station (ue)
satNumber = 19;
referenceSatellite = "Satellite "+ string(satNumber) ;
rows = matches(intvls.Source, referenceSatellite);

% find specific acces-in times for the defined satellite
% startEndTimes = table2array(intvls(rows, ["StartTime", "EndTime"]));
% startTime = startEndTimes(1);
% endTime = startEndTimes(2);

%% Read tle to calculate specific satellite position and speed
% sampleTime = 60;
% sc = satelliteScenario(startTime,endTime,sampleTime);
% constellationInterval = satellite(sc, tleFile);
% gs = groundStation(sc, ueStation(1), ueStation(2));
% acInterval = access(constellationInterval,gs);

tleStruct = tleread('leoSatelliteConstellation.tle');
for i = 1:15
    orbitTime = startTime+minutes(i);
    [satPos(:,i),satVelocity(:,i)] = propagateOrbit(orbitTime, ...
                                    tleStruct(satNumber), ...
                                    "OutputCoordinateFrame","fixed-frame");
    [~,~,distanceSatToUe(i)] = aer(gsUE, constellation(satNumber), orbitTime);
end
c = physconst("LightSpeed");
toaSatToUe = distanceSatToUe / c ;

TOAs = [toaSatToUe(1),toaSatToUe(5),toaSatToUe(10),toaSatToUe(15)];

TDOAs = [TOAs(1)-TOAs(2),      % 1 2
         TOAs(1)-TOAs(3),      % 1 3 
         TOAs(1)-TOAs(4),      % 1 4
         TOAs(2)-TOAs(3),      % 2 3
         TOAs(2)-TOAs(4),      % 2 4
         TOAs(3)-TOAs(4)];     % 3 4 
             

satPosxyz=satPos(:,[1,5,10,15])';

% Visualization of hyperboloids for TDOA differences
figure;
hold on;
title('TDOA-Based Hyperboloids in 3D Space');
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
grid on;

% Number of grid points for mesh
n_points = 100;  % Adjust grid resolution for better performance

% Create meshgrid for X, Y coordinates (Z is fixed for slices)
[x_range, y_range, z_range] = meshgrid(linspace(-7e6, 7e6, n_points), linspace(-7e6, 7e6, n_points), linspace(-7e6, 7e6, n_points));

% Set Z values for slices (you can modify this to get a range of slices)
z_slices = linspace(-7e6, 7e6, 10);  % Define multiple Z-slices to visualize hyperboloids

% Define the index of the satellite pairs for TDOA visualization
index = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];  % Index for pairs (1-2, 1-3, 1-4, etc.)

for k = 1:6  % Visualize all six TDOA measurements (1-2, 1-3, 1-4, 2-3, 2-4, 3-4)
    % Select the two satellites for the current pair
    sat1 = satPosxyz(index(k,1), :);  % First satellite in the pair
    sat2 = satPosxyz(index(k,2), :);  % Second satellite in the pair
    
    % Distance difference for this TDOA
    delta_d = TDOAs(k) * c;  % Convert TDOA back to distance

   
    % Calculate the Euclidean distances from all points in the grid to the two satellites
    d1 = sqrt((x_range - sat1(1)).^2 + (y_range - sat1(2)).^2 + (z_range - sat1(3)).^2);
    d2 = sqrt((x_range - sat2(1)).^2 + (y_range - sat2(2)).^2 + (z_range - sat2(3)).^2);
    
    % Hyperboloid equation: |d1 - d2| = delta_d
    hyperboloid_eq = abs(d1 - d2) - delta_d;


    % Plot the hyperboloid using contour3 for the current Z-slice
    contour3(x_range, y_range, hyperboloid_eq, [0, 0], 'LineWidth', 2);
  

    % Plot the hyperboloid using surf, contour, or other visualization method
    surf(x_range, y_range, hyperboloid_eq, 'EdgeColor', 'none', 'FaceAlpha', 1); % Plot the surface directly
        

end

% Plot satellite positions in 3D space
scatter3(satPosxyz(:, 1), satPosxyz(:, 2), satPosxyz(:, 3), 50, 'r', 'filled');
text(satPosxyz(:, 1), satPosxyz(:, 2), satPosxyz(:, 3), ...
    compose('Sat %d', 1:4), 'VerticalAlignment', 'bottom');

% Plot the UE position in 3D
scatter3(ueStationECEF(1), ueStationECEF(2), ueStationECEF(3), 50, 'b', 'filled');
text(ueStationECEF(1), ueStationECEF(2), ueStationECEF(3), 'UE', 'VerticalAlignment', 'bottom');

hold off;
