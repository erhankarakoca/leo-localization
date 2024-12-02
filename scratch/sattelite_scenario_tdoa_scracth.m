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

TDOAs = [toaSatToUe(1)-toaSatToUe(5),       % 1 2
         toaSatToUe(1)-toaSatToUe(10),      % 1 3 
         toaSatToUe(1)-toaSatToUe(15),      % 1 4
         toaSatToUe(5)-toaSatToUe(10),      % 2 3
         toaSatToUe(5)-toaSatToUe(15),      % 2 4
         toaSatToUe(10)-toaSatToUe(15)];     % 3 4 
             

satPosxyz=satPos(:,[1,5,10,15])';

% Visualize hyperboloids for the first three TDOAs
figure;
hold on;
title('TDOA-Based Hyperboloids in 3D Space');
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
grid on;

% Plot each hyperboloid
n_points = 500; % Number of grid points for mesh
[X, Y] = meshgrid(linspace(-7e6, 7e6, n_points), linspace(-7e6, 7e6, n_points));

index = [1,2 ; 1,3 ; 1,4; 2,3; 2,4; 3,4];
for k = 1:6 % Visualize first three TDOAs
    sat1 = satPosxyz(index(k,1), :);  % First satellite in the pair
    sat2 = satPosxyz(index(k,2), :);  % Second satellite in the pair
    
    % Distance difference for this TDOA
    delta_d = TDOAs(k) * c; % Convert TDOA back to distance
    
    % Hyperboloid equation
    Z = sqrt((sqrt((X - sat1(1)).^2 + (Y - sat1(2)).^2 ) - ...
              sqrt((X - sat2(1)).^2 + (Y - sat2(2)).^2 ) - delta_d).^2);
    
    % Plot the hyperboloid
    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

% Plot the satellite positions
scatter3(satPosxyz(:, 1), satPosxyz(:, 2), satPosxyz(:, 3), 50, 'r', 'filled');
text(satPosxyz(:, 1), satPosxyz(:, 2), satPosxyz(:, 3), ...
    compose('Sat %d', 1:4), 'VerticalAlignment', 'bottom');

% Plot the UE position
scatter3(ueStationECEF(1), ueStationECEF(2), ueStationECEF(3), 50, 'b', 'filled');
text(ueStationECEF(1), ueStationECEF(2), ueStationECEF(3), 'UE', 'VerticalAlignment', 'bottom');

hold off;
