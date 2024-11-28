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
play(satscene);

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
for i = 1:4
    orbitTime = startTime+minutes(i);
    [satPos(:,i),satVelocity(:,i)] = propagateOrbit(orbitTime, ...
                                    tleStruct(satNumber), ...
                                    "OutputCoordinateFrame","fixed-frame");
    [~,~,distanceSatToUe(i)] = aer(gsUE, constellation(satNumber), orbitTime);
end
c = physconst("LightSpeed");
toaSatToUe = distanceSatToUe / c ;

TDOAs = [toaSatToUe(1)-toaSatToUe(2), 
         toaSatToUe(1)-toaSatToUe(3),
         toaSatToUe(1)-toaSatToUe(4)];

satPosxyz=satPos';
% [x,y,z] = tdoa_localization(satPosxyz, TDOAs);

% Example input (satellite positions, initial guess for UE position, TDOAs)
% sat_coords = [x1, y1, z1; x2, y2, z2; x3, y3, z3; ...];  % Each row is [x, y, z] of a satellite
% TDOAs = [tdoa_1, tdoa_2, tdoa_3, ...];  % TDOAs for each satellite pair

% Initial guess for UE position (start with a guess)
initial_guess = [0, 0, 0];  % [x, y, z]

% Call lsqnonlin to find the UE position
options = optimset('Display', 'off');
[UE_position, ~] = lsqnonlin(@(P_UE) tdoa_objective(P_UE, satPosxyz, TDOAs, c), initial_guess, [], [], options);

disp('Estimated UE position:');
disp(UE_position);

function [x_UE, y_UE, z_UE] = tdoa_localization(sat_coords, TDOAs)
    % sat_coords: 4x3 matrix of satellite coordinates in ECEF
    % TDOAs: vector of TDOA values between satellites (e.g., [TDOA12, TDOA13, TDOA14])
    
    % Speed of light
    c =  physconst("LightSpeed");  % meters per second
    
    % Objective function to minimize (sum of squared residuals)
    objective_function = @(P_UE) sum( ...
        (sqrt((P_UE(1) - sat_coords(:,1)).^2 + (P_UE(2) - sat_coords(:,2)).^2 + (P_UE(3) - sat_coords(:,3)).^2) - ...
         sqrt((P_UE(1) - sat_coords(2:end,1)).^2 + (P_UE(2) - sat_coords(2:end,2)).^2 + (P_UE(3) - sat_coords(2:end,3)).^2) - ...
         c * TDOAs).^2);
    
    % Initial guess for UE position (e.g., near the first satellite)
    initial_guess = [0, 0, 0];
    
    % Solve using nonlinear least squares
    options = optimoptions('lsqnonlin', 'Display', 'off');
    [UE_position, ~] = lsqnonlin(objective_function, initial_guess, [], [], options);
    
    % Extract UE coordinates
    x_UE = UE_position(1);
    y_UE = UE_position(2);
    z_UE = UE_position(3);
end


function [error] = tdoa_objective(P_UE, sat_coords, TDOAs, c)
    % P_UE: Estimated UE position [x, y, z]
    % sat_coords: Matrix of satellite coordinates (each row is a satellite [x, y, z])
    % TDOAs: Time difference of arrival values
    % c: Speed of light (constant)

    % Calculate the distance from the UE to each satellite
    distances = sqrt((P_UE(1) - sat_coords(:,1)).^2 + ...
                     (P_UE(2) - sat_coords(:,2)).^2 + ...
                     (P_UE(3) - sat_coords(:,3)).^2);

    % Calculate the TDOA errors (distance differences)
    % Compute the pairwise differences
    distance_diff = distances(1) - distances(2:end);  % Difference with the first satellite
    
    % Calculate the expected time differences based on distances
    expected_TDOAs = distance_diff / c;  % Convert to time difference
    
    % Compute the error (difference between observed and expected TDOAs)
    error = (expected_TDOAs - TDOAs).^2;
    
    % Sum the squared errors
    error = sum(error);
end



