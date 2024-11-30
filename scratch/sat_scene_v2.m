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

TOAs = [toaSatToUe(1),toaSatToUe(5),toaSatToUe(8),toaSatToUe(12),toaSatToUe(15)];

TDOAs = [TOAs(1)-TOAs(2),      % 1 2
         TOAs(1)-TOAs(3),      % 1 3 
         TOAs(1)-TOAs(4),      % 1 4
         TOAs(2)-TOAs(3),      % 2 3
         TOAs(2)-TOAs(4),      % 2 4
         TOAs(3)-TOAs(4)];     % 3 4 
             

satPosxyz=satPos(:,[1,5,8,12,15])';

radii = TOAs * c;  % Convert TOAs to radii

% Define objective function for intersection point
objective_func = @(p) arrayfun(@(i) ...
    sqrt((p(1) - satPosxyz(i, 1))^2 + (p(2) - satPosxyz(i, 2))^2 + (p(3) - satPosxyz(i, 3))^2) - radii(i), ...
    1:length(radii));

% Initial guess for intersection point (e.g., origin or weighted center of satellites)
initial_guess = mean(satPosxyz, 1);

% Solve using nonlinear least squares
options = optimoptions('lsqnonlin', 'Display', 'iter');
intersection_point = lsqnonlin(objective_func, initial_guess, [], [], options);

% Actual UE position (example coordinates; update with actual values if available)
actual_UE_position = ueStationECEF;  % Replace `ueStationECEF` with actual UE ECEF coordinates

% Calculate error
error_vector = intersection_point - actual_UE_position;
error_magnitude = norm(error_vector);


% Display results
disp('Intersection Point (Estimated UE Location):');
disp(intersection_point);
disp('Actual UE Position:');
disp(actual_UE_position);
disp(error_magnitude)

% Visualization
figure;
hold on;
title('TOA-Based Spheres, Intersection, and Actual UE Position');
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
grid on;
axis equal;

% Plot spheres
n_points = 1000;
[theta, phi] = meshgrid(linspace(0, 2*pi, n_points), linspace(0, pi, n_points));
for i = 1:length(radii)
    % Sphere coordinates
    x_sphere = radii(i) * sin(phi) .* cos(theta) + satPosxyz(i, 1);
    y_sphere = radii(i) * sin(phi) .* sin(theta) + satPosxyz(i, 2);
    z_sphere = radii(i) * cos(phi) + satPosxyz(i, 3);
    surf(x_sphere, y_sphere, z_sphere, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Plot satellite positions
scatter3(satPosxyz(:, 1), satPosxyz(:, 2), satPosxyz(:, 3), 100, 'r', 'filled');
text(satPosxyz(:, 1), satPosxyz(:, 2), satPosxyz(:, 3), ...
    compose('Sat %d', 1:length(radii)), 'VerticalAlignment', 'bottom');

% Plot intersection point
scatter3(intersection_point(1), intersection_point(2), intersection_point(3), 200, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
text(intersection_point(1), intersection_point(2), intersection_point(3), 'Intersection', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Plot actual UE position
scatter3(actual_UE_position(1), actual_UE_position(2), actual_UE_position(3), 200, 'o', 'MarkerEdgeColor', 'g', 'LineWidth', 2);
text(actual_UE_position(1), actual_UE_position(2), actual_UE_position(3), 'UE', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

hold off;