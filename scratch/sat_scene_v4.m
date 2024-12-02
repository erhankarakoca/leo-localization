clear all;
close all;
clc;

%% Create a satellite scenario
startTime = datetime(2020, 05, 04, 18,45,50);
% stopTime = datetime(2020, 05, 04, 19,02,20);
stopTime = datetime(2020, 05, 04, 23,02,20);
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

%% Play the sat scene defined on tle with ue 
play(satscene);

%% Any satellite can be selected from the viewer 
% Be just sure it it can access to the groun station (ue)
satNumber = 19;
referenceSatellite = "Satellite "+ string(satNumber) ;
rows = matches(intvls.Source, referenceSatellite);

tleStruct = tleread('leoSatelliteConstellation.tle');
for i = 1:30
    orbitTime = startTime+minutes(i);
    [satPos(:,i),satVelocity(:,i)] = propagateOrbit(orbitTime, ...
                                    tleStruct(satNumber), ...
                                    "OutputCoordinateFrame","fixed-frame");
    [~,~,distanceSatToUe(i)] = aer(gsUE, constellation(satNumber), orbitTime);
end
c = physconst("LightSpeed");
toaSatToUe = distanceSatToUe / c;

selectedTimeInstants = [1,7,13,22,30];

TOAs = toaSatToUe(selectedTimeInstants);

% Get all combinations of indices for pairs (i, j) where i < j
pairs = nchoosek(1:length(TOAs), 2);

% Calculate TDOA for each pair
TDOAs = arrayfun(@(row) TOAs(pairs(row, 1)) - TOAs(pairs(row, 2)), 1:size(pairs, 1));

% Add Gaussian noise to TDOA measurements (simulate clock error)
sigma = 1e-6;  % Standard deviation of the Gaussian noise (in seconds)
TDOA_noise = sigma * randn(size(TDOAs));  % Gaussian noise
TDOAs_with_noise = TDOAs + TDOA_noise;  % Add noise to TDOA

% Display results
disp('TDOA vector with noise (all combinations):');
disp(TDOAs_with_noise);

satPosxyz=satPos(:,(selectedTimeInstants))';

% True UE position (replace with actual coordinates)
actual_UE_position = ueStationECEF;  % Replace this variable with the actual UE position

% Compute radii differences for TDOAs
radii_differences = TDOAs_with_noise * c;  % Convert TDOAs to distances

% Objective function for localization
objective_func = @(p) arrayfun(@(k) ...
    abs(sqrt((p(1) - satPosxyz(pairs(k, 1), 1))^2 + (p(2) - satPosxyz(pairs(k, 1), 2))^2 + (p(3) - satPosxyz(pairs(k, 1), 3))^2) - ...
        sqrt((p(1) - satPosxyz(pairs(k, 2), 1))^2 + (p(2) - satPosxyz(pairs(k, 2), 2))^2 + (p(3) - satPosxyz(pairs(k, 2), 3))^2) - ...
        radii_differences(k)), 1:size(pairs, 1));

% Initial guess for UE position
initial_guess = mean(satPosxyz, 1);

% Solve for UE position using optimization
options = optimoptions('lsqnonlin', 'Display', 'iter');
estimated_UE_position = lsqnonlin(objective_func, initial_guess, [], [], options);

% Calculate error
localization_error = norm(estimated_UE_position - actual_UE_position);

% Display results
disp('Estimated UE Position:');
disp(estimated_UE_position);
disp('Actual UE Position:');
disp(actual_UE_position);
disp('Localization Error (meters):');
disp(localization_error);

%% 3D Visualization of Hyperboloids with Projections
figure;
hold on;
title('TDOA Localization: 3D Hyperboloids and 2D Projections');
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
grid on;
axis equal;

% Number of grid points for the 3D space
n_points = 100;
[x_range, y_range, z_range] = meshgrid(linspace(-10e6, 10e6, n_points), ...
                                       linspace(-10e6, 10e6, n_points), ...
                                       linspace(0, 10e6, n_points));

% Plot hyperboloids for each TDOA
for k = 1:size(pairs, 1)
    sat1 = satPosxyz(pairs(k, 1), :);  % First satellite in the pair
    sat2 = satPosxyz(pairs(k, 2), :);  % Second satellite in the pair

    % Calculate distances for each point in 3D space
    d1 = sqrt((x_range - sat1(1)).^2 + (y_range - sat1(2)).^2 + (z_range - sat1(3)).^2);
    d2 = sqrt((x_range - sat2(1)).^2 + (y_range - sat2(2)).^2 + (z_range - sat2(3)).^2);

    % Compute the hyperboloid for the current TDOA
    hyperboloid = abs(d1 - d2) - radii_differences(k);

    % Visualize the hyperboloid as an isosurface
    iso_val = 0;  % The zero level of the hyperboloid equation
    p = patch(isosurface(x_range, y_range, z_range, hyperboloid, iso_val));
    set(p, 'FaceColor', rand(1, 3), 'EdgeColor', 'none', 'FaceAlpha', 0.3);  % Random color

end

% Add satellite positions
scatter3(satPosxyz(:, 1), satPosxyz(:, 2), satPosxyz(:, 3), 100, 'r', 'filled');
text(satPosxyz(:, 1), satPosxyz(:, 2), satPosxyz(:, 3), ...
    compose('Sat %d', 1:size(satPosxyz, 1)), 'VerticalAlignment', 'bottom');

% Add estimated UE position
scatter3(estimated_UE_position(1), estimated_UE_position(2), estimated_UE_position(3), ...
    200, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
text(estimated_UE_position(1), estimated_UE_position(2), estimated_UE_position(3), ...
    'Estimated UE', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Add actual UE position
scatter3(actual_UE_position(1), actual_UE_position(2), actual_UE_position(3), ...
    200, 'o', 'MarkerEdgeColor', 'g', 'LineWidth', 2);
text(actual_UE_position(1), actual_UE_position(2), actual_UE_position(3), ...
    'Actual UE', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Add lighting and view adjustments
camlight;
lighting gouraud;
view(3);  % 3D view

hold off;

%% 2D Visualization of TDOA Localization (Hyperbolas)
figure;
hold on;
title('TDOA Localization with 2D Hyperbolas');
xlabel('X (meters)');
ylabel('Y (meters)');
grid on;
axis equal;

% Number of grid points for visualization
n_points = 1000;
[x_range, y_range] = meshgrid(linspace(-7e6, 7e6, n_points), linspace(-7e6, 7e6, n_points));

% Assume a constant Z-coordinate (e.g., mean of satellite heights)
z_constant = ueStationECEF(3);

% Plot hyperbolas for each TDOA
for k = 1:size(pairs, 1)
    sat1 = satPosxyz(pairs(k, 1), :);  % First satellite in the pair
    sat2 = satPosxyz(pairs(k, 2), :);  % Second satellite in the pair

    % Calculate distances
    d1 = sqrt((x_range - sat1(1)).^2 + (y_range - sat1(2)).^2 + (z_constant - sat1(3)).^2);
    d2 = sqrt((x_range - sat2(1)).^2 + (y_range - sat2(2)).^2 + (z_constant - sat2(3)).^2);

    % Compute the hyperbola for the current TDOA
    hyperbola = abs(d1 - d2) - radii_differences(k);

    % Visualize the hyperbola as a contour plot
    contour(x_range, y_range, hyperbola, [0, 0]);
end

% Add satellite positions
scatter(satPosxyz(:, 1), satPosxyz(:, 2), 100, 'r', 'filled');
text(satPosxyz(:, 1), satPosxyz(:, 2), compose('Sat %d', 1:size(satPosxyz, 1)), ...
    'VerticalAlignment', 'bottom');

% Add estimated UE position
scatter(estimated_UE_position(1), estimated_UE_position(2), 200, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
text(estimated_UE_position(1), estimated_UE_position(2), 'Estimated UE', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Add actual UE position
scatter(actual_UE_position(1), actual_UE_position(2), 200, 'o', 'MarkerEdgeColor', 'g', 'LineWidth', 2);
text(actual_UE_position(1), actual_UE_position(2), 'Actual UE', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

hold off;