clear all;
close all;
clc;
%% Create a satellite scenario
startTime = datetime(2020, 05, 04, 18,45,50);
stopTime = datetime(2020, 05, 04, 19,02,20);
sampleTime = 10;
satscene = satelliteScenario(startTime,stopTime,sampleTime);

% Add satellites from TLE file.
tleFile = "leoSatelliteConstellation.tle";
constellation = satellite(satscene, tleFile);

% Define ue to be ground station Lat Long Alt;
ueStationLLA = [40.786648, 29.449502, 182];
% ECEF correspondings
ueStationECEF = lla2ecef(ueStationLLA);
% Add to the satellite scene
gsUE = groundStation(satscene, ...
                     "Latitude",  ueStationLLA(1), ...
                     "Longitude", ueStationLLA(2), ...
                     "Altitude",  ueStationLLA(3));


%% Find the access intervals
ac = access(constellation,gsUE);
intvls = accessIntervals(ac);

%% Extract accessed satellites for specific date-time
% Generate a random date-time within the startTime and stopTime
totalSamples = seconds(stopTime - startTime) / sampleTime;

%% Generate a random sample index (integer between 0 and totalSamples - 1)
randomSampleIndex = randi([0, totalSamples - 1]);

%% Calculate random date-time based on sampleTime
randomDateTime = startTime + seconds(randomSampleIndex * sampleTime);
randomDateTime = datetime(randomDateTime, 'TimeZone', 'UTC');
%% Initialize an array to hold accessed satellites
accessedSatellites = [];

% Loop through the intervals to check if the random date-time is within any access interval
for i = 1:height(intvls)
    % Extract start and end time from the table
    accessStartTime = intvls{i, 4}; % 4th column: Access start date-time
    accessEndTime = intvls{i, 5};   % 5th column: Access end date-time
    
    % Check if the random date-time is within the access interval
    if randomDateTime >= accessStartTime && randomDateTime <= accessEndTime
        % If true, append the satellite name (1st column) to the array
        accessedSatellites = [accessedSatellites; intvls{i, 1}];
    end
end

%% Display the results
if isempty(accessedSatellites)
    disp('No satellites accessed at the random date-time.');
else
    fprintf('Random Date-Time: %s\n', datestr(randomDateTime));
    disp('Satellites accessed at this time:');
    disp(accessedSatellites);
end
%% Play the sat scene defined on tle with ue 
% play(satscene);

%% Any satellite can be selected from the viewer 
% Be just sure it it can access to the groun station (ue)


tleStruct = tleread('leoSatelliteConstellation.tle');
% rows = matches(tleStruct.Name, accessedSatellites);

% Extract satellite names from tleStruct
satelliteNamesInTLE = {tleStruct.Name}';
% Find indices in tleStruct corresponding to accessedSatellites
indicesInTLE = find(matches(string(satelliteNamesInTLE), accessedSatellites));

orbitTime = randomDateTime;

% Access the corresponding TLE data for the accessed satellites
accessedTLEStruct = tleStruct(indicesInTLE);

% Extract orbit propagation and calculate parameters
[accessedSatPositions, accessedSatVelocities] = propagateOrbit(orbitTime, ...
                                                               accessedTLEStruct, ...
                                                               "OutputCoordinateFrame", "fixed-frame");

% Compute distances to UE for all accessed satellites
[accessedSatAzimuths, accessedSatElevations, accessedSatDistances] = aer(gsUE, ...
                                                                         constellation(indicesInTLE), ...
                                                                         orbitTime);

accessedSatPositions = squeeze(accessedSatPositions);
accessedSatPositions = accessedSatPositions';
accessedSatVelocities = squeeze(accessedSatVelocities);

accessedSatDistances = squeeze(accessedSatDistances);
accessedSatAzimuths = squeeze(accessedSatAzimuths);
accessedSatElevations = squeeze(accessedSatElevations);

c = physconst("LightSpeed");
TOAs = accessedSatDistances / c ;

%% Initial guess for UE position
% initial_guess = mean(satPosxyz, 1);
initial_guess = lla2ecef([39.284593, 33.421097, 887]);
% initial_guess = [0,0,0];
%%
% Define GDOP threshold (tune as needed)
gdopThreshold = 5;

% Initialize storage for selected subsets
% selectedCombinations = [];

% Iterate over all subsets of satellite positions
numSats = size(accessedSatPositions, 1);
allCombinations = nchoosek(1:numSats, 4); % For 4 satellites at a time

for i = 1:size(allCombinations, 1)
    subset = allCombinations(i, :);
    subsetPositions = accessedSatPositions(subset, :);
    gdop(i) = calculateGDOP(subsetPositions, initial_guess);
    % if gdop(i) <= gdopThreshold
    %     selectedCombinations = [selectedCombinations; subset]; %#ok<AGROW>
    % end
end

% Use selected combinations for TDOA calculations
% disp('Selected satellite combinations based on GDOP:');
% disp(selectedCombinations);

[~, index] = min(gdop);
selectedSatIndices = allCombinations(index, :);
selectedSatPositions = accessedSatPositions(selectedSatIndices, :);

% Compute TOAs for selected satellites
selectedTOAs = TOAs(selectedSatIndices);

% Recompute pairs for selected satellites
pairs = nchoosek(1:length(selectedTOAs), 2);

% Calculate TDOA for each pair
TDOAs = arrayfun(@(row) selectedTOAs(pairs(row, 1)) - selectedTOAs(pairs(row, 2)), ...
                 1:size(pairs, 1));

%%
% True UE position (replace with actual coordinates)
actual_UE_position = ueStationECEF;  % Replace this variable with the actual UE position

% Compute radii differences for TDOAs
radii_differences = TDOAs * c;  % Convert TDOAs to distances


%% Objective Functions
objective_func = @(p) ((arrayfun(@(k) ...
    abs(sqrt((p(1) - selectedSatPositions(pairs(k, 1), 1))^2 + ...
             (p(2) - selectedSatPositions(pairs(k, 1), 2))^2 + ...
             (p(3) - selectedSatPositions(pairs(k, 1), 3))^2) - ...
        ...
        sqrt((p(1) - selectedSatPositions(pairs(k, 2), 1))^2 + ...
             (p(2) - selectedSatPositions(pairs(k, 2), 2))^2 + ...
             (p(3) - selectedSatPositions(pairs(k, 2), 3))^2) - ...
              radii_differences(k)), 1:size(pairs, 1))));

%% Define bounds for ECEF coordinates (in meters)
earth_radius = 6.371e6; % Approximate Earth radius in meters
alt_min = 0; % Minimum altitude (e.g., Dead Sea, below sea level)
alt_max = 1e3; % Maximum altitude (e.g., low Earth orbit)

% Calculate bounds for ECEF coordinates
lower_bound = [(earth_radius + alt_min) * -1, ...
               (earth_radius + alt_min) * -1, ...
               (earth_radius + alt_min) * -1];

upper_bound = [(earth_radius + alt_max), ...
               (earth_radius + alt_max), ...
               (earth_radius + alt_max)];

%% Set Eq. Solver
options = optimoptions('lsqnonlin', 'Display', 'iter');
% Solve for UE position using lsqnonlin with bounds
estimated_UE_position = lsqnonlin(objective_func, initial_guess, lower_bound, upper_bound, options);

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
% Colormap Selection
cmap = lines(size(pairs, 1)); % Generate a set of distinct colors based on 'lines' colormap.

% Colormap for satellites
sat_cmap = lines(size(selectedSatPositions, 1)); % Generate a set of distinct colors for satellites

figure;
hold on;
title('TDOA Localization: 3D Hyperboloids');
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
grid on;
axis equal;

% Number of grid points for the 3D space
n_points = 100;
[x_range, y_range, z_range] = meshgrid(linspace(-12e6, 12e6, n_points), ...
                                       linspace(-12e6, 12e6, n_points), ...
                                       linspace(0, 12e6, n_points));

% Legend entries for TDOA_ij pairs
legend_entries_3D = cell(size(pairs, 1) + size(selectedSatPositions, 1) + 2, 1);

% Plot hyperboloids for each TDOA
for k = 1:size(pairs, 1)
    sat1 = selectedSatPositions(pairs(k, 1), :);  % First satellite in the pair
    sat2 = selectedSatPositions(pairs(k, 2), :);  % Second satellite in the pair

    % Calculate distances for each point in 3D space
    d1 = sqrt((x_range - sat1(1)).^2 + (y_range - sat1(2)).^2 + (z_range - sat1(3)).^2);
    d2 = sqrt((x_range - sat2(1)).^2 + (y_range - sat2(2)).^2 + (z_range - sat2(3)).^2);

    % Compute the hyperboloid for the current TDOA
    hyperboloid = abs(d1 - d2) - abs(radii_differences(k));

    % Visualize the hyperboloid as an isosurface
    iso_val = 0;  % The zero level of the hyperboloid equation
    p = patch(isosurface(x_range, y_range, z_range, hyperboloid, iso_val));
    set(p, 'FaceColor', cmap(k, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Add legend entry for the TDOA pair
    legend_entries_3D{k} = sprintf('TDOA_{%d,%d}', pairs(k, 1), pairs(k, 2));
end

% Plot satellite positions
for i = 1:size(selectedSatPositions, 1)
    scatter3(selectedSatPositions(i, 1), selectedSatPositions(i, 2), selectedSatPositions(i, 3), ...
        100, sat_cmap(i, :), 'filled');
    legend_entries_3D{size(pairs, 1) + i} = sprintf('Satellite %d', i);
end

% Plot estimated UE position
scatter3(estimated_UE_position(1), estimated_UE_position(2), estimated_UE_position(3), ...
    200, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
legend_entries_3D{end-1} = 'Estimated UE Position';

% Plot actual UE position
scatter3(actual_UE_position(1), actual_UE_position(2), actual_UE_position(3), ...
    200, 'o', 'MarkerEdgeColor', 'g', 'LineWidth', 2);
legend_entries_3D{end} = 'Actual UE Position';

% Add lighting and view adjustments
camlight;
lighting gouraud;
view(3);  % 3D view

% Add legend
legend(legend_entries_3D, 'Location', 'bestoutside');

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
[x_range, y_range] = meshgrid(linspace(-12e6, 12e6, n_points), linspace(-12e6, 12e6, n_points));

% Assume a constant Z-coordinate (e.g., mean of satellite heights)
z_constant = ueStationECEF(3);

% Legend entries for TDOA_ij pairs
legend_entries_2D = cell(size(pairs, 1) + size(selectedSatPositions, 1) + 2, 1);

% Plot hyperbolas for each TDOA
for k = 1:size(pairs, 1)
    sat1 = selectedSatPositions(pairs(k, 1), :);  % First satellite in the pair
    sat2 = selectedSatPositions(pairs(k, 2), :);  % Second satellite in the pair

    % Calculate distances
    d1 = sqrt((x_range - sat1(1)).^2 + (y_range - sat1(2)).^2 + (z_constant - sat1(3)).^2);
    d2 = sqrt((x_range - sat2(1)).^2 + (y_range - sat2(2)).^2 + (z_constant - sat2(3)).^2);

    % Define the hyperbola
    hyperbola = abs(d1 - d2) - abs(radii_differences(k));

    % Plot the hyperbola (contour where hyperbola = 0)
    contour(x_range, y_range, hyperbola, [0, 0], 'Color', cmap(k, :), 'LineWidth', 1);

    % Add legend entry for the TDOA pair
    legend_entries_2D{k} = sprintf('TDOA_{%d,%d}', pairs(k, 1), pairs(k, 2));
end

% Plot satellite positions in 2D
for i = 1:size(selectedSatPositions, 1)
    scatter(selectedSatPositions(i, 1), selectedSatPositions(i, 2), 100, sat_cmap(i, :), 'filled');
    legend_entries_2D{size(pairs, 1) + i} = sprintf('Satellite %d', i);
end

% Plot estimated UE position in 2D
scatter(estimated_UE_position(1), estimated_UE_position(2), 200, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
legend_entries_2D{end-1} = 'Estimated UE Position';

% Plot actual UE position in 2D
scatter(actual_UE_position(1), actual_UE_position(2), 200, 'o', 'MarkerEdgeColor', 'g', 'LineWidth', 2);
legend_entries_2D{end} = 'Actual UE Position';

legend(legend_entries_2D, 'Location', 'bestoutside');
hold off;


function gdop = calculateGDOP(satPositions, uePosition)
    % Compute unit vectors from UE to satellites
    numSats = size(satPositions, 1);
    H = zeros(numSats, 4);
    for i = 1:numSats
        vec = satPositions(i, :) - uePosition;
        range = norm(vec);
        H(i, 1:3) = vec / range;
        H(i, 4) = 1; % Clock bias component
    end
    % GDOP is derived from the trace of (H'*H)^-1
    Q = inv(H' * H);
    gdop = sqrt(trace(Q));
end
