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

rng(42)
%% Find the access intervals
ac = access(constellation,gsUE);
accessIntervals = accessIntervals(ac);

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
for i = 1:height(accessIntervals)
    % Extract start and end time from the table
    accessStartTime = accessIntervals{i, 4}; % 4th column: Access start date-time
    accessEndTime = accessIntervals{i, 5};   % 5th column: Access end date-time
    
    % Check if the random date-time is within the access interval
    if randomDateTime >= accessStartTime && randomDateTime <= accessEndTime
        % If true, append the satellite name (1st column) to the array
        accessedSatellites = [accessedSatellites; accessIntervals{i, 1}];
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

% Tune the data type as desired
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
initialGuess = [lla2ecef([39.284593, 33.421097, 887]), 0];
% initial_guess = [initial_guess, 0];
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
    [gdop(i), ~] = calculateGDOP(subsetPositions, initialGuess(1:3));
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
pairs= pairs(1:4,:);
% Calculate TDOA for each pair
TDOAs = arrayfun(@(row) selectedTOAs(pairs(row, 1)) - selectedTOAs(pairs(row, 2)), ...
                 1:size(pairs, 1));
%TDOAs= TDOAs(1:4);
%% Adding TDOA clock error terms
% Define clock error statistics (in seconds)
meanTOAClockError = 0; 
stdTOAClockError = 1e-6;
varTOAClockError = stdTOAClockError^2;

varTOADifferenceClockError = 2*varTOAClockError;

stdTDOAError = sqrt(varTOADifferenceClockError);

% Generate random clock error samples
clockErrors = normrnd(meanTOAClockError, stdTDOAError, size(TDOAs));

TDOAswithError = TDOAs + clockErrors;
%%
% True UE position (replace with actual coordinates)
actualUEPosition = ueStationECEF;  % Replace this variable with the actual UE position

% Compute radii differences for TDOAs
radiiDifferenceGroundTruth = TDOAs * c; 
radiiDifferenceswithError = TDOAswithError * c;  % Convert TDOAs to distances


%% Objective Functions
lambda = 1; % Regularization weight

% ObjectiveRegularizedClockErrorVariable = @(p) arrayfun(@(k) ...
%     abs(sqrt((p(1) - selectedSatPositions(pairs(k, 1), 1))^2 + ...
%              (p(2) - selectedSatPositions(pairs(k, 1), 2))^2 + ...
%              (p(3) - selectedSatPositions(pairs(k, 1), 3))^2) - ...
%         sqrt((p(1) - selectedSatPositions(pairs(k, 2), 1))^2 + ...
%              (p(2) - selectedSatPositions(pairs(k, 2), 2))^2 + ...
%              (p(3) - selectedSatPositions(pairs(k, 2), 3))^2) - ...
%         radiiDifferenceswithError(k) +  lambda *p(4)*c ) , 1:size(pairs, 1));

ObjectiveRegularizedClockErrorVariable = @(p) (arrayfun(@(k) ...
    abs((sqrt((p(1) - selectedSatPositions(pairs(k, 1), 1))^2 + ...
          (p(2) - selectedSatPositions(pairs(k, 1), 2))^2 + ...
          (p(3) - selectedSatPositions(pairs(k, 1), 3))^2) - ...
      sqrt((p(1) - selectedSatPositions(pairs(k, 2), 1))^2 + ...
           (p(2) - selectedSatPositions(pairs(k, 2), 2))^2 + ...
           (p(3) - selectedSatPositions(pairs(k, 2), 3))^2) - ...
       radiiDifferenceswithError(k))) / stdTDOAError, 1:size(pairs, 1)));

ObjectiveEstimationErrorTDOA = @(p) arrayfun(@(k) ...
    abs(sqrt((p(1) - selectedSatPositions(pairs(k, 1), 1))^2 + ...
             (p(2) - selectedSatPositions(pairs(k, 1), 2))^2 + ...
             (p(3) - selectedSatPositions(pairs(k, 1), 3))^2) - ...
        sqrt((p(1) - selectedSatPositions(pairs(k, 2), 1))^2 + ...
             (p(2) - selectedSatPositions(pairs(k, 2), 2))^2 + ...
             (p(3) - selectedSatPositions(pairs(k, 2), 3))^2) - ...
        radiiDifferenceswithError(k)), 1:size(pairs, 1));

ObjectiveGroundTruthTDOA = @(p) arrayfun(@(k) ...
    abs(sqrt((p(1) - selectedSatPositions(pairs(k, 1), 1))^2 + ...
             (p(2) - selectedSatPositions(pairs(k, 1), 2))^2 + ...
             (p(3) - selectedSatPositions(pairs(k, 1), 3))^2) - ...
        sqrt((p(1) - selectedSatPositions(pairs(k, 2), 1))^2 + ...
             (p(2) - selectedSatPositions(pairs(k, 2), 2))^2 + ...
             (p(3) - selectedSatPositions(pairs(k, 2), 3))^2) - ...
        radiiDifferenceGroundTruth(k)), 1:size(pairs, 1));


%% Define bounds for ECEF coordinates (in meters)
earthRadius = 6.371e6; % Approximate Earth radius in meters
altitudeMin = 0; % Minimum altitude (e.g., Dead Sea, below sea level)
altitudeMax = 1e3; % Maximum altitude (e.g., low Earth orbit)

% Calculate bounds for ECEF coordinates
lowerBound = [(earthRadius + altitudeMin) * -1, ...
               (earthRadius + altitudeMin) * -1, ...
               (earthRadius + altitudeMin) * -1];

upperBound = [(earthRadius + altitudeMax), ...
               (earthRadius + altitudeMax), ...
               (earthRadius + altitudeMax)];

%% Set Eq. Solver
options = optimoptions('lsqnonlin', 'Display', 'iter');
% Solve for UE position using lsqnonlin with bounds
estUEPosClockErrorVar = lsqnonlin(ObjectiveRegularizedClockErrorVariable, ...
                                                    initialGuess, ...
                                                    lowerBound, ...
                                                    upperBound, ...
                                                    options);
estUEPosTDOAError = lsqnonlin(ObjectiveEstimationErrorTDOA, ...
                                    initialGuess(1:3), ...
                                    lowerBound, ...
                                    upperBound, ...
                                    options);

estUEPosGroundTruthTDOA = lsqnonlin(ObjectiveGroundTruthTDOA, ...
                                                    initialGuess(1:3), ...
                                                    lowerBound, ...
                                                    upperBound, ...
                                                    options);
% Calculate error
localization_error_with_radii_error = norm(estUEPosTDOAError(1:3) - actualUEPosition);
localization_error_with_radii_true = norm(estUEPosGroundTruthTDOA(1:3) - actualUEPosition);
localization_error_with_clock_error = norm(estUEPosClockErrorVar(1:3) - actualUEPosition);

% Display results
disp('Estimated UE Position with Ground Truth TDOA Result:');
disp(estUEPosGroundTruthTDOA);

disp('Estimated UE Position with TDOA Error Result:');
disp(estUEPosTDOAError);

disp('Estimated UE Position with Clock Error Variable Result:');
disp(estUEPosClockErrorVar);

disp('Actual UE Position:');
disp(actualUEPosition);

disp('Localization Error with ground truth TDOA(meters):');
disp(localization_error_with_radii_true);

disp('Localization Error with TDOA error (meters):');
disp(localization_error_with_radii_error);

disp('Localization Error witch clock error variable introduced (meters):');
disp(localization_error_with_clock_error);

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
    hyperboloid = abs(d1 - d2) - abs(radiiDifferenceswithError(k));

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
scatter3(estUEPosTDOAError(1), estUEPosTDOAError(2), estUEPosTDOAError(3), ...
    200, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
legend_entries_3D{end-1} = 'Estimated UE Position';

% Plot actual UE position
scatter3(actualUEPosition(1), actualUEPosition(2), actualUEPosition(3), ...
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
    hyperbola = abs(d1 - d2) - abs(radiiDifferenceswithError(k));

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
scatter(estUEPosTDOAError(1), estUEPosTDOAError(2), 200, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
legend_entries_2D{end-1} = 'Estimated UE Position';

% Plot actual UE position in 2D
scatter(actualUEPosition(1), actualUEPosition(2), 200, 'o', 'MarkerEdgeColor', 'g', 'LineWidth', 2);
legend_entries_2D{end} = 'Actual UE Position';

legend(legend_entries_2D, 'Location', 'bestoutside');
hold off;

%% Add Estimation Uncertainty


% TDOA error standard deviation (sqrt(2) * sigma_toa_dist)
stdTDOADistError = stdTDOAError*c;

% Confidence level (95%)
confidence_scale = 2; % ~2 sigma for 95% confidence

% Generate uncertainty ellipsoids for each TDOA
figure;
hold on;
for i = 1:size(selectedSatPositions, 1)
    % Define covariance matrix for TDOA-based uncertainty
    % Assuming isotropic uncertainty for simplicity
    cov_matrix = stdTDOADistError^2 * eye(3);
    
    % Generate ellipsoid centered at the estimated UE position
    [X, Y, Z] = ellipsoid(estUEPosTDOAError(1), estUEPosTDOAError(2), estUEPosTDOAError(3), ...
                          confidence_scale * sqrt(cov_matrix(1,1)), ...
                          confidence_scale * sqrt(cov_matrix(2,2)), ...
                          confidence_scale * sqrt(cov_matrix(3,3)), 50);
                      
    % Plot ellipsoid
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Plot estimated UE position
plot3(estUEPosTDOAError(1), estUEPosTDOAError(2), estUEPosTDOAError(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Estimated Position');

% Plot actual UE position
plot3(ueStationECEF(1), ueStationECEF(2), ueStationECEF(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Actual Position');

xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Uncertainty Regions with Estimated and Actual UE Positions');
legend show;
grid on;
view(3);
hold off;

%% Jacobian computation for better uncertainty estimation 

referenceSat = selectedSatPositions(1, :);

% Initialize Jacobian matrix (TDOA rows x 3 columns for x, y, z)
G = zeros(4 - 1, 3);

% Compute Jacobian rows for each satellite pair
for i = 2:4
    % Position vectors for the current satellite and reference satellite
    sat_i = selectedSatPositions(i, :);
    ref_sat = referenceSat;

    % Distance vectors from UE to satellites
    d_i = norm(sat_i - initialGuess(1:3)); % Distance to satellite i
    d_ref = norm(ref_sat - initialGuess(1:3)); % Distance to reference satellite

    % Jacobian row for satellite i
    G(i - 1, :) = (sat_i - initialGuess(1:3)) / d_i - (ref_sat - initialGuess(1:3)) / d_ref;
end

% Display the Jacobian matrix
disp('Jacobian matrix (G):');
disp(G);

% TDOA covariance matrix
covTDOADist = stdTDOADistError^2 * eye(size(G,1));
% Compute position covariance matrix
JacobianMat = inv(G)*covTDOADist*inv(G)'; % Transform TDOA errors to position errors
% Eigenvalue decomposition of position covariance matrix
[U, S, ~] = svd(JacobianMat); % U: eigenvectors, S: eigenvalues

% Scale eigenvalues for 95% confidence
confidence_scale = 2; % 95% confidence interval
radii = confidence_scale * sqrt(diag(S)); % Radii of uncertainty ellipsoid

% Generate ellipsoid data
[X, Y, Z] = ellipsoid(0, 0, 0, radii(1), radii(2), radii(3), 50);

% Rotate ellipsoid to align with eigenvectors
ellipsoid_points = [X(:), Y(:), Z(:)] * U'; % Align with eigenvectors
X_rot = reshape(ellipsoid_points(:, 1), size(X));
Y_rot = reshape(ellipsoid_points(:, 2), size(Y));
Z_rot = reshape(ellipsoid_points(:, 3), size(Z));


% Plot uncertainty ellipsoid
figure;
surf(X_rot + estUEPosTDOAError(1), Y_rot + estUEPosTDOAError(2), Z_rot + estUEPosTDOAError(3), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;

% Plot estimated and actual positions
scatter3(estUEPosTDOAError(1), estUEPosTDOAError(2), estUEPosTDOAError(3), ...
    100, 'r', 'filled', 'DisplayName', 'Estimated Position');
scatter3(ueStationECEF(1), ueStationECEF(2), ueStationECEF(3), ...
    100, 'b', 'filled', 'DisplayName', 'Actual Position');

title('3D Position Uncertainty Region');
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('show');
grid on;
axis equal;

%% Helper Functions
function [gdop, Q] = calculateGDOP(satPositions, uePosition)
    % Compute unit vectors from UE to satellites
    numSats = size(satPositions, 1);
    H = zeros(numSats, 3);
    for i = 1:numSats
        vec = satPositions(i, :) - uePosition;
        range = norm(vec);
        H(i, 1:3) = vec / range;
        % H(i, 4) = 1; % Clock bias component
    end
    % GDOP is derived from the trace of (H'*H)^-1
    Q = inv(H' * H);
    gdop = sqrt(trace(Q));
end
