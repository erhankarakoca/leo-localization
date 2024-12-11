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
% initialGuess = [lla2ecef([39.284593, 33.421097, 887])
initialGuess = [lla2ecef([39.284593, 33.421097, 887]), 0];
% initial_guess = [initial_guess, 0];
% initial_guess = [0,0,0];
%%
% Define GDOP threshold (tune as needed)
gdopThreshold = 5;

% Initialize storage for selected subsets
% selectedCombinations = [];

% Iterate over all subsets of satellite positions
numSatsTotal = size(accessedSatPositions, 1);
numUsedSats = 4 ; 
allCombinations = nchoosek(1:numSatsTotal, numUsedSats); % For 4 satellites at a time

for i = 1:size(allCombinations, 1)
    subset = allCombinations(i, :);
    subsetPositions = accessedSatPositions(subset, :);
    gdop(i) = calculateGDOP(subsetPositions, initialGuess(1:3));
    % if gdop(i) <= gdopThreshold
    %     selectedCombinations = [selectedCombinations; subset]; %#ok<AGROW>
    % end
end

% Use selected combinations for TDOA calculations
% disp('Selected satellite combinations based on GDOP:');
% disp(selectedCombinations);

%% Descending order to choose different combinations
[~, GDOPAscendedIndexes]=sort(gdop, "ascend");
AscendedMultipleSatIndices = allCombinations(GDOPAscendedIndexes, : );

NumberofCombinations = 3; 

SelectedAscendedMultipleSatsIndicesMat = AscendedMultipleSatIndices(1:NumberofCombinations,:);
SelectedAscendedMultipleSatsIndices = reshape(SelectedAscendedMultipleSatsIndicesMat.',1,[]);
SelectedAscendedMultipleSatPositions = accessedSatPositions(SelectedAscendedMultipleSatsIndices,:);

% Compute TOAs for selected satellites
selectedTOAs = TOAs(SelectedAscendedMultipleSatsIndices);

TDOAs = zeros(6,NumberofCombinations);

% Define TDOA clock error terms
meanTOAClockError = 0; 
stdTOAClockError = 1e-6;
varTOAClockError = stdTOAClockError^2;

varTOADifferenceClockError = 2*varTOAClockError;

stdTDOAError = sqrt(varTOADifferenceClockError);

% Generate random clock error samples
clockErrors = normrnd(meanTOAClockError, stdTDOAError, size(TDOAs));



%% Define bounds for ECEF coordinates (in meters)
earthRadius = 6.371e6; % Approximate Earth radius in meters
altitudeMin = 0; % Minimum altitude (e.g., Dead Sea, below sea level)
altitudeMax = 1e3; % Maximum altitude (e.g., low Earth orbit)

% Calculate bounds for ECEF coordinates
lowerBound = [(earthRadius + altitudeMin) * -1, ...
               (earthRadius + altitudeMin) * -1, ...
               (earthRadius + altitudeMin) * -1, ...
                    ];

upperBound = [(earthRadius + altitudeMax), ...
               (earthRadius + altitudeMax), ...
               (earthRadius + altitudeMax), ...
              ];

options = optimoptions('lsqnonlin', 'Display', 'iter');

actualUEPosition = ueStationECEF;  


%% 
for i = 1 : NumberofCombinations
    
    pairs(:,:,i) = nchoosek(1:length(selectedTOAs(4*(i-1)+1:i*4)), 2) + 4*(i-1);
    currentPairs = pairs(:, :, i);
    
    currentTDOAs = arrayfun(@(row) selectedTOAs(currentPairs(row, 1)) - ...
                                      selectedTOAs(currentPairs(row, 2)), ...
                            1:size(currentPairs, 1))';
    
    TDOAs(:,i) = currentTDOAs;
    TDOAswithError(:,i) = currentTDOAs + clockErrors(:,i);

    % Compute radii differences for TDOAs
    radiiDifferenceGroundTruth = TDOAs(:,i) * c; 
    radiiDifferenceswithError = TDOAswithError(:,i) * c;  % Convert TDOAs to distances


    %% Objective Functions
    ObjectiveEstimationErrorTDOA = @(p) arrayfun(@(k) ...
        abs(sqrt((p(1) - SelectedAscendedMultipleSatPositions(currentPairs(k, 1), 1))^2 + ...
                 (p(2) - SelectedAscendedMultipleSatPositions(currentPairs(k, 1), 2))^2 + ...
                 (p(3) - SelectedAscendedMultipleSatPositions(currentPairs(k, 1), 3))^2) - ...
            sqrt((p(1) - SelectedAscendedMultipleSatPositions(currentPairs(k, 2), 1))^2 + ...
                 (p(2) - SelectedAscendedMultipleSatPositions(currentPairs(k, 2), 2))^2 + ...
                 (p(3) - SelectedAscendedMultipleSatPositions(currentPairs(k, 2), 3))^2) - ...
            radiiDifferenceswithError(k)), 1:size(currentPairs, 1));
    
    ObjectiveGroundTruthTDOA = @(p) arrayfun(@(k) ...
        abs(sqrt((p(1) - SelectedAscendedMultipleSatPositions(currentPairs(k, 1), 1))^2 + ...
                 (p(2) - SelectedAscendedMultipleSatPositions(currentPairs(k, 1), 2))^2 + ...
                 (p(3) - SelectedAscendedMultipleSatPositions(currentPairs(k, 1), 3))^2) - ...
            sqrt((p(1) - SelectedAscendedMultipleSatPositions(currentPairs(k, 2), 1))^2 + ...
                 (p(2) - SelectedAscendedMultipleSatPositions(currentPairs(k, 2), 2))^2 + ...
                 (p(3) - SelectedAscendedMultipleSatPositions(currentPairs(k, 2), 3))^2) - ...
            radiiDifferenceGroundTruth(k)), 1:size(pairs, 1));
    
    
    % Solve using lsqnonlin for the current combination
    estUEPosTDOAError = lsqnonlin(ObjectiveEstimationErrorTDOA, ...
                                        initialGuess(1:3), ...
                                        lowerBound(1:3), ...
                                        upperBound(1:3), ...
                                        options);

    estUEPosGroundTruthTDOA = lsqnonlin(ObjectiveGroundTruthTDOA, ...
                                                        initialGuess(1:3), ...
                                                        lowerBound(1:3), ...
                                                        upperBound(1:3), ...
                                                        options);
    
    
    % Calculate localization errors for the current combination
    localizationErrorsWithError(i) = norm(estUEPosTDOAError(1:3) - actualUEPosition);
    localizationErrorsGroundTruth(i) = norm(estUEPosGroundTruthTDOA(1:3) - actualUEPosition);

    % Store the results
    estimatedPositionsWithError(i, :) = estUEPosTDOAError;
    estimatedPositionsGroundTruth(i, :) = estUEPosGroundTruthTDOA;


   
end


% Display results for all combinations
for i = 1:NumberofCombinations
    fprintf('Combination %d:\n', i);
    fprintf('Estimated UE Position with Ground Truth TDOA Result: [%.3f, %.3f, %.3f]\n', ...
            estimatedPositionsGroundTruth(i, :));
    fprintf('Estimated UE Position with TDOA Error Result: [%.3f, %.3f, %.3f]\n', ...
            estimatedPositionsWithError(i, :));
    fprintf('Localization Error with Ground Truth TDOA (meters): %.3f\n', ...
            localizationErrorsGroundTruth(i));
    fprintf('Localization Error with TDOA Error (meters): %.3f\n', ...
            localizationErrorsWithError(i));
end

disp('Actual UE Position:');
disp(actualUEPosition);


%% 2D Visualization of TDOA Localization (Hyperbolas)

% % Colormap Selection
% cmap = lines(size(pairs, 1)); % Generate a set of distinct colors based on 'lines' colormap.
% 
% % Colormap for satellites
% sat_cmap = lines(size(selectedSatPositions, 1)); % Generate a set of distinct colors for satellites
% 
% figure;
% hold on;
% title('TDOA Localization with 2D Hyperbolas');
% xlabel('X (meters)');
% ylabel('Y (meters)');
% grid on;
% axis equal;
% 
% % Number of grid points for visualization
% n_points = 1000;
% [x_range, y_range] = meshgrid(linspace(-12e6, 12e6, n_points), linspace(-12e6, 12e6, n_points));
% 
% % Assume a constant Z-coordinate (e.g., mean of satellite heights)
% z_constant = ueStationECEF(3);
% 
% % Legend entries for TDOA_ij pairs
% legend_entries_2D = cell(size(pairs, 1) + size(selectedSatPositions, 1) + 2, 1);
% 
% % Plot hyperbolas for each TDOA
% for k = 1:size(pairs, 1)
%     sat1 = selectedSatPositions(pairs(k, 1), :);  % First satellite in the pair
%     sat2 = selectedSatPositions(pairs(k, 2), :);  % Second satellite in the pair
% 
%     % Calculate distances
%     d1 = sqrt((x_range - sat1(1)).^2 + (y_range - sat1(2)).^2 + (z_constant - sat1(3)).^2);
%     d2 = sqrt((x_range - sat2(1)).^2 + (y_range - sat2(2)).^2 + (z_constant - sat2(3)).^2);
% 
%     % Define the hyperbola
%     hyperbola = abs(d1 - d2) - abs(radiiDifferenceswithError(k));
% 
%     % Plot the hyperbola (contour where hyperbola = 0)
%     contour(x_range, y_range, hyperbola, [0, 0], 'Color', cmap(k, :), 'LineWidth', 1);
% 
%     % Add legend entry for the TDOA pair
%     legend_entries_2D{k} = sprintf('TDOA_{%d,%d}', pairs(k, 1), pairs(k, 2));
% end
% 
% % Plot satellite positions in 2D
% for i = 1:size(selectedSatPositions, 1)
%     scatter(selectedSatPositions(i, 1), selectedSatPositions(i, 2), 100, sat_cmap(i, :), 'filled');
%     legend_entries_2D{size(pairs, 1) + i} = sprintf('Satellite %d', i);
% end
% 
% % Plot estimated UE position in 2D
% scatter(estUEPosTDOAError(1), estUEPosTDOAError(2), 200, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
% legend_entries_2D{end-1} = 'Estimated UE Position';
% 
% % Plot actual UE position in 2D
% scatter(actualUEPosition(1), actualUEPosition(2), 200, 'o', 'MarkerEdgeColor', 'g', 'LineWidth', 2);
% legend_entries_2D{end} = 'Actual UE Position';
% 
% legend(legend_entries_2D, 'Location', 'bestoutside');
% hold off;


%%

% TDOA error standard deviation (sqrt(2) * sigma_toa_dist)
stdTDOADistError = stdTDOAError * c;

% Confidence level (95%)
confidence_scale = 2; % ~2 sigma for 95% confidence

% Generate uncertainty ellipsoids for each combination
figure;
hold on;

for comb = 1:NumberofCombinations
    % Current estimated position
    currentEstPos = estimatedPositionsWithError(comb, :);
    
    % Define covariance matrix for TDOA-based uncertainty
    % Assuming isotropic uncertainty for simplicity
    cov_matrix = stdTDOADistError^2 * eye(3);
    
    % Generate ellipsoid centered at the estimated UE position
    [X, Y, Z] = ellipsoid(currentEstPos(1), currentEstPos(2), currentEstPos(3), ...
                          confidence_scale * sqrt(cov_matrix(1, 1)), ...
                          confidence_scale * sqrt(cov_matrix(2, 2)), ...
                          confidence_scale * sqrt(cov_matrix(3, 3)), 50);
                      
    % Plot ellipsoid
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
         'DisplayName', sprintf('Uncertainty Area (Combination %d)', comb));
end

% Plot estimated UE positions
for comb = 1:NumberofCombinations
    plot3(estimatedPositionsWithError(comb, 1), ...
          estimatedPositionsWithError(comb, 2), ...
          estimatedPositionsWithError(comb, 3), ...
          'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
          'DisplayName', sprintf('Estimated Position (Combination %d)', comb));
end

% Plot actual UE position
plot3(ueStationECEF(1), ueStationECEF(2), ueStationECEF(3), ...
      'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', ...
      'DisplayName', 'Actual Position');

xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Uncertainty Regions with Estimated and Actual UE Positions');
legend show;
grid on;
view(3);
hold off;



%% Helper Functions
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
