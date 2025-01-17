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

NumberofCombinations = 10; 

SatCombinationsMat = AscendedMultipleSatIndices(1:NumberofCombinations,:);
SatCombinationsIndices = reshape(SatCombinationsMat.',1,[]);
SatPositionsCombinationsCombined = accessedSatPositions(SatCombinationsIndices,:);

% Compute TOAs for selected satellites
selectedTOAs = TOAs(SatCombinationsIndices);

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
        abs(sqrt((p(1) - SatPositionsCombinationsCombined(currentPairs(k, 1), 1))^2 + ...
                 (p(2) - SatPositionsCombinationsCombined(currentPairs(k, 1), 2))^2 + ...
                 (p(3) - SatPositionsCombinationsCombined(currentPairs(k, 1), 3))^2) - ...
            sqrt((p(1) - SatPositionsCombinationsCombined(currentPairs(k, 2), 1))^2 + ...
                 (p(2) - SatPositionsCombinationsCombined(currentPairs(k, 2), 2))^2 + ...
                 (p(3) - SatPositionsCombinationsCombined(currentPairs(k, 2), 3))^2) - ...
            radiiDifferenceswithError(k)), 1:size(currentPairs, 1));
    
    ObjectiveGroundTruthTDOA = @(p) arrayfun(@(k) ...
        abs(sqrt((p(1) - SatPositionsCombinationsCombined(currentPairs(k, 1), 1))^2 + ...
                 (p(2) - SatPositionsCombinationsCombined(currentPairs(k, 1), 2))^2 + ...
                 (p(3) - SatPositionsCombinationsCombined(currentPairs(k, 1), 3))^2) - ...
            sqrt((p(1) - SatPositionsCombinationsCombined(currentPairs(k, 2), 1))^2 + ...
                 (p(2) - SatPositionsCombinationsCombined(currentPairs(k, 2), 2))^2 + ...
                 (p(3) - SatPositionsCombinationsCombined(currentPairs(k, 2), 3))^2) - ...
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

%%

% TDOA error standard deviation (sqrt(2) * sigma_toa_dist)
stdTDOADistError = stdTDOAError * c;

% Confidence level (95%)
confidence_scale = 2; % ~2 sigma for 95% confidence

% Loop through each combination
% Plot uncertainty ellipsoid
figure;
for comb_idx = 1:NumberofCombinations
    % Get the current satellite combination
    selected_combination = SatCombinationsMat(comb_idx, :);
    
    % Extract the positions of the selected satellites
    selectedSatPositionsCombination = accessedSatPositions(selected_combination,:);
    
    % Reference satellite (use the first satellite in the combination)
    referenceSat = selectedSatPositionsCombination(1, :);
    
    % Initialize Jacobian matrix (TDOA rows x 3 columns for x, y, z)
    G = zeros(length(selected_combination) - 1, 3);
    
    % Compute Jacobian rows for each satellite pair in the current combination
    for i = 2:length(selected_combination)
        % Position vectors for the current satellite and reference satellite
        sat_i = selectedSatPositionsCombination(i, :);
        ref_sat = referenceSat;
    
        % Distance vectors from UE to satellites
        d_i = norm(sat_i - estimatedPositionsWithError(comb_idx, :)); % Distance to satellite i
        d_ref = norm(ref_sat - estimatedPositionsWithError(comb_idx, :)); % Distance to reference satellite
    
        % Jacobian row for satellite i
        G(i - 1, :) = (sat_i - estimatedPositionsWithError(comb_idx, :)) / d_i - (ref_sat - estimatedPositionsWithError(comb_idx, :)) / d_ref;
    end
    
    % Display the Jacobian matrix for the current combination
    % disp(['Jacobian matrix (G) for combination ', num2str(comb_idx), ':']);
    % disp(G);
    
    % TDOA covariance matrix
    covTDOADist = stdTDOADistError^2 * eye(size(G,1));
    
    % Compute position covariance matrix
    JacobianMat = inv(G) * covTDOADist * inv(G)'; % Transform TDOA errors to position errors
    
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
    
    
    surf(X_rot + estimatedPositionsWithError(comb_idx, 1), Y_rot + estimatedPositionsWithError(comb_idx, 2), Z_rot + estimatedPositionsWithError(comb_idx, 3), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none');
    hold on;
    
    % Plot estimated and actual positions
    scatter3(estimatedPositionsWithError(comb_idx, 1), estimatedPositionsWithError(comb_idx, 2), estimatedPositionsWithError(comb_idx, 3), ...
        100, 'r', 'filled', 'DisplayName', 'Estimated Position');
    scatter3(ueStationECEF(1), ueStationECEF(2), ueStationECEF(3), ...
        100, 'b', 'filled', 'DisplayName', 'Actual Position');
    
    title(['3D Position Uncertainty Region for Combination ', num2str(comb_idx)]);
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    legend('show');
    grid on;
    axis equal;
end


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
