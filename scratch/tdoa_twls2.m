clear all;
close all;
clc;

%% Satellite Scenario Setup
startTime = datetime(2020, 05, 04, 18, 45, 50);
stopTime = datetime(2020, 05, 04, 19, 02, 20);
sampleTime = 10;
satscene = satelliteScenario(startTime, stopTime, sampleTime);

% Add satellites from TLE file
tleFile = "leoSatelliteConstellation.tle";
constellation = satellite(satscene, tleFile);

% Define ground station (UE)
ueStationLLA = [40.786648, 29.449502, 182]; % [Lat, Long, Alt]
ueStationECEF = lla2ecef(ueStationLLA); % Convert to ECEF
gsUE = groundStation(satscene, ...
    "Latitude", ueStationLLA(1), ...
    "Longitude", ueStationLLA(2), ...
    "Altitude", ueStationLLA(3));

c = physconst("LightSpeed");

%% Access Intervals and Random Time Selection
ac = access(constellation, gsUE);
accessIntervals = accessIntervals(ac);

% Random sample date-time within the simulation period
totalSamples = seconds(stopTime - startTime) / sampleTime;
randomSampleIndex = randi([0, totalSamples - 1]);
randomDateTime = startTime + seconds(randomSampleIndex * sampleTime);
randomDateTime = datetime(randomDateTime, 'TimeZone', 'UTC');

% Get accessed satellites at the selected time
accessedSatellites = [];
for i = 1:height(accessIntervals)
    accessStartTime = accessIntervals{i, 4};
    accessEndTime = accessIntervals{i, 5};
    if randomDateTime >= accessStartTime && randomDateTime <= accessEndTime
        accessedSatellites = [accessedSatellites; accessIntervals{i, 1}];
    end
end

if isempty(accessedSatellites)
    disp('No satellites accessed at the random date-time.');
    return;
else
    fprintf('Random Date-Time: %s\n', datestr(randomDateTime));
    disp('Satellites accessed at this time:');
    disp(accessedSatellites);
end

% Propagate satellite positions
tleStruct = tleread(tleFile);
satelliteNamesInTLE = {tleStruct.Name}';
indicesInTLE = find(matches(string(satelliteNamesInTLE), accessedSatellites));
accessedTLEStruct = tleStruct(indicesInTLE);
[accessedSatPositions, ~] = propagateOrbit(randomDateTime, ...
    accessedTLEStruct, ...
    "OutputCoordinateFrame", "fixed-frame");
accessedSatPositions = squeeze(accessedSatPositions)';

% Compute TOAs
accessedSatDistances = vecnorm(accessedSatPositions - ueStationECEF, 2, 2);
TOAs = accessedSatDistances / c;

%% Subset Selection and TWLS Localization
numSats = size(accessedSatPositions, 1);
if numSats < 4
    disp('Not enough satellites for localization.');
    return;
end

% Select subsets of 4 satellites
allCombinations = nchoosek(1:numSats, 4);
estimatedPositions = zeros(size(allCombinations, 1), 3);
localizationErrors = zeros(size(allCombinations, 1), 1);

% TWLS for each subset
for i = 1:size(allCombinations, 1)
    subset = allCombinations(i, :);
    subsetPositions = accessedSatPositions(subset, :);
    subsetTOAs = TOAs(subset);

    % Compute TDOAs relative to the first satellite in the subset
    referenceTOA = subsetTOAs(1);
    tdoaMeasurements = subsetTOAs - referenceTOA;
    tdoaMeasurements = tdoaMeasurements(2:end); % Remove reference satellite

    % Compute GDOP and construct weights
    gdopValues = calculateGDOP(subsetPositions, ueStationECEF);
    weights = 1 ./ (gdopValues.^2); % Inverse GDOP squared
    weightMatrix = diag(weights(2:end)); % Exclude reference satellite

    % Ensure dimensions match
    fprintf('Subset %d: tdoaMeasurements size = [%d, %d], weightMatrix size = [%d, %d]\\n', ...
        i, size(tdoaMeasurements), size(weightMatrix));

    % TWLS Localization
    [estPosition, ~] = taylor_wls_tdoa(subsetPositions, tdoaMeasurements / c, ...
        ueStationECEF, 100, 1e-6, c, weightMatrix);

    % Store results
    estimatedPositions(i, :) = estPosition;
    localizationErrors(i) = norm(estPosition - ueStationECEF);
end


% Find best estimate (minimum localization error)
[~, bestIndex] = min(localizationErrors);
bestEstimate = estimatedPositions(bestIndex, :);

%% Display Results
disp('Best Estimate of User Position (ECEF):');
disp(bestEstimate);
disp('True Position (ECEF):');
disp(ueStationECEF);
disp('Localization Error (meters):');
disp(localizationErrors(bestIndex));

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

function [p_est, residual] = taylor_wls_tdoa(satPositions, tdoaMeasurements, initialGuess, maxIter, tol, c, W)
    % Implements Taylor WLS for TDOA-based localization with GDOP-based weights
    %
    % Inputs:
    %   satPositions - Nx3 matrix of satellite positions
    %   tdoaMeasurements - Measured TDOAs (N-1x1 vector)
    %   initialGuess - Initial UE position [x, y, z]
    %   maxIter - Maximum number of iterations
    %   tol - Convergence tolerance
    %   c - Speed of light
    %   W - Weight matrix (N-1xN-1 diagonal matrix)
    %
    % Outputs:
    %   p_est - Estimated UE position [x, y, z]
    %   residual - Final residual norm

    % Initialize variables
    p_est = initialGuess;
    numSats = size(satPositions, 1);
    residual = inf; % Initialize residual norm
    iter = 0; % Iteration counter

    % Ensure dimensions match
    if length(tdoaMeasurements) ~= numSats - 1
        error('tdoaMeasurements must have size of (numSats - 1).');
    end
    if size(W, 1) ~= numSats - 1 || size(W, 2) ~= numSats - 1
        error('Weight matrix W must be (numSats - 1)x(numSats - 1).');
    end

    % Main iteration loop
    while iter < maxIter && residual > tol
        iter = iter + 1;

        % Compute distances from current estimate
        distances = vecnorm(satPositions - p_est, 2, 2);

        % Linearize the TDOA equations
        H = zeros(numSats - 1, 3); % Jacobian matrix
        delta = zeros(numSats - 1, 1); % Residual vector

        for i = 2:numSats
            % Gradient of TDOA equation
            H(i-1, :) = (p_est - satPositions(i, :)) / distances(i) - ...
                        (p_est - satPositions(1, :)) / distances(1);

            % TDOA residual
            delta(i-1) = (distances(i) - distances(1)) - c * tdoaMeasurements(i-1);
        end

        % Solve for position update using Weighted Least Squares
        delta_p = (H' * W * H) \ (H' * W * delta);

        % Update position estimate
        p_est = p_est - delta_p';

        % Update residual norm
        residual = norm(delta);
    end

    % Convergence check
    if residual > tol
        warning('Taylor WLS did not converge within the maximum number of iterations.');
    end
end

