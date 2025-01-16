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

ueStationLLA = [40.786648, 29.449502, 182];
ueStationECEF = lla2ecef(ueStationLLA);

gsUE = groundStation(satscene, ...
                     "Latitude",  ueStationLLA(1), ...
                     "Longitude", ueStationLLA(2), ...
                     "Altitude",  ueStationLLA(3));

c = physconst("LightSpeed");
%% Find the access intervals
ac = access(constellation,gsUE);
accessIntervals = accessIntervals(ac);

totalSamples = seconds(stopTime - startTime) / sampleTime;

randomSampleIndex = randi([0, totalSamples - 1]);

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

%% Any satellite can be selected from the viewer 

tleStruct = tleread('leoSatelliteConstellation.tle');

satelliteNamesInTLE = {tleStruct.Name}';
indicesInTLE = find(matches(string(satelliteNamesInTLE), accessedSatellites));

orbitTime = randomDateTime;

accessedTLEStruct = tleStruct(indicesInTLE);

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

TOAs = accessedSatDistances / c ;

%% Initial guess for UE position
% initial_guess = mean(satPosxyz, 1);
initialGuess = lla2ecef([39.284593, 33.421097, 887]);
% initial_guess = [initial_guess, 0];
% initial_guess = [0,0,0];
%%

% Iterate over all subsets of satellite positions
numSats = size(accessedSatPositions, 1);
allCombinations = nchoosek(1:numSats, 4); % For 4 satellites at a time

for i = 1:size(allCombinations, 1)
    subset = allCombinations(i, :);
    subsetPositions = accessedSatPositions(subset, :);
    [gdop(i), ~] = calculateGDOP(subsetPositions, initialGuess(1:3));
end

[~, index] = min(gdop);
selectedSatIndices = allCombinations(index, :);
selectedSatPositions = accessedSatPositions(selectedSatIndices, :);


selectedTOAs = TOAs(selectedSatIndices);

pairs = nchoosek(1:length(selectedTOAs), 2);

TDOAs = arrayfun(@(row) selectedTOAs(pairs(row, 1)) - selectedTOAs(pairs(row, 2)), ...
                 1:size(pairs, 1));

%% Adding TDOA clock error terms
meanTOAClockError = 0; 
stdTOAClockError = 1e-6;
varTOAClockError = stdTOAClockError^2;
varTOADifferenceClockError = 2 * varTOAClockError;
stdTDOAError = sqrt(varTOADifferenceClockError);

% Generate random clock error samples
clockErrors = normrnd(meanTOAClockError, stdTDOAError, size(TDOAs));
TDOAswithError = TDOAs + clockErrors;

% Compute radii differences for TDOAs
radiiDifferenceGroundTruth = TDOAs * c; 
radiiDifferenceswithError = TDOAswithError * c;  % Convert TDOAs to distances

% True UE position (replace with actual coordinates)
actualUEPosition = ueStationECEF;  % Replace this variable with the actual UE position

%% Define bounds for ECEF coordinates (in meters)
earthRadius = 6.371e6; % Approximate Earth radius in meters
altitudeMin = 0; % Minimum altitude
altitudeMax = 1e3; % Maximum altitude
lowerBound = [(earthRadius + altitudeMin) * -1, ...
               (earthRadius + altitudeMin) * -1, ...
               (earthRadius + altitudeMin) * -1];
upperBound = [(earthRadius + altitudeMax), ...
               (earthRadius + altitudeMax), ...
               (earthRadius + altitudeMax)];

%% Compute Covariance Matrices for TDOA pairs
numPairs = size(pairs, 1);
covarianceMatrices = zeros(3, 3, numPairs);

for k = 1:numPairs
    % Placeholder: Use GDOP or simulation to generate realistic covariances
    covarianceMatrices(:, :, k) = (stdTDOAError^2) * eye(3);
end

%% Define the Cost Function
costFunction = @(xUE) sum(arrayfun(@(k) ...
    (1 / trace(covarianceMatrices(:, :, k))) * ...
    (norm(xUE - selectedSatPositions(pairs(k, 1), :)) - ...
     norm(xUE - selectedSatPositions(pairs(k, 2), :)) - ...
     radiiDifferenceswithError(k))^2, ...
    1:numPairs));
%% Non Linear Least Square Solution
ObjectiveEstimationErrorTDOA = @(p) arrayfun(@(k) ...
    abs(sqrt((p(1) - selectedSatPositions(pairs(k, 1), 1))^2 + ...
             (p(2) - selectedSatPositions(pairs(k, 1), 2))^2 + ...
             (p(3) - selectedSatPositions(pairs(k, 1), 3))^2) - ...
        sqrt((p(1) - selectedSatPositions(pairs(k, 2), 1))^2 + ...
             (p(2) - selectedSatPositions(pairs(k, 2), 2))^2 + ...
             (p(3) - selectedSatPositions(pairs(k, 2), 3))^2) - ...
        radiiDifferenceswithError(k)), 1:size(pairs, 1));

options = optimoptions('lsqnonlin', 'Display', 'iter');

estUEPosTDOAError = lsqnonlin(ObjectiveEstimationErrorTDOA, ...
                                    initialGuess(1:3), ...
                                    lowerBound, ...
                                    upperBound, ...
                                    options);
localizationErrorNLS = norm(estUEPosTDOAError(1:3) - actualUEPosition);


%% Solve the Optimization Problem
% initialGuess = mean(selectedSatPositions, 1); % Use centroid as an initial guess

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

[optimalPosition, fval, exitFlag] = fmincon(costFunction, initialGuess(1:3), [], [], [], [], ...
                                            lowerBound, upperBound, [], options);

% Display Results
disp('Optimal User Position (ECEF):');
disp(optimalPosition);

% Compute Error Relative to True Position
if exist('actualUEPosition', 'var')
    localizationError = norm(optimalPosition - actualUEPosition);
    disp('Localization Error (meters):');
    disp(localizationError);
end



% TWLS
[p_est, residual] = taylor_wls_tdoa(selectedSatPositions, TDOAswithError, pairs, covarianceMatrices, ...
                                    initialGuess, 10000, 1e-9, c);

% Display results
disp('NonLin LS Localization Error (meters):');
disp(localizationErrorNLS);

localizationErrorTWLS = norm(p_est - actualUEPosition);
    disp('TWLS Localization Error (meters):');
    disp(localizationErrorTWLS);

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

function [p_est, residual] = taylor_wls_tdoa(selectedSatPositions, TDOAs, pairs, covarianceMatrices, initialGuess, maxIter, tol, c)
    % Implements Taylor Weighted Least Squares (TWLS) for TDOA localization
    %
    % Inputs:
    %   selectedSatPositions - 4x3 matrix of satellite positions
    %   TDOAs - 6x1 vector of TDOA measurements
    %   pairs - 6x2 matrix of satellite pairs
    %   covarianceMatrices - 3x3x6 tensor of covariance matrices for each pair
    %   initialGuess - Initial UE position [x, y, z]
    %   maxIter - Maximum number of iterations
    %   tol - Convergence tolerance
    %   c - Speed of light
    %
    % Outputs:
    %   p_est - Estimated UE position [x, y, z]
    %   residual - Final residual norm

    % Initialize variables
    p_est = initialGuess;
    residual = inf; % Initialize residual norm
    iter = 0; % Iteration counter
    numPairs = size(pairs, 1);

    % Main iterative loop
    while iter < maxIter && residual > tol
        iter = iter + 1;

        % Initialize Jacobian matrix and residual vector
        H = zeros(numPairs, 3);
        delta = zeros(numPairs, 1);

        % Compute Jacobian and residuals for each pair
        for k = 1:numPairs
            sat_i = selectedSatPositions(pairs(k, 1), :);
            sat_j = selectedSatPositions(pairs(k, 2), :);

            % Distances from current estimate to satellites
            d_i = norm(p_est - sat_i);
            d_j = norm(p_est - sat_j);

            % Jacobian row for pair (i, j)
            H(k, :) = (p_est - sat_i) / d_i - (p_est - sat_j) / d_j;

            % Residual for pair (i, j)
            delta(k) = (d_i - d_j) - c * TDOAs(k);
        end

        % Construct the weight matrix
        W = zeros(numPairs, numPairs);
        for k = 1:numPairs
            W(k, k) = 1 / trace(covarianceMatrices(:, :, k)); % Use covariance trace as weight
        end

        % Solve for position update using Weighted Least Squares
        delta_p = (H' * W * H) \ (H' * W * delta);

        % Update position estimate
        p_est = p_est - delta_p';

        % Update residual norm
        residual = norm(delta);

        % Check convergence
        if residual < tol
            break;
        end
    end

    % Convergence check
    if residual > tol
        warning('TWLS did not converge within the maximum number of iterations.');
    end
end
