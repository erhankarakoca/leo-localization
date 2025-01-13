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

% Constants
c = physconst('LightSpeed'); % Speed of light in meters/second


% Compute TOAs for selected satellites
selectedTOAs = TOAs(selectedSatIndices);

% Recompute pairs for selected satellites
pairs = nchoosek(1:length(selectedTOAs), 2);

% Calculate TDOA for each pair
TDOAs = arrayfun(@(row) selectedTOAs(pairs(row, 1)) - selectedTOAs(pairs(row, 2)), ...
                 1:size(pairs, 1));

%% Adding TDOA clock error terms
% Define clock error statistics (in seconds)
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

%% Solve the Optimization Problem
% initialGuess = mean(selectedSatPositions, 1); % Use centroid as an initial guess
initialGuess = [lla2ecef([39.284593, 33.421097, 887]), 0];

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