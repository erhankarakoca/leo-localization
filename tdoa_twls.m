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
altitudeMax = 1e5; % Maximum altitude
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

%% Non Linear Least Square Solution
ObjectiveEstimationErrorTDOA = @(p) arrayfun(@(k) ...
    abs(sqrt((p(1) - selectedSatPositions(pairs(k, 1), 1))^2 + ...
             (p(2) - selectedSatPositions(pairs(k, 1), 2))^2 + ...
             (p(3) - selectedSatPositions(pairs(k, 1), 3))^2) - ...
        sqrt((p(1) - selectedSatPositions(pairs(k, 2), 1))^2 + ...
             (p(2) - selectedSatPositions(pairs(k, 2), 2))^2 + ...
             (p(3) - selectedSatPositions(pairs(k, 2), 3))^2) - ...
        radiiDifferenceswithError(k)), 1:size(pairs, 1));

options = optimoptions('lsqnonlin');

estUEPosTDOAError = lsqnonlin(ObjectiveEstimationErrorTDOA, ...
                                    initialGuess(1:3), ...
                                    lowerBound, ...
                                    upperBound, ...
                                    options);
localizationErrorNLS = norm(estUEPosTDOAError(1:3) - actualUEPosition);

%% TWLS
[p_est, residual] = TWLSoptimizer(selectedSatPositions, TDOAswithError, pairs, covarianceMatrices, ...
                                    initialGuess, 10, 1e-9, c);

% Display results
disp('NonLin LS Localization Error (meters):');
disp(localizationErrorNLS);

localizationErrorTWLS = norm(p_est - actualUEPosition);
disp('TWLS Localization Error (meters):');
disp(localizationErrorTWLS);


