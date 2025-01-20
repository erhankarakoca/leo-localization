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

carrierFreq = 10e9; % 10 GHz

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

%% GDOP based satellite selection
% Iterate over all subsets of satellite positions
numSats = size(accessedSatPositions, 1);
allCombinations = nchoosek(1:numSats, 4); % For 4 satellites at a time
gdop = zeros(size(allCombinations, 1),1);
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


%% Non Linear Least Square Solution

[estUEPosTDOAError, localizationErrorNLS] = NLSoptimizer( ...
    selectedSatPositions, pairs, radiiDifferenceswithError, initialGuess, lowerBound, upperBound, ueStationECEF);

disp('NonLin LS Localization Error (meters):');
disp(localizationErrorNLS);

[estUEPosTDOAError2, localizationErrorNLS] = NLSoptimizer( ...
    selectedSatPositions, pairs, radiiDifferenceswithError, estUEPosTDOAError, lowerBound, upperBound, ueStationECEF);

disp('NonLin LS Localization Error (meters):');
disp(localizationErrorNLS);

%% FDOA 
stdFOAError = 0; % Standard deviation of FOA noise in Hz

% Compute FOAs
FOAs = zeros(size(selectedSatPositions, 1), 1);
for i = 1:size(selectedSatPositions, 1)
    r_i = selectedSatPositions(i, :) - estUEPosTDOAError; 
    v_i = accessedSatVelocities(:,i); 
    FOAs(i) = (carrierFreq / c) * dot(v_i, r_i) / norm(r_i); % Doppler shift (FOA)
end

% Add noise to FOA
FOAsNoisy = FOAs + normrnd(0, stdFOAError, size(FOAs));

% Initialize FDOA array
FDOAs = zeros(size(pairs, 1), 1);

% Compute FDOA for each pair of satellites
for k = 1:size(pairs, 1)
    sat1 = pairs(k, 1);
    sat2 = pairs(k, 2);
    FDOAs(k) = FOAsNoisy(sat1) - FOAsNoisy(sat2); % FDOA as difference of FOAs
end

% Add noise to FDOA
stdFDOAError = 0; % Standard deviation of FDOA noise in Hz
FDOAsNoisy = FDOAs + normrnd(0, stdFDOAError, size(FDOAs));


% Combine TDOA and FDOA in the optimization
objective = @(p) objectiveTDOAFDOA(p, selectedSatPositions, accessedSatVelocities.', pairs, TDOAswithError, FDOAsNoisy, carrierFreq, c);

% Optimize
options = optimoptions('lsqnonlin');
[estimatedPosition, residual] = lsqnonlin(objective, initialGuess, lowerBound, upperBound, options);

% Display Results
disp('Estimated Position:');
disp(estimatedPosition);

localizationErrorTDOAFDOA = norm(estimatedPosition - ueStationECEF);
    disp('TDOA-FDOA Localization Error (meters):');
    disp(localizationErrorTDOAFDOA);

function error = objectiveTDOAFDOA(p, selectedSatPositions, accessedSatVelocities, pairs, TDOAs, FDOAs, f_c, c)
    error = zeros(size(pairs, 1) * 2, 1); 
    
    for k = 1:size(pairs, 1)
        sat1 = pairs(k, 1);
        sat2 = pairs(k, 2);
    
        r1 = selectedSatPositions(sat1, :) - p(1:3);
        r2 = selectedSatPositions(sat2, :) - p(1:3);
        v1 = accessedSatVelocities(sat1, :);
        v2 = accessedSatVelocities(sat2, :);
    
        d1 = norm(r1);
        d2 = norm(r2);
    
        % Compute TDOA residual
        tdoa_pred = (d1 - d2) / c;
        error(k) = TDOAs(k) - tdoa_pred;
    
        % Compute FDOA residual
        fdoa_pred = (f_c / c) * (dot(v1, r1) / d1 - dot(v2, r2) / d2);
        error(size(pairs, 1) + k) = FDOAs(k) - fdoa_pred;
    end
end



