clear all; close all; clc;

%% Satellite Scenario Setup
startTime = datetime(2020, 05, 04, 18,45,50);
stopTime = datetime(2020, 05, 04, 19,02,20);
sampleTime = 10;
satscene = satelliteScenario(startTime, stopTime, sampleTime);
c = physconst("LightSpeed");

% Add satellites from TLE
tleFile = "leoSatelliteConstellation.tle";
constellation = satellite(satscene, tleFile);

% UE Ground Station
ueStationLLA = [40.786648, 29.449502, 182];
ueStationECEF = lla2ecef(ueStationLLA);
gsUE = groundStation(satscene, "Latitude", ueStationLLA(1), ...
                     "Longitude", ueStationLLA(2), "Altitude", ueStationLLA(3));

%% Monte Carlo Parameters
numTrials = 1000;            % Number of Monte Carlo trials
sigma_clock = 0;         % Clock jitter std (1 µs)
sigma_TDoA = 0;          % TDoA measurement noise std (5 µs)
sigma_total = sqrt(2 * sigma_clock^2 + sigma_TDoA^2); % Combined noise

% Preallocate results
errorsNLS = zeros(numTrials, 1);
errorsHuber = zeros(numTrials, 1);

%% Monte Carlo Loop
parfor mc = 1:numTrials
    rng(mc); % Set unique seed for each trial (for reproducibility)
    
    % Randomly select time and find accessed satellites
    totalSamples = seconds(stopTime - startTime) / sampleTime;
    randomSampleIndex = randi([0, totalSamples - 1]);
    randomDateTime = startTime + seconds(randomSampleIndex * sampleTime);
    randomDateTime = datetime(randomDateTime, 'TimeZone', 'UTC');
    
    % Find accessed satellites
    ac = access(constellation, gsUE);
    accessIntervalsTable = accessIntervals(ac); % Correctly extract access intervals
    
    accessedSatellites = [];
    for i = 1:height(accessIntervalsTable)
        accessStartTime = accessIntervalsTable{i, 4};
        accessEndTime = accessIntervalsTable{i, 5};
        if randomDateTime >= accessStartTime && randomDateTime <= accessEndTime
            accessedSatellites = [accessedSatellites; accessIntervalsTable{i, 1}];
        end
    end
    
    % Skip trials with <4 satellites
    if numel(accessedSatellites) < 4
        continue;
    end
    
    % Propagate satellite positions
    tleStruct = tleread('leoSatelliteConstellation.tle');
    satelliteNamesInTLE = {tleStruct.Name}';
    indicesInTLE = find(matches(string(satelliteNamesInTLE), accessedSatellites));
    [accessedSatPositions, ~] = propagateOrbit(randomDateTime, tleStruct(indicesInTLE), ...
                                               "OutputCoordinateFrame", "fixed-frame");
    accessedSatPositions = squeeze(accessedSatPositions)';
    
    % Compute distances and TOAs
    [~, ~, accessedSatDistances] = aer(gsUE, constellation(indicesInTLE), randomDateTime);
    TOAs = squeeze(accessedSatDistances) / c;
    
    % Generate TDOA measurements with noise
    pairs = nchoosek(1:numel(TOAs), 2);
    TDOAs = TOAs(pairs(:,1)) - TOAs(pairs(:,2));
    noise = normrnd(0, sigma_total, size(TDOAs)); % Combined noise
    TDOAswithError = TDOAs + noise;
    
    % GDOP-based satellite selection (best 4 satellites)
    numSats = size(accessedSatPositions, 1);
    allCombinations = nchoosek(1:numSats, 4);
    gdop = zeros(size(allCombinations, 1), 1);
    for i = 1:size(allCombinations, 1)
        subset = allCombinations(i, :);
        subsetPositions = accessedSatPositions(subset, :);
        [gdop(i), ~] = calculateGDOP(subsetPositions, ueStationECEF);
    end
    [~, idx] = min(gdop);
    selectedSatIndices = allCombinations(idx, :);
    selectedSatPositions = accessedSatPositions(selectedSatIndices, :);
    
    % Recompute pairs for selected satellites
    selectedPairs = nchoosek(1:4, 2); % Pairs for 4 selected satellites
    selectedTDOAswithError = TDOAswithError(selectedPairs(:,1) + selectedPairs(:,2) - 1);
    
    % Define optimization functions
    initialGuess = lla2ecef([39.284593, 33.421097, 887]); % Initial guess
    earthRadius = 6.371e6;
    lb = -earthRadius * ones(1, 3);
    ub = earthRadius * ones(1, 3);
    options = optimoptions('lsqnonlin', 'Display', 'off');
    
    % Normal NLS
    ObjectiveNLS = @(p) arrayfun(@(k) ...
        (norm(p - selectedSatPositions(selectedPairs(k,1), :)) - ...
         norm(p - selectedSatPositions(selectedPairs(k,2), :)) - ...
         selectedTDOAswithError(k) * c) , 1:size(selectedPairs,1));
    estNLS = lsqnonlin(ObjectiveNLS, initialGuess, lb, ub, options);
    
    % Huber Loss
    k = 1.5; % Huber parameter
    huberLoss = @(r) (abs(r) <= k) .* r.^2 + (abs(r) > k) .* (2*k*abs(r) - k^2);
    ObjectiveHuber = @(p) arrayfun(@(k) ...
        sqrt(huberLoss((norm(p - selectedSatPositions(selectedPairs(k,1), :)) - ...
                        norm(p - selectedSatPositions(selectedPairs(k,2), :)) - ...
                        selectedTDOAswithError(k) * c) )), 1:size(selectedPairs,1));
    estHuber = lsqnonlin(ObjectiveHuber, initialGuess, lb, ub, options);
    
    % Compute errors
    errorsNLS(mc) = norm(estNLS - ueStationECEF);
    errorsHuber(mc) = norm(estHuber - ueStationECEF);
end

%% Results Analysis
% Remove skipped trials (errors = 0)
errorsNLS = errorsNLS(errorsNLS ~= 0);
errorsHuber = errorsHuber(errorsHuber ~= 0);

% Compute RMSE
rmseNLS = sqrt(mean(errorsNLS.^2));
rmseHuber = sqrt(mean(errorsHuber.^2));

% Plot error distributions
figure;
boxplot([errorsNLS, errorsHuber], 'Labels', {'Normal NLS', 'Huber Loss'});
ylabel('Localization Error (meters)');
title('Error Distributions (Monte Carlo)');

% Display RMSE
disp('=== Monte Carlo Results ===');
disp(['Normal NLS RMSE: ', num2str(rmseNLS), ' meters']);
disp(['Huber Loss RMSE: ', num2str(rmseHuber), ' meters']);