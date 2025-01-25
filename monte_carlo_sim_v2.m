%% Clear Workspace
clear; close all; clc;

%% Satellite Scenario Setup (UTC Time Zone)
startTime = datetime(2020, 5, 4, 18, 45, 50, 'TimeZone', 'UTC');
stopTime = datetime(2020, 5, 4, 19, 2, 20, 'TimeZone', 'UTC');
sampleTime = 10;
satscene = satelliteScenario(startTime, stopTime, sampleTime);
c = physconst("LightSpeed");

% Add satellites from TLE
tleFile = "leoSatelliteConstellation.tle";
constellation = satellite(satscene, tleFile);

% UE Ground Station (Ankara)
ueStationLLA = [40.786648, 29.449502, 182];
ueStationECEF = lla2ecef(ueStationLLA);
gsUE = groundStation(satscene, "Latitude", ueStationLLA(1), ...
                    "Longitude", ueStationLLA(2), "Altitude", ueStationLLA(3));

%% Monte Carlo Parameters
numTrials = 100;            % Reduced for testing
sigma_clock = 1e-9;         % Clock jitter std (1 ns)
sigma_TDoA = 1e-9;          % TDoA measurement noise std (1 ns)
sigma_total = sqrt(2*sigma_clock^2 + sigma_TDoA^2);

% Preallocate results
errorsNLS = zeros(numTrials, 1);
errorsLM = zeros(numTrials, 1);      % Levenberg-Marquardt
errorsPSO = zeros(numTrials, 1);     % Particle Swarm
errorsHuber = zeros(numTrials, 1);

%% Helper Functions
function gdop = calculateGDOP(satPositions, truePos)
    H = [(satPositions - truePos)./vecnorm(satPositions - truePos, 2, 2), -ones(size(satPositions,1),1)];
    gdop = sqrt(trace(inv(H'*H)));
end

function initialGuess = tdoaClosedForm(positions, tdoaMeas, c)
    % Linearized least squares solution
    refPos = positions(1,:);
    numPairs = size(tdoaMeas,1);
    
    A = zeros(numPairs,3);
    b = zeros(numPairs,1);
    
    for k = 1:numPairs
        j = k+1; % Pair with reference
        rij = tdoaMeas(k)*c;
        A(k,:) = 2*(positions(j,:) - refPos);
        b(k) = sum(positions(j,:).^2) - sum(refPos.^2) - rij^2 - 2*rij*norm(refPos);
    end
    
    initialGuess = (A \ b)';
end

%% Enhanced Monte Carlo Loop
parfor mc = 1:numTrials
    rng(mc); % Reproducible randomness
    
    try
        % Generate random UTC time
        totalSamples = seconds(stopTime - startTime)/sampleTime;
        randomSampleIndex = randi([0, totalSamples-1]);
        randomTime = startTime + seconds(randomSampleIndex*sampleTime);
        
        % Find accessible satellites
        ac = access(constellation, gsUE);
        accessTable = accessIntervals(ac);
        validSats = {};
        for i = 1:height(accessTable)
            if randomTime >= accessTable{i,4} && randomTime <= accessTable{i,5}
                validSats = [validSats; accessTable{i,1}];
            end
        end
        if numel(validSats) < 4, continue; end
        
        % Propagate satellite positions
        tleStruct = tleread('leoSatelliteConstellation.tle');
        satIndices = find(ismember({tleStruct.Name}, validSats));
        [satPositions, ~] = propagateOrbit(randomTime, tleStruct(satIndices), ...
                                          "OutputCoordinateSystem", "fixed-frame");
        satPositions = squeeze(satPositions)';
        
        % GDOP-based satellite selection
        numSats = size(satPositions,1);
        combs = nchoosek(1:numSats,4);
        gdops = arrayfun(@(i) calculateGDOP(satPositions(combs(i,:),:), ueStationECEF), 1:size(combs,1));
        [~,bestIdx] = min(gdops);
        selectedSats = combs(bestIdx,:);
        selectedPos = satPositions(selectedSats,:);
        
        % Generate measurements
        [~,~,distances] = aer(gsUE, constellation(satIndices(selectedSats)), randomTime);
        TOAs = squeeze(distances)/c;
        pairs = nchoosek(1:4,2);
        trueTDOA = TOAs(pairs(:,1)) - TOAs(pairs(:,2));
        noisyTDOA = trueTDOA + normrnd(0, sigma_total, size(trueTDOA));
        
        % Closed-form initial guess
        initialGuess = tdoaClosedForm(selectedPos, noisyTDOA, c);
        
        % Optimization setup
        earthRadius = 6.371e6;
        lb = -earthRadius*ones(1,3);
        ub = earthRadius*ones(1,3);
        
        % 1. Standard NLS
        objNLS = @(p) arrayfun(@(k) (norm(p - selectedPos(pairs(k,1),:)) - ...
                                    norm(p - selectedPos(pairs(k,2),:)) - ...
                                    noisyTDOA(k)*c)/sigma_total, 1:size(pairs,1));
        options = optimoptions('lsqnonlin', 'Display','off', 'MaxIterations',1000);
        estNLS = lsqnonlin(objNLS, initialGuess, lb, ub, options);
        
        % 2. Levenberg-Marquardt
        optionsLM = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', ...
                               'Display','off', 'MaxIterations',1000);
        estLM = lsqnonlin(objNLS, initialGuess, lb, ub, optionsLM);
        
        % 3. Particle Swarm Optimization
        objPSO = @(p) sum(objNLS(p).^2);
        optionsPSO = optimoptions('particleswarm', 'Display','off', ...
                                'SwarmSize', 100, 'MaxIterations', 200);
        estPSO = particleswarm(objPSO, 3, lb, ub, optionsPSO);
        
        % 4. Huber Loss with fminunc
        k = 1.345; % Huber parameter
        huberObj = @(p) sum((abs(objNLS(p)) <= k).*0.5.*objNLS(p).^2 + ...
                           (abs(objNLS(p)) > k).*(k*abs(objNLS(p)) - 0.5*k^2));
        optionsHuber = optimoptions('fminunc', 'Algorithm','quasi-newton', ...
                                  'Display','off', 'MaxIterations', 1000);
        estHuber = fminunc(huberObj, initialGuess, optionsHuber);
        
        % Store results
        errorsNLS(mc) = norm(estNLS - ueStationECEF);
        errorsLM(mc) = norm(estLM - ueStationECEF);
        errorsPSO(mc) = norm(estPSO - ueStationECEF);
        errorsHuber(mc) = norm(estHuber - ueStationECEF);
        
    catch ME
        fprintf('Error in trial %d: %s\n', mc, ME.message);
        continue;
    end
end

%% Post-processing
validTrials = errorsNLS ~= 0;
errorsNLS = errorsNLS(validTrials);
errorsLM = errorsLM(validTrials);
errorsPSO = errorsPSO(validTrials);
errorsHuber = errorsHuber(validTrials);

% Calculate RMSE
rmse = @(x) sqrt(mean(x.^2));
fprintf('=== Localization Performance ===\n');
fprintf('Standard NLS RMSE: %.2f m\n', rmse(errorsNLS));
fprintf('Levenberg-Marquardt RMSE: %.2f m\n', rmse(errorsLM));
fprintf('Particle Swarm RMSE: %.2f m\n', rmse(errorsPSO));
fprintf('Huber Robust RMSE: %.2f m\n', rmse(errorsHuber));

% Visualization
figure;
boxplot([errorsNLS, errorsLM, errorsPSO, errorsHuber], ...
       'Labels', {'NLS', 'LM', 'PSO', 'Huber'});
ylabel('Localization Error (m)');
title('Optimization Algorithm Comparison');
grid on;