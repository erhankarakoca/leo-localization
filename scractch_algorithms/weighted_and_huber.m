%% TDoA Localization: Comparing Standard NLS, Weighted NLS, and Robust NLS
clc; clear; close all;

%% 1. Simulation Parameters
c = physconst("lightspeed");                    % Speed of light (m/s)
sigma_TDoA = 2e-9;%1e-9;          % TDoA noise std (5 ns)
sigma_clock = 1e-6;         % Clock jitter std (1 ns)
num_sats = 5;               % Number of satellites
numTrials = 1000;            % Monte Carlo trials
rng(42);                    % Seed for reproducibility  

% True UE position (fixed across trials)
true_P = [1000; 2000; 3000];

% Generate satellites in a tetrahedral configuration (good geometry)
% S = 1e7 * [1, 1, 1; -1, -1, 1; -1, 1, -1; 1, -1, -1]; 
S = 1e7 * [1, 1, 1; -1, -1, 1; -1, 1, -1; 1, -1, -1; 0, 0, 1; 0, 1, 0]; 

% Total noise variance (TDoA + 2*clock jitter)
sigma_total = sqrt(sigma_TDoA^2 + 2 * sigma_clock^2);

%% 2. Preallocate Results
errors_standard = zeros(numTrials, 1);
errors_weighted = zeros(numTrials, 1);
errors_robust = zeros(numTrials, 1);

%% 3. Monte Carlo Loop
parfor trial = 1:numTrials
    % Simulate TDoA measurements with clock jitter
    TDoA_meas = simulateTDoA(true_P, S, sigma_TDoA, sigma_clock, c);
    
    % Initial guess (centroid of satellites)
    initial_guess = [mean(S(:,1)), mean(S(:,2)), mean(S(:,3))]';
    % initial_guess = [0,0,0]';
    % --- Method 1: Standard NLS (no noise modeling) ---
    costFunc_standard = @(x) computeResiduals_standard(x, S, TDoA_meas, c);
    x_opt_standard = lsqnonlin(costFunc_standard, initial_guess, [], [], optimoptions('lsqnonlin', 'Display', 'off'));
    errors_standard(trial) = norm(x_opt_standard - true_P);
    
    % --- Method 2: Weighted NLS (accounts for total noise) ---
    costFunc_weighted = @(x) computeResiduals_weighted(x, S, TDoA_meas, sigma_total, c);
    x_opt_weighted = lsqnonlin(costFunc_weighted, initial_guess, [], [], optimoptions('lsqnonlin', 'Display', 'off'));
    errors_weighted(trial) = norm(x_opt_weighted - true_P);
    
    % --- Method 3: Robust NLS (Huber loss) ---
    costFunc_robust = @(x) computeResiduals_robust(x, S, TDoA_meas, sigma_total, c);
    x_opt_robust = lsqnonlin(costFunc_robust, initial_guess, [], [], optimoptions('lsqnonlin', 'Display', 'off'));
    errors_robust(trial) = norm(x_opt_robust - true_P);
    
    fprintf('Trial %d/%d\n', trial, numTrials);
end

%% 4. Compute RMSE
rmse_standard = sqrt(mean(errors_standard.^2));
rmse_weighted = sqrt(mean(errors_weighted.^2));
rmse_robust = sqrt(mean(errors_robust.^2));

%% 5. Plot Results
figure;
boxplot([errors_standard, errors_weighted, errors_robust], ...
    'Labels', {'Standard NLS', 'Weighted NLS', 'Robust NLS'});
ylabel('Position Error (meters)');
title('TDoA Localization Performance Comparison');

fprintf('=== RMSE Results ===\n');
fprintf('Standard NLS: %.2f meters\n', rmse_standard);
fprintf('Weighted NLS: %.2f meters\n', rmse_weighted);
fprintf('Robust NLS: %.2f meters\n', rmse_robust);

%% 6. Helper Functions
function TDoA_meas = simulateTDoA(true_P, S, sigma_TDoA, sigma_clock, c)
    num_sats = size(S, 1);
    tau = zeros(num_sats, 1);
    for i = 1:num_sats
        delta_i = sigma_clock * randn(); % Clock jitter per satellite
        tau(i) = norm(true_P - S(i,:)) / c + delta_i;
    end
    TDoA_meas = (tau(2:end) - tau(1)) + sigma_TDoA * randn(num_sats-1, 1);
end

function residuals = computeResiduals_standard(x, S, TDoA_meas, c)
    P = x(1:3);
    residuals = zeros(length(TDoA_meas), 1);
    for i = 1:length(TDoA_meas)
        sat_j = i + 1;
        tau_j = norm(P - S(sat_j,:)) / c;
        tau_1 = norm(P - S(1,:)) / c;
        residuals(i) = (tau_j - tau_1 - TDoA_meas(i));
    end
end

function residuals = computeResiduals_weighted(x, S, TDoA_meas, sigma_total, c)
    P = x(1:3);
    residuals = zeros(length(TDoA_meas), 1);
    for i = 1:length(TDoA_meas)
        sat_j = i + 1;
        tau_j = norm(P - S(sat_j,:)) / c;
        tau_1 = norm(P - S(1,:)) / c;
        residuals(i) = (tau_j - tau_1 - TDoA_meas(i)) / sigma_total;
    end
end

function residuals = computeResiduals_robust(x, S, TDoA_meas, sigma_total, c)
    % Huber loss parameters (k = 1.345 for 95% Gaussian efficiency)
    k = 1.345;
    P = x(1:3);
    residuals = zeros(length(TDoA_meas), 1);
    for i = 1:length(TDoA_meas)
        sat_j = i + 1;
        tau_j = norm(P - S(sat_j,:)) / c;
        tau_1 = norm(P - S(1,:)) / c;
        r = (tau_j - tau_1 - TDoA_meas(i)) / sigma_total;
        % Apply Huber weighting
        if abs(r) <= k
            residuals(i) = r;
        else
            residuals(i) = sign(r) * sqrt(2 * k * abs(r) - k^2);
        end
    end
end