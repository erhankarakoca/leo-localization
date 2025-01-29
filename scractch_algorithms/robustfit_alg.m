%% TDoA Localization: Comparing Standard NLS, Weighted NLS, Robust NLS (Huber), and RobustFit
clc; clear; close all;

%% 1. Simulation Parameters
c = 3e8;                    % Speed of light (m/s)
sigma_TDoA = 5e-9;          % TDoA noise std (5 ns)
sigma_clock = 1e-9;         % Clock jitter std (1 ns)
num_sats = 5;               % Use 6 satellites (now 5 TDoA measurements)
numTrials = 1000;            % Monte Carlo trials
maxIter = 10;               % Max iterations for RobustFit
rng(42);                    % Seed for reproducibility

% True UE position (fixed across trials)
true_P = [1000; 2000; 3000];

% Generate satellites in a diverse 3D configuration
S = 1e7 * [1, 1, 1; -1, -1, 1; -1, 1, -1; 1, -1, -1; 0, 0, 1; 0, 1, 0]; 

% Total noise variance (TDoA + 2*clock jitter)
sigma_total = sqrt(sigma_TDoA^2 + 2 * sigma_clock^2);

%% 2. Preallocate Results
errors_standard = zeros(numTrials, 1);
errors_weighted = zeros(numTrials, 1);
errors_robust = zeros(numTrials, 1);
errors_robustfit = zeros(numTrials, 1);

%% 3. Monte Carlo Loop
parfor trial = 1:numTrials
    % Simulate TDoA measurements with clock jitter
    TDoA_meas = simulateTDoA(true_P, S, sigma_TDoA, sigma_clock, c);
    
    % Initial guess (centroid of satellites)
    initial_guess = [mean(S(:,1)), mean(S(:,2)), mean(S(:,3))]';
    
    % --- Method 1: Standard NLS (no noise modeling) ---
    costFunc_standard = @(x) computeResiduals_standard(x, S, TDoA_meas, c);
    x_opt_standard = lsqnonlin(costFunc_standard, initial_guess, [], [], optimoptions('lsqnonlin', 'Display', 'off'));
    errors_standard(trial) = norm(x_opt_standard - true_P);
    
    % --- Method 2: Weighted NLS (accounts for total noise) ---
    costFunc_weighted = @(x) computeResiduals_weighted(x, S, TDoA_meas, sigma_total, c);
    x_opt_weighted = lsqnonlin(costFunc_weighted, initial_guess, [], [], optimoptions('lsqnonlin', 'Display', 'off'));
    errors_weighted(trial) = norm(x_opt_weighted - true_P);
    
    % --- Method 3: Robust NLS (Huber loss via lsqnonlin) ---
    costFunc_robust = @(x) computeResiduals_robust(x, S, TDoA_meas, sigma_total, c);
    x_opt_robust = lsqnonlin(costFunc_robust, initial_guess, [], [], optimoptions('lsqnonlin', 'Display', 'off'));
    errors_robust(trial) = norm(x_opt_robust - true_P);
    
    % --- Method 4: Iterative RobustFit (with linearization) ---
    P_robustfit = iterativeRobustFit(TDoA_meas, S, initial_guess, sigma_total, c, maxIter);
    errors_robustfit(trial) = norm(P_robustfit - true_P);
    
    fprintf('Trial %d/%d\n', trial, numTrials);
end

%% 4. Compute RMSE
rmse_standard = sqrt(mean(errors_standard.^2));
rmse_weighted = sqrt(mean(errors_weighted.^2));
rmse_robust = sqrt(mean(errors_robust.^2));
rmse_robustfit = sqrt(mean(errors_robustfit.^2));

%% 5. Plot Results
figure;
boxplot([errors_standard, errors_weighted, errors_robust, errors_robustfit], ...
    'Labels', {'Standard NLS', 'Weighted NLS', 'Robust NLS (Huber)', 'RobustFit'});
ylabel('Position Error (meters)');
title('TDoA Localization Performance Comparison');

fprintf('=== RMSE Results ===\n');
fprintf('Standard NLS: %.2f meters\n', rmse_standard);
fprintf('Weighted NLS: %.2f meters\n', rmse_weighted);
fprintf('Robust NLS (Huber): %.2f meters\n', rmse_robust);
fprintf('RobustFit: %.2f meters\n', rmse_robustfit);

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
    k = 1.345; % Huber parameter for 95% Gaussian efficiency
    P = x(1:3);
    residuals = zeros(length(TDoA_meas), 1);
    for i = 1:length(TDoA_meas)
        sat_j = i + 1;
        tau_j = norm(P - S(sat_j,:)) / c;
        tau_1 = norm(P - S(1,:)) / c;
        r = (tau_j - tau_1 - TDoA_meas(i)) / sigma_total;
        residuals(i) = (abs(r) <= k) * r + (abs(r) > k) * (sign(r) * sqrt(2 * k * abs(r) - k^2));
    end
end

function P_opt = iterativeRobustFit(TDoA_meas, S, initial_guess, sigma_total, c, maxIter)
    P_opt = initial_guess;
    tol = 1e-3;
    for iter = 1:maxIter
        % Compute Jacobian and residuals
        [res, J] = computeResidualsAndJacobian(P_opt, S, TDoA_meas, sigma_total, c);
        
        % Add small regularization to avoid rank deficiency
        J_reg = [J; eye(3) * 1e-6]; % Augment Jacobian
        res_reg = [res; zeros(3, 1)]; % Augment residuals
        
        % Solve using robustfit (no intercept, bisquare weighting)
        [delta_P, ~] = robustfit(J_reg, res_reg, 'bisquare', [], 'off');
        
        % Update position
        P_opt = P_opt + delta_P;
        
        % Check convergence
        if norm(delta_P) < tol
            break;
        end
    end
end

function [res, J] = computeResidualsAndJacobian(P, S, TDoA_meas, sigma_total, c)
    num_meas = length(TDoA_meas);
    res = zeros(num_meas, 1);
    J = zeros(num_meas, 3); % Jacobian has 3 columns (x, y, z gradients)
    
    for i = 1:num_meas
        sat_j = i + 1;
        
        % Ensure satellite positions are column vectors
        S_j = S(sat_j, :)'; % Transpose to column vector (3x1)
        S_1 = S(1, :)';     % Transpose to column vector (3x1)
        
        % Compute distances (avoid division by zero with a small epsilon)
        eps = 1e-3;
        d_j = max(norm(P - S_j), eps);
        d_1 = max(norm(P - S_1), eps);
        
        % Predicted TDoA
        tau_pred = (d_j - d_1) / c;
        
        % Residual (scaled by total noise)
        res(i) = (TDoA_meas(i) - tau_pred) / sigma_total;
        
        % Compute gradients (3x1 column vectors)
        grad_j = (P - S_j) / d_j;
        grad_1 = (P - S_1) / d_1;
        
        % Jacobian row: transpose to 1x3 row vector
        J(i, :) = ((grad_j - grad_1) / (c * sigma_total))';
    end
end