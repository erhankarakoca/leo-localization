%% Monte Carlo Simulation: FGO vs. NLS for TDoA Localization
clc; clear; close all;

%% 1. Simulation Parameters
c = 3e8;                   % Speed of light (m/s)
sigma_clock = 1e-9;        % Oscillator noise std (1 ns)
sigma_TDoA = 5e-9;         % TDoA measurement noise std (5 ns)
num_sats = 4;              % Number of LEO satellites
numTrials = 100;           % Number of Monte Carlo trials
rng(1);                   % Seed for reproducibility

% True UE position (fixed across trials)
true_P = [1000; 2000; 3000]; % [x, y, z] (meters)
true_delta = 1;              % True clock offset

% Adversarial satellite geometry (clustered)
S = 1e7 * [randn(1,3); randn(1,3) + 0.1; randn(1,3) + 0.2; randn(1,3) + 0.3]; 

%% 2. Preallocate Results Storage
errors_NLS = zeros(numTrials, 1);  % Position errors for NLS
errors_FGO = zeros(numTrials, 1);  % Position errors for FGO
delta_errors_NLS = zeros(numTrials, 1); % Clock offset errors for NLS
delta_errors_FGO = zeros(numTrials, 1); % Clock offset errors for FGO

%% 3. Monte Carlo Loop
for trial = 1:numTrials
    % Simulate TDoA measurements with new noise
    tau_true = zeros(num_sats, 1);
    true_delta = sigma_clock * randn(); % Vary delta across trials

    for i = 1:num_sats
        tau_true(i) = norm(true_P - S(i,:)) / c + true_delta;
    end
    TDoA_measurements = (tau_true(2:end) - tau_true(1)) ...
        + sqrt(2*sigma_clock) * randn(num_sats-1,1) ...  % Oscillator noise
        + sigma_TDoA * randn(num_sats-1,1);      % TDoA noise
    
    % Poor initial guess (far from true position)
    initial_guess = [1e4; 1e4; 1e4; 0];  % Bad initialization
    
    % Solve with NLS (no oscillator prior)
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off');
    x_opt_NLS = lsqnonlin(@(x) computeResiduals_NLS(x, S, TDoA_measurements,sigma_TDoA, c), ...
        initial_guess, [], [], options);
    errors_NLS(trial) = norm(x_opt_NLS(1:3) - true_P);
    delta_errors_NLS(trial) = abs(x_opt_NLS(4) - true_delta);
    
    % Solve with FGO (with oscillator prior)
    x_opt_FGO = lsqnonlin(@(x) computeResiduals_FGO(x, S, TDoA_measurements,sigma_TDoA, sigma_clock, c), ...
        initial_guess, [], [], options);
    errors_FGO(trial) = norm(x_opt_FGO(1:3) - true_P);
    delta_errors_FGO(trial) = abs(x_opt_FGO(4) - true_delta);
    
    % Display progress
    fprintf('Trial %d/%d completed.\n', trial, numTrials);
end

%% 4. Compute RMSE and Statistics
rmse_NLS = sqrt(mean(errors_NLS.^2));
rmse_FGO = sqrt(mean(errors_FGO.^2));
mean_delta_NLS = mean(delta_errors_NLS);
mean_delta_FGO = mean(delta_errors_FGO);

%%
figure;
boxplot([errors_NLS, errors_FGO], 'Labels', {'NLS', 'FGO'});
ylabel('Position Error (meters)');
title('FGO vs. NLS (Adversarial Conditions)');

fprintf('=== Results ===\n');
fprintf('NLS RMSE: %.2f meters\n', rmse_NLS);
fprintf('FGO RMSE: %.2f meters\n', rmse_FGO);

%% 5. Plot Results
% Position RMSE comparison
figure;
bar([rmse_NLS, rmse_FGO]);
set(gca, 'XTickLabel', {'Standard NLS', 'FGO'});
ylabel('Position RMSE (meters)');
title('RMSE Comparison (Monte Carlo Simulation)');

% Clock offset error comparison
figure;
bar([mean_delta_NLS, mean_delta_FGO]);
set(gca, 'XTickLabel', {'Standard NLS', 'FGO'});
ylabel('Mean Clock Offset Error (seconds)');
title('Clock Offset Error Comparison');

% Histogram of position errors
figure;
histogram(errors_NLS, 'BinWidth', 10, 'FaceColor', 'b', 'EdgeColor', 'none'); hold on;
histogram(errors_FGO, 'BinWidth', 10, 'FaceColor', 'r', 'EdgeColor', 'none');
xlabel('Position Error (meters)');
ylabel('Frequency');
legend('NLS', 'FGO');
title('Distribution of Position Errors');

%% 6. Residual Functions
function residuals = computeResiduals_NLS(x, S, TDoA_meas,sigma_TDoA, c)
    P = x(1:3);
    delta = x(4);
    residuals = zeros(length(TDoA_meas), 1);
    for i = 1:length(TDoA_meas)
        sat_j = i + 1;
        tau_j = norm(P - S(sat_j,:)) / c + delta;
        tau_1 = norm(P - S(1,:)) / c + delta;
        residuals(i) = (tau_j - tau_1 - TDoA_meas(i))/sigma_TDoA;
    end
end

function residuals = computeResiduals_FGO(x, S, TDoA_meas,sigma_TDoA, sigma_clock, c)
    P = x(1:3);
    delta = x(4);
    num_meas = length(TDoA_meas);
    residuals = zeros(num_meas + 1, 1);
    for i = 1:num_meas
        sat_j = i + 1;
        tau_j = norm(P - S(sat_j,:)) / c + delta;
        tau_1 = norm(P - S(1,:)) / c + delta;
        residuals(i) = (tau_j - tau_1 - TDoA_meas(i))/sigma_TDoA;
    end
    residuals(end) = delta / sigma_clock;
end