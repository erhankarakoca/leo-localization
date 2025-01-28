% 3D TDOA Localization Validation with Monte Carlo Simulations

clear; clc; close all;

%% System Parameters
anchors = [0 0 0;   % Anchor positions (4 anchors in tetrahedron configuration)
           1 0 0;
           0 1 0;
           0 0 1];  
source_true = [0.3, 0.4, 0.5];  % True source position
c = 3e8;            % Speed of light (m/s)
SNR_dB = 20;        % Signal-to-Noise Ratio (dB)
B = 1e6;            % Bandwidth (Hz)
numMonteCarlo = 10000; % Number of Monte Carlo simulations

%% Derived Parameters
SNR_linear = 10^(SNR_dB/10);       % Convert SNR to linear scale
var_TOA = 1/(8*pi^2*B^2*SNR_linear); % TOA measurement variance
std_TOA = sqrt(var_TOA);           % TOA standard deviation

%% Preallocate Error Storage
errors_LLS = zeros(numMonteCarlo, 1);
errors_NLS = zeros(numMonteCarlo, 1);

%% Monte Carlo Simulation
for mc = 1:numMonteCarlo
    % Generate true TOA measurements
    true_dist = sqrt(sum((anchors - source_true).^2, 2));
    true_TOA = true_dist/c;
    
    % Add Gaussian noise to TOA measurements
    noisy_TOA = true_TOA + std_TOA*randn(size(true_TOA));
    
    % Calculate TDOA measurements (relative to first anchor)
    tdoa_meas = noisy_TOA(2:end) - noisy_TOA(1);
    
    %% Localization Methods
    % Method 1: Linear Least Squares (LLS)
    [est_LLS] = LLS_localization(anchors, tdoa_meas, c);
    
    % Method 2: Nonlinear Least Squares (NLS)
    [est_NLS] = NLS_localization(anchors, tdoa_meas, c, est_LLS);
    
    % Calculate errors
    errors_LLS(mc) = norm(est_LLS - source_true);
    errors_NLS(mc) = norm(est_NLS - source_true);
end

%% Performance Analysis
% Calculate RMSE
rmse_LLS = sqrt(mean(errors_LLS.^2));
rmse_NLS = sqrt(mean(errors_NLS.^2));

% Display results
fprintf('Localization RMSE Results:\n');
fprintf('LLS Method: %.4f meters\n', rmse_LLS);
fprintf('NLS Method: %.4f meters\n', rmse_NLS);

% Plot error distributions
figure;
histogram(errors_LLS, 'Normalization', 'pdf', 'BinWidth', 0.01);
hold on;
histogram(errors_NLS, 'Normalization', 'pdf', 'BinWidth', 0.01);
xlabel('Localization Error (m)');
ylabel('Probability Density');
title('Error Distribution Comparison');
legend('Linear LS', 'Nonlinear LS');
grid on;

%% Localization Functions
function [est_pos] = LLS_localization(anchors, tdoa_meas, c)
    % Linear Least Squares localization
    num_anchors = size(anchors, 1);
    A = [];
    b = [];
    ref_anchor = anchors(1,:);
    
    for i = 2:num_anchors
        x_i = anchors(i,1); y_i = anchors(i,2); z_i = anchors(i,3);
        x1 = ref_anchor(1); y1 = ref_anchor(2); z1 = ref_anchor(3);
        
        % Calculate matrix components
        K = x_i^2 + y_i^2 + z_i^2 - (x1^2 + y1^2 + z1^2);
        tau = tdoa_meas(i-1);
        
        A = [A; 2*(x_i - x1), 2*(y_i - y1), 2*(z_i - z1)];
        b = [b; (c*tau)^2 - K + 2*c*tau*norm(ref_anchor)];
    end
    
    % Solve linear system
    est_pos = (A\b)';
end

function [est_pos] = NLS_localization(anchors, tdoa_meas, c, init_guess)
    % Nonlinear Least Squares localization
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                          'Algorithm', 'levenberg-marquardt');
    
    cost_func = @(pos) tdoa_residuals(pos, anchors, tdoa_meas, c);
    est_pos = lsqnonlin(cost_func, init_guess, [], [], options);
end

function residuals = tdoa_residuals(pos, anchors, tdoa_meas, c)
    % Calculate TDOA residuals
    distances = sqrt(sum((anchors - pos).^2, 2));
    ref_dist = distances(1);
    tdoa_pred = (distances(2:end) - ref_dist)/c;
    residuals = tdoa_pred - tdoa_meas;
end