% Enhanced 3D TDOA Localization with Clock Instability Analysis

clear; clc; close all;

%% System Parameters
anchors = [0 0 0;   % Anchor positions
           1 0 0;
           0 1 0;
           0 0 1]*1e6;  
source_true = [0.3, 0.4, 0.5]*1e3;  % True source position
c = 3e8;            % Speed of light (m/s)
SNR_dB = 10;        % Signal-to-Noise Ratio (dB)
B = 1e6;            % Bandwidth (Hz)
clock_jitter_std = 1e-9; % 1 ns standard deviation
numMonteCarlo = 1000; 

%% Derived Parameters
SNR_linear = 10^(SNR_dB/10);
var_TOA = 1/(8*pi^2*B^2*SNR_linear);
std_TOA = sqrt(var_TOA);
total_TOA_std = sqrt(std_TOA^2 + clock_jitter_std^2);

fprintf('Error Components:\n');
fprintf('Thermal noise std: %.2e s\n', std_TOA);
fprintf('Clock jitter std:  %.2e s\n', clock_jitter_std);
fprintf('Combined std:      %.2e s\n\n', total_TOA_std);

%% Preallocate Error Storage
errors_LLS = zeros(numMonteCarlo, 1);
errors_NLS = zeros(numMonteCarlo, 1);

%% Monte Carlo Simulation
for mc = 1:numMonteCarlo
    % True TOA measurements
    true_dist = sqrt(sum((anchors - source_true).^2, 2));
    true_TOA = true_dist/c;
    
    % Add combined noise (thermal + clock jitter)
    noisy_TOA = true_TOA + std_TOA*randn(size(true_TOA)) + ...
                clock_jitter_std*randn(size(true_TOA));
    
    % Calculate TDOA measurements
    tdoa_meas = noisy_TOA(2:end) - noisy_TOA(1);
    
    %% Localization Methods
    [est_LLS] = LLS_localization(anchors, tdoa_meas, c);
    [est_NLS] = NLS_localization(anchors, tdoa_meas, c, est_LLS);
    
    errors_LLS(mc) = norm(est_LLS - source_true);
    errors_NLS(mc) = norm(est_NLS - source_true);
end

%% Performance Analysis
rmse_LLS = sqrt(mean(errors_LLS.^2));
rmse_NLS = sqrt(mean(errors_NLS.^2));

fprintf('Localization RMSE Results:\n');
fprintf('LLS Method: %.4f meters\n', rmse_LLS);
fprintf('NLS Method: %.4f meters\n', rmse_NLS);

figure;
hold on;
histogram(errors_LLS, 'Normalization', 'pdf', 'BinWidth', 0.01);
histogram(errors_NLS, 'Normalization', 'pdf', 'BinWidth', 0.01);
xlabel('Localization Error (m)');
ylabel('Probability Density');
title('Error Distribution with Clock Instability');
legend('Linear LS', 'Nonlinear LS');
grid on;

%% Localization Functions (Same as previous)
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