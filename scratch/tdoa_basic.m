% Constants
c = 3e8; % Speed of light (m/s)
Ti = 20; % Time interval between localization signals (s)
num_instants = 10; % Number of localization signals
leo_velocity = [7500, 0, 0]; % LEO satellite velocity (m/s) [vx, vy, vz]
leo_initial_position = [0, 0, 500e3]; % Initial position of LEO (x, y, z in meters)
ue_position = [1e4, 2e4, 0]; % Ground UE position (x, y, z in meters)
noise_std = 0; % Standard deviation for TDOA noise (seconds)

% Initialize variables
leo_positions = zeros(num_instants, 3);
toa = zeros(1, num_instants);

% LEO Position Updates
leo_positions(1, :) = leo_initial_position;
for i = 2:num_instants
    leo_positions(i, :) = leo_positions(i-1, :) + leo_velocity * Ti;
end

% TOA Calculation
for i = 1:num_instants
    distance = norm(leo_positions(i, :) - ue_position); % Distance between LEO and UE
    toa(i) = distance / c; % Time of arrival
end

% TDOA Calculation with Noise
tdoa = diff(toa); % TDOA between successive signals
tdoa_noisy = tdoa + noise_std * randn(size(tdoa)); % Add Gaussian noise

% UE Localization using TDOA
% Set up least-squares localization
% TDOA equation: ||UE - LEO_2|| - ||UE - LEO_1|| = c * TDOA
A = zeros(num_instants - 1, 3);
b = zeros(num_instants - 1, 1);

for i = 1:(num_instants - 1)
    L1 = leo_positions(i, :);
    L2 = leo_positions(i+1, :);
    dL = L2 - L1;
    A(i, :) = dL / norm(L2 - ue_position); % Approximation for small TDOA
    b(i) = c * tdoa_noisy(i) - (norm(L2 - ue_position) - norm(L1 - ue_position));
end

% Solve for UE position using least squares
ue_estimated = pinv(A) * b;

% Display Results
disp('True UE Position:');
disp(ue_position);
disp('Estimated UE Position:');
disp(ue_estimated);

% Visualization
figure;
plot3(leo_positions(:, 1), leo_positions(:, 2), leo_positions(:, 3), '-o', 'LineWidth', 2);
hold on;
plot3(ue_position(1), ue_position(2), ue_position(3), 'r*', 'MarkerSize', 10);
plot3(ue_estimated(1), ue_estimated(2), ue_estimated(3), 'gs', 'MarkerSize', 10);
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('LEO Satellite Trajectory and UE Localization');
legend('LEO Trajectory', 'True UE Position', 'Estimated UE Position');


% Adaptive Hyperbola Plotting for All Instants
figure;
hold on;
grid on;
axis equal;

% Define a grid for visualization
x_range = linspace(-1e5, 1e5, 500); % Adjust range as needed
y_range = linspace(-1e5, 1e5, 500);
[X, Y] = meshgrid(x_range, y_range);

% Iterate over each pair of TDOA instants to plot hyperbolas
for i = 1:(num_instants - 1)
    % LEO positions for the pair
    L1 = leo_positions(i, 1:2);   % Project to 2D (x, y)
    L2 = leo_positions(i+1, 1:2);
    
    % Distances for grid points
    R1 = sqrt((X - L1(1)).^2 + (Y - L1(2)).^2);
    R2 = sqrt((X - L2(1)).^2 + (Y - L2(2)).^2);
    
    % Hyperbolic equation
    hyperbola = abs(R1 - R2) - c * abs(tdoa_noisy(i));
    
    % Plot contour where hyperbolic equation is approximately zero
    contour(X, Y, hyperbola, [0 0], 'LineWidth', 2, ...
        'DisplayName', sprintf('Hyperbola %d', i));
end

% Plot LEO positions
plot(leo_positions(:, 1), leo_positions(:, 2), 'ko', 'MarkerSize', 8, ...
    'DisplayName', 'LEO Positions');

% Plot true and estimated UE positions
plot(ue_position(1), ue_position(2), 'r*', 'MarkerSize', 10, ...
    'DisplayName', 'True UE Position');
plot(ue_estimated(1), ue_estimated(2), 'gs', 'MarkerSize', 10, ...
    'DisplayName', 'Estimated UE Position');

% Labels and legend
xlabel('X (m)');
ylabel('Y (m)');
title('Hyperbolas for Adaptive TDOA Localization');
legend('show');