% Clock error statistics
sigma_clock = 1e-6; % Example: 1 microsecond
c = 3e8; % Speed of light (m/s)
clock_error_std = c * sigma_clock; % Standard deviation in distance

% Example geometry matrix (replace with actual satellite geometry matrix)
G = [3, 0.5, 0.2; 0.5, 2, 0.3; 0.2, 0.3, 1.5]; % Example covariance-like matrix

% Covariance matrix due to clock error
C = clock_error_std^2 * inv(G);

% Eigen decomposition of covariance matrix
[V, D] = eig(C); % V: Eigenvectors, D: Diagonal matrix with eigenvalues

% Semi-principal axes of the ellipsoid
a = sqrt(D(1,1)); % Semi-axis along X
b = sqrt(D(2,2)); % Semi-axis along Y
c = sqrt(D(3,3)); % Semi-axis along Z

% Generate ellipsoid in 3D space
[x, y, z] = ellipsoid(0, 0, 0, a, b, c, 50); % 50 points for smooth surface

% Rotate ellipsoid using eigenvectors
ellipsoid_points = [x(:), y(:), z(:)] * V'; % Rotate using eigenvectors
x_rot = reshape(ellipsoid_points(:,1), size(x)) + estimated_position(1);
y_rot = reshape(ellipsoid_points(:,2), size(y)) + estimated_position(2);
z_rot = reshape(ellipsoid_points(:,3), size(z)) + estimated_position(3);

% Plotting
figure;
hold on;
surf(x_rot, y_rot, z_rot, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Ellipsoid
plot3(estimated_position(1), estimated_position(2), estimated_position(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % UE position
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('3D Uncertainty Ellipsoid for Clock Error');
legend('Uncertainty Ellipsoid', 'Estimated Position');
axis equal;
hold off;
