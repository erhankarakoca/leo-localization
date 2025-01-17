function [p_est, residual] = TWLSoptimizer(selectedSatPositions, TDOAs, pairs, covarianceMatrices, initialGuess, maxIter, tol, c)
    % Implements Taylor Weighted Least Squares (TWLS) for TDOA localization
    %
    % Inputs:
    %   selectedSatPositions - 4x3 matrix of satellite positions
    %   TDOAs - 6x1 vector of TDOA measurements
    %   pairs - 6x2 matrix of satellite pairs
    %   covarianceMatrices - 3x3x6 tensor of covariance matrices for each pair
    %   initialGuess - Initial UE position [x, y, z]
    %   maxIter - Maximum number of iterations
    %   tol - Convergence tolerance
    %   c - Speed of light
    %
    % Outputs:
    %   p_est - Estimated UE position [x, y, z]
    %   residual - Final residual norm

    % Initialize variables
    p_est = initialGuess;
    residual = inf; 
    iter = 0; 
    numPairs = size(pairs, 1);

    while iter < maxIter && residual > tol
        iter = iter + 1;

        % Initialize Jacobian matrix and residual vector
        J = zeros(numPairs, 3);
        delta = zeros(numPairs, 1);

        % Compute Jacobian and residuals for each pair
        for k = 1:numPairs
            sat_i = selectedSatPositions(pairs(k, 1), :);
            sat_j = selectedSatPositions(pairs(k, 2), :);

            % Distances from current estimate to satellites
            d_i = norm(p_est - sat_i);
            d_j = norm(p_est - sat_j);

            % Jacobian row for pair (i, j)
            J(k, :) = (p_est - sat_i) / d_i - (p_est - sat_j) / d_j;

            % Residual for pair (i, j)
            delta(k) = (d_i - d_j) - c * TDOAs(k);
        end

        % Construct the weight matrix
        W = zeros(numPairs, numPairs);
        for k = 1:numPairs
            W(k, k) = 1 / trace(covarianceMatrices(:, :, k)); % Use covariance trace as weight
        end

        % Solve for position update using Weighted Least Squares
        delta_p = (J' * W * J) \ (J' * W * delta);

        % Update position estimate
        p_est = p_est - delta_p';
        residual = norm(delta);

        if residual < tol
            break;
        end
    end

    if residual > tol
        warning('TWLS did not converge within the maximum number of iterations.');
    end
end
