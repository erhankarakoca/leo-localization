function [estUEPosTDOAError, localizationErrorNLS] = NLSoptimizer(selectedSatPositions, pairs, radiiDifferenceswithError, initialGuess, lowerBound, upperBound, actualUEPosition)
    % estimateTDOALocalization
    % Estimates the UE position using TDOA measurements and nonlinear least squares.
    %
    % Inputs:
    %   selectedSatPositions - Nx3 matrix of satellite positions
    %   pairs - Px2 matrix of satellite index pairs
    %   radiiDifferenceswithError - Px1 vector of TDOA radii differences with noise
    %   initialGuess - 1x3 vector of the initial guess for UE position [x, y, z]
    %   lowerBound - 1x3 vector of lower bounds for UE position [x_min, y_min, z_min]
    %   upperBound - 1x3 vector of upper bounds for UE position [x_max, y_max, z_max]
    %   actualUEPosition - 1x3 vector of the true UE position [x, y, z] (for error calculation)
    %
    % Outputs:
    %   estUEPosTDOAError - Estimated UE position [x, y, z]
    %   localizationErrorNLS - Localization error (Euclidean distance between estimated and actual positions)

    % Define the objective function for TDOA estimation error
    function error = objectiveTDOA(p)
        % Compute the residuals for all TDOA pairs
        numPairs = size(pairs, 1);
        error = zeros(numPairs, 1);
        for k = 1:numPairs
            % Get the satellite indices for the current pair
            sat1 = pairs(k, 1);
            sat2 = pairs(k, 2);
            
            % Compute the distance differences for the current pair
            dist1 = sqrt((p(1) - selectedSatPositions(sat1, 1))^2 + ...
                         (p(2) - selectedSatPositions(sat1, 2))^2 + ...
                         (p(3) - selectedSatPositions(sat1, 3))^2);
            dist2 = sqrt((p(1) - selectedSatPositions(sat2, 1))^2 + ...
                         (p(2) - selectedSatPositions(sat2, 2))^2 + ...
                         (p(3) - selectedSatPositions(sat2, 3))^2);
            
            % Calculate the residual for this pair
            error(k) = abs(dist1 - dist2 - radiiDifferenceswithError(k));
        end
    end

    % Optimization options for lsqnonlin
    options = optimoptions('lsqnonlin');

    % Solve the nonlinear least squares problem
    estUEPosTDOAError = lsqnonlin(@objectiveTDOA, initialGuess, lowerBound, upperBound, options);

    % Calculate the localization error (distance to the true position)
    localizationErrorNLS = norm(estUEPosTDOAError - actualUEPosition);
end
