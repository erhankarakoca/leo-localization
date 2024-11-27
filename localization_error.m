function error = localization_error(true_position, estimated_position)
    % Compute localization error (Euclidean distance)
    error = norm(true_position - estimated_position);
end
