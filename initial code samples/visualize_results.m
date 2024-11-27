function visualize_results(satellite_positions, doppler_values, localization_error)
    % Visualize Satellite Orbits
    figure;
    hold on;
    num_sats = size(satellite_positions, 1);
    for sat = 1:num_sats
        plot3(satellite_positions(sat, :, 1), satellite_positions(sat, :, 2), ...
              satellite_positions(sat, :, 3), 'LineWidth', 1.5);
    end
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('Satellite Orbits');
    grid on;
    legend(arrayfun(@(x) sprintf('Sat %d', x), 1:num_sats, 'UniformOutput', false));
    hold off;

    % Visualize Doppler Shifts
    figure;
    plot(doppler_values', 'LineWidth', 1.5);
    xlabel('Time Step'); ylabel('Doppler Shift (Hz)');
    title('Doppler Shift for Each Satellite');
    legend(arrayfun(@(x) sprintf('Sat %d', x), 1:num_sats, 'UniformOutput', false));

    % Display Localization Error
    fprintf('Localization Error: %.2f meters\n', localization_error);
end
