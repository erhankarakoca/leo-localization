function toa = toa_calculation(sat_positions, receivers, dt)
    % Compute Time of Arrival (TOA) for signals from satellites to receivers
    c = 3e8; % Speed of light in m/s
    num_sats = size(sat_positions, 1);
    num_receivers = size(receivers, 1);
    num_time_steps = size(sat_positions, 2);

    toa = zeros(num_receivers, num_time_steps);

    for t = 1:num_time_steps
        for r = 1:num_receivers
            for sat = 1:num_sats
                % Calculate distance between satellite and receiver
                sat_pos = squeeze(sat_positions(sat, t, :))';
                distance = norm(sat_pos - receivers(r, :));
                
                % Calculate TOA
                toa(r, t) = distance / c;
            end
        end
    end
end
