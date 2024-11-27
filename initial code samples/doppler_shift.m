function doppler = doppler_shift(sat_positions, ue_position, carrier_freq, dt)
    c = 3e8; % Speed of light
    num_sats = size(sat_positions, 1);
    time_steps = size(sat_positions, 2);

    doppler = zeros(num_sats, time_steps);
    for sat = 1:num_sats
        for t = 1:time_steps - 1
            sat_pos = squeeze(sat_positions(sat, t, :))';
            sat_next_pos = squeeze(sat_positions(sat, t+1, :))';
            velocity = (sat_next_pos - sat_pos) / dt;

            relative_velocity = dot(velocity, (sat_pos - ue_position)) / norm(sat_pos - ue_position);
            doppler(sat, t) = (relative_velocity / c) * carrier_freq;
        end
    end
end
