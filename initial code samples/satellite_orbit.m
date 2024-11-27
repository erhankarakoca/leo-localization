function positions = satellite_orbit(num_satellites, altitude, sim_time, dt)
    % Calculate satellite positions over time
    earth_radius = 6371e3; % Earth radius in meters
    orbit_radius = earth_radius + altitude;
    time_steps = 0:dt:sim_time;

    positions = zeros(num_satellites, length(time_steps), 3);
    for sat = 1:num_satellites
        angle_offset = 2 * pi * (sat - 1) / num_satellites; % Evenly spaced
        for t = 1:length(time_steps)
            theta = 2 * pi * time_steps(t) / (90 * 60) + angle_offset; % 90-min orbit
            positions(sat, t, :) = [orbit_radius * cos(theta), ...
                                    orbit_radius * sin(theta), ...
                                    0]; % Assuming equatorial orbit
        end
    end
end
