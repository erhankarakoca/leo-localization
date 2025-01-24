% Define the satellite scenario
startTime = datetime(2025, 1, 20, 0, 0, 0);
stopTime = startTime + hours(1); % Simulate for 1 hour
sampleTime = 10; % Time step in seconds
sc = satelliteScenario(startTime, stopTime, sampleTime);

% Define constellation parameters
altitude = 1200e3; % Altitude in meters
inclination = 88; % Near-polar orbit for global coverage
numPlanes = 36; % Number of orbital planes
satsPerPlane = 20; % Number of satellites per plane
raanSeparation = 360 / numPlanes; % RAAN separation between planes
phasingOffset = 360 / satsPerPlane; % Phasing offset within each plane

% Add satellites to the scenario
for plane = 0:(numPlanes-1)
    raan = plane * raanSeparation; % RAAN for the current plane
    for sat = 0:(satsPerPlane-1)
        % Compute phasing of satellites in the plane
        meanAnomaly = mod(sat * phasingOffset, 360); % Offset in degrees

        % Orbital elements
        semiMajorAxis = 6371e3 + altitude; % Earth's radius + altitude
        eccentricity = 0; % Circular orbit
        inclinationRad = deg2rad(inclination); % Convert inclination to radians
        raanRad = deg2rad(raan); % Convert RAAN to radians
        meanAnomalyRad = deg2rad(meanAnomaly); % Convert mean anomaly to radians
        argumentOfPeriapsis = 0; % Not relevant for circular orbits

        % Add satellite using Keplerian elements
        satellite(sc, semiMajorAxis, eccentricity, inclinationRad, ...
                  raanRad, argumentOfPeriapsis, meanAnomalyRad);
    end
end

% Visualize the constellation
viewer = satelliteScenarioViewer(sc);
play(sc);
