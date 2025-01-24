%% Create a satellite scenario
% startTime = datetime(2020, 05, 04, 18,45,50);
% stopTime = datetime(2020, 05, 04, 19,02,20);
sampleTime = 10;
satscene = satelliteScenario('AutoSimulate', true);

% Add satellites from TLE file.
tleFile = "starlink.tle";
constellation = satellite(satscene, "allLeoSatellites.tle");
% numSatellites = numel(constellation);
ueStationLLA = [40.786648, 29.449502, 182];
ueStationECEF = lla2ecef(ueStationLLA);

minElevationAngle = 30; % degrees
gsUE = groundStation(satscene, ...
                     "Latitude",  ueStationLLA(1), ...
                     "Longitude", ueStationLLA(2), ...
                     "Altitude",  ueStationLLA(3), ...
                     MinElevationAngle=minElevationAngle);
ac = access(constellation,gsUE);

c = physconst("LightSpeed");
v = satelliteScenarioViewer(satscene,'ShowDetails',false);

play(satscene)