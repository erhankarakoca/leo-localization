clear all;
close all;
clc;
%% Create a satellite scenario
startTime = datetime(2020, 05, 04, 18,45,50);
stopTime = datetime(2020, 05, 04, 19,02,20);
sampleTime=10;
satscene = satelliteScenario(startTime,stopTime,sampleTime);

% Add satellites from TLE file.
tleFile = "leoSatelliteConstellation.tle";
constellation = satellite(satscene, tleFile);

% Define ue to be ground station Lat Long Alt;
ueStationLLA = [40.786648, 29.449502, 182];
gsUE = groundStation(satscene, ...
                     "Latitude",  ueStationLLA(1), ...
                     "Longitude", ueStationLLA(2), ...
                     "Altitude",  ueStationLLA(3));


%% Find the access intervals
ac = access(constellation,gsUE);
intvls = accessIntervals(ac);
% play(satscene)

%% Play the sat scene defined on tle with ue 
% play(satscene);

%% Any satellite can be selected from the viewer 
% Be just sure it it can access to the groun station (ue)
satNumber = 19;
referenceSatellite = "Satellite "+ string(satNumber) ;
rows = matches(intvls.Source, referenceSatellite);

% find specific acces-in times for the defined satellite
% startEndTimes = table2array(intvls(rows, ["StartTime", "EndTime"]));
% startTime = startEndTimes(1);
% endTime = startEndTimes(2);

%% Read tle to calculate specific satellite position and speed
% sampleTime = 60;
% sc = satelliteScenario(startTime,endTime,sampleTime);
% constellationInterval = satellite(sc, tleFile);
% gs = groundStation(sc, ueStation(1), ueStation(2));
% acInterval = access(constellationInterval,gs);

tleStruct = tleread('leoSatelliteConstellation.tle');
for i = 1:16
    orbitTime = startTime+minutes(i);
    [satPos(:,i),satVelocity(:,i)] = propagateOrbit(orbitTime, tleStruct(satNumber));
    [~,~,distanceSatToUe(i)] = aer(gsUE, constellation(satNumber), orbitTime);
end

toaSatToUe = distanceSatToUe / physconst("LightSpeed");

tdoa = toaSatToUe(end)-(toaSatToUe(1:end-1));


