clear; clc; close all;

% Simulation Parameters
num_satellites = 10; % Number of satellites
orbit_altitude = 550e3; % Altitude in meters
carrier_frequency = 5e9; % Frequency in Hz
ue_position = [0, 0, 0]; % UE position (x, y, z) in meters
simulation_time = 50; % Simulation duration in seconds
time_step = 0.1; % Time step in seconds

% Generate Satellite Orbits
satellite_positions = satellite_orbit(num_satellites, orbit_altitude, simulation_time, time_step);

% Calculate Doppler Shift for Each Satellite
doppler_values = doppler_shift(satellite_positions, ue_position, carrier_frequency, time_step);

% Compute Time-of-Arrival (TOA)
toa_values = toa_calculation(satellite_positions, ue_position, time_step);

% Estimate Localization and Error
localization_error = localization_error(toa_values, doppler_values);

% Visualize Results
visualize_results(satellite_positions, doppler_values, localization_error);
