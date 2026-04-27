clear; close all; clc;

addpath('src');

%% Setup model

model = setup_model();

%% Simulate trajectory

m_sim = 300;
X_sim = simulate_trajectory(model, m_sim);

figure;
plot(X_sim(1,:), X_sim(4,:), 'b-', 'LineWidth', 1.5); hold on;
xlabel('x_1');
ylabel('x_2');
title('Simulated Trajectory');
grid on;
axis equal;

%% Load data

S = load('data/stations.mat');
pos_vec = S.pos_vec;

R = load('data/RSSI-measurements.mat');
Y = extract_measurements(R);

%% SIS

N_sis = 10000;
hist_times = [0 5 10 25 50];

sis = run_sis(model, Y, pos_vec, N_sis, hist_times);

plot_results('sis', sis, pos_vec, hist_times);

disp('Selected ESS values:')
disp(table(hist_times(:), sis.ESS(hist_times(:)+1), ...
    'VariableNames', {'n','ESS'}))

%% SISR

N_sisr = 10000;

sisr = run_sisr(model, Y, pos_vec, N_sisr);

plot_results('sisr', sisr, pos_vec, []);

%% Estimate unknown observation noise standard deviation

R_unknown = load('data/RSSI-measurements-unknown-sigma.mat');
Y_unknown = extract_measurements(R_unknown);

varsigma_grid = 0.2:0.1:2.9;
N_calib = 5000;

calib = estimate_varsigma_grid(model, Y_unknown, pos_vec, varsigma_grid, N_calib);

disp('Estimated observation noise standard deviation:')
disp(calib.varsigma_hat)

plot_results('calibration', calib, pos_vec, []);