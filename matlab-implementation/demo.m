% =======================================================
% DEMO: Two-Factor Decomposition on Example FC and SC Data
% =======================================================

clear; clc;

% Load example data: [subjects × regions × regions]
load('example_data_connectivity.mat');  % contains data_FC and data_SC

%% Step 1: Vectorize upper triangle for each subject
example_FC = vectorizeUpperTriangle(data_FC);  % size: [subjects × edges]
example_SC = vectorizeUpperTriangle(data_SC);  % size: [subjects × edges]

%% Step 2: FC Decomposition
FC = atanh(example_FC);
[mu_FC, alpha_FC, beta_FC, eta_FC, varpi_FC, residual_FC] = ...
    compute_mu_alpha_beta_eta_varpi(FC, 100, 1e-6);

% Print variance breakdown for FC
print_variance_decomposition(FC, alpha_FC, beta_FC, eta_FC, varpi_FC, residual_FC)

%% Step 3: SC Decomposition
SC = log(example_SC + 1);
[mu_SC, alpha_SC, beta_SC, eta_SC, varpi_SC, residual_SC] = ...
    compute_mu_alpha_beta_eta_varpi(SC, 100, 1e-6);

% Print variance breakdown for SC
print_variance_decomposition(SC, alpha_SC, beta_SC, eta_SC, varpi_SC, residual_SC)

%% Step 4: Correlation of Effects Between FC and SC
print_effect_correlations(alpha_FC, beta_FC, eta_FC, varpi_FC, ...
                          alpha_SC, beta_SC, eta_SC, varpi_SC);

