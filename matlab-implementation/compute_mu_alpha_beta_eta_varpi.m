function [mu, alpha, beta, eta, varpi, residual] = compute_mu_alpha_beta_eta_varpi(data, max_iter, tol)
    % COMPUTE_MU_ALPHA_BETA_ETA_VARPI estimates main and interaction effects
    %
    %   [mu, alpha, beta, eta, varpi, residual] = compute_mu_alpha_beta_eta_varpi(data, max_iter, tol)
    %
    %   Input:
    %     data        - Input connectivity matrix (subjects × edges or regions)
    %     max_iter    - Maximum number of iterations for optimization
    %     tol         - Convergence tolerance
    %
    %   Output:
    %     mu          - Grand mean (scalar)
    %     alpha       - Main edge effects (1 × edges)
    %     beta        - Main subject effects (subjects × 1)
    %     eta         - Scaled interaction vector for edges (1 × edges)
    %     varpi       - Scaled interaction vector for subjects (subjects × 1)
    %     residual    - Residual term (subjects × edges)

    % Step 1: Compute main effects
    mu = mean(data, 'all');
    alpha = mean(data - mu, 1);
    beta = mean(data - mu - alpha, 2);
    epsilon = data - mu - alpha - beta;

    % Step 2: Random initialization of interaction terms
    a = 0.01; 
    b = 5;
    variance = a + (b - a) * rand(2, 1);  % [var_eta; var_varpi]
    eta_init = variance(1) * randn(1, size(data, 2));
    varpi_init = variance(2) * randn(size(data, 1), 1);

    % Step 3: Optimize interaction terms
    result = NR_optimization(epsilon, eta_init, varpi_init, max_iter, tol);

    % Step 4: Rescale interaction terms
    eta_raw = result.final_eta;
    varpi_raw = result.final_varpi;
    
    eta_std = std(eta_raw(:));
    eta = (eta_raw - mean(eta_raw(:))) / eta_std;
    varpi = varpi_raw * eta_std;

    % Step 5: Compute residual (data - predicted)
    predicted = mu + alpha + beta + varpi * eta;
    residual = data - predicted;
end
