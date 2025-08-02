function print_variance_decomposition(data, alpha, beta, eta, varpi, residual)
    % PRINT_VARIANCE_DECOMPOSITION - Print sum of squares decomposition and proportions
    %
    % Inputs:
    %   data     - Original matrix (subjects × edges)
    %   alpha    - Edge effects (1 × edges)
    %   beta     - Subject effects (subjects × 1)
    %   eta      - Interaction vector for edges (1 × edges)
    %   varpi    - Interaction vector for subjects (subjects × 1)
    %   residual - Residual matrix (subjects × edges)
    %
    % Output:
    %   Prints variance decomposition (SSA, SSB, SSI, SSE, SST) and their proportions

    % Step 1: Total sum of squares
    grand_mean = mean(data(:));
    SST = sum((data(:) - grand_mean).^2);

    % Step 2: SSA (Edge/Alpha effects)
    SSA = sum((alpha(:) - mean(alpha(:))).^2) * size(data, 1);

    % Step 3: SSB (Subject/Beta effects)
    SSB = sum((beta(:) - mean(beta(:))).^2) * size(data, 2);

    % Step 4: SSI (Interaction effect: varpi × eta)
    interaction_matrix = varpi * eta;
    interaction_mean = mean(interaction_matrix(:));
    SSI = sum((interaction_matrix(:) - interaction_mean).^2);

    % Step 5: SSE (Residual)
    residual_centered = residual - mean(residual(:));
    SSE = sum(residual_centered(:).^2);

    % Step 6: Proportions
    prop_SSA = SSA / SST;
    prop_SSB = SSB / SST;
    prop_SSI = SSI / SST;
    prop_SSE = SSE / SST;

    % Step 7: Print results
    fprintf('\n========== Variance Decomposition ==========\n');
    fprintf('Total Sum of Squares (SST):       %.4f\n', SST);
    fprintf('SSA (Edge/Alpha Effect):          %.4f (%.2f%%)\n', SSA, 100 * prop_SSA);
    fprintf('SSB (Subject/Beta Effect):        %.4f (%.2f%%)\n', SSB, 100 * prop_SSB);
    fprintf('SSI (Interaction: Varpi × Eta):   %.4f (%.2f%%)\n', SSI, 100 * prop_SSI);
    fprintf('SSE (Residual/Error):             %.4f (%.2f%%)\n', SSE, 100 * prop_SSE);
    fprintf('--------------------------------------------\n');
    fprintf('Total Explained Variance:         %.2f%%\n', 100 * (prop_SSA + prop_SSB + prop_SSI));
    fprintf('Residual Unexplained Variance:    %.2f%%\n', 100 * prop_SSE);
    fprintf('============================================\n\n');
end
