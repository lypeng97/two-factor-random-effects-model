function print_effect_correlations(alpha_FC, beta_FC, eta_FC, varpi_FC, ...
                                   alpha_SC, beta_SC, eta_SC, varpi_SC)
% PRINT_EFFECT_CORRELATIONS - Prints correlation between FC and SC effects
% Flips eta_FC and varpi_FC if rho_eta is negative

    rho_eta = corr(eta_FC(:), eta_SC(:));
    if rho_eta > 0
        eta_FC   = -eta_FC;
        varpi_FC = -varpi_FC;
    end

    rho_alpha  = corr(alpha_FC(:), alpha_SC(:));
    rho_beta   = corr(beta_FC(:), beta_SC(:));
    rho_eta    = corr(eta_FC(:), eta_SC(:));     % now guaranteed e 0
    rho_varpi  = corr(varpi_FC(:), varpi_SC(:)); % may now be negative

    fprintf('\nFC-SC Effect Correlations\n');
    fprintf('%s\n', repmat('=', 1, 60));
    fprintf('%-37s %s\n', 'Effect Type', 'Correlation');
    fprintf('%s\n', repmat('-', 1, 60));

    print_corr('Correlation between main edge effects', rho_alpha);
    print_corr('Correlation between main subject effects', rho_beta);
    print_corr('Correlation between interaction edge effects', rho_eta);
    print_corr('Correlation between interaction subject effects', rho_varpi);

    fprintf('%s\n\n', repmat('-', 1, 60));
end

function print_corr(label, val)
    if isnan(val)
        val_str = 'NaN';
    else
        val_str = sprintf('%8.4f', val);
    end
    fprintf('%-37s %10s\n', label, val_str);
end

