"""
Random effects models parameter estimation.
"""

import numpy as np


def estimate_global_mean(connectivity_data):
    """Estimate global mean μ."""
    return np.nanmean(connectivity_data)


def estimate_main_edge_effects(connectivity_data, mu_hat):
    """Estimate main edge effects α_{ij}."""
    return np.nanmean(connectivity_data, axis=0) - mu_hat


def estimate_main_subject_effects(connectivity_data, mu_hat):
    """Estimate main subject effects β^s."""
    return np.nanmean(connectivity_data, axis=1) - mu_hat


def enforce_identifiability_constraints(eta, varpi, n_edges, center_effects=True, normalize_eta=True):
    """Enforce identifiability constraints and variance normalization."""
    eta_constrained = eta.copy()
    varpi_constrained = varpi.copy()
    
    # Constraints
    if center_effects:
        eta_constrained = eta_constrained - np.mean(eta_constrained)
        varpi_constrained = varpi_constrained - np.mean(varpi_constrained)
    
    # Variance normalization
    if normalize_eta:
        eta_squared_sum = np.sum(eta_constrained**2)
        if eta_squared_sum > 1e-12:
            scale_factor = np.sqrt(n_edges / eta_squared_sum)
            eta_constrained = eta_constrained * scale_factor
            varpi_constrained = varpi_constrained / scale_factor
    
    return eta_constrained, varpi_constrained


def estimate_interaction_effects_svd(residual_matrix, n_edges, enforce_constraints=True):
    """Estimate interaction effects η_{ij}, ϖ^s using SVD rank-1 approximation."""
    # SVD decomposition
    u, singular_values, vh = np.linalg.svd(residual_matrix, full_matrices=False)
    
    s1 = singular_values[0]
    u1 = u[:, 0] 
    v1 = vh[0, :]     
    
    varpi_init = np.sqrt(s1) * u1  
    eta_init = np.sqrt(s1) * v1 
    
    # Apply constraints
    if enforce_constraints:
        eta_hat, varpi_hat = enforce_identifiability_constraints(
            eta_init, varpi_init, n_edges
        )
    else:
        eta_hat, varpi_hat = eta_init, varpi_init
    
    rank1_approximation = np.outer(varpi_hat, eta_hat)
    
    return eta_hat, varpi_hat, rank1_approximation


def fit_random_effects_model(connectivity_data, enforce_constraints=True):
    """Fit complete random effects model."""
    n_subjects, n_edges = connectivity_data.shape
    
    mu_hat = estimate_global_mean(connectivity_data)
    alpha_hat = estimate_main_edge_effects(connectivity_data, mu_hat)
    beta_hat = estimate_main_subject_effects(connectivity_data, mu_hat)
    
    additive_effects = mu_hat + alpha_hat[None, :] + beta_hat[:, None]
    residuals_additive = connectivity_data - additive_effects

    eta_hat, varpi_hat, rank1_approximation = estimate_interaction_effects_svd(
        residuals_additive, n_edges, enforce_constraints
    )
    
    epsilon_hat = residuals_additive - rank1_approximation
    fitted_values = additive_effects + rank1_approximation
    
    return {
        'mu': mu_hat,
        'alpha': alpha_hat,
        'beta': beta_hat, 
        'eta': eta_hat,
        'varpi': varpi_hat,
        'epsilon': epsilon_hat,
        'fitted_values': fitted_values,
        'rank1_interaction': rank1_approximation,
        'additive_effects': additive_effects
    }


def ensure_positive_correlation_varpi(fc_results, sc_results):
    """Ensure positive correlation between FC and SC interaction subject effects for identifiability (Optional)."""
    varpi_fc = fc_results['varpi']
    varpi_sc = sc_results['varpi']
    
    valid_idx = ~(np.isnan(varpi_fc) | np.isnan(varpi_sc))
    if np.sum(valid_idx) < 2:
        return 0.0
    
    rho_varpi = np.corrcoef(varpi_fc[valid_idx], varpi_sc[valid_idx])[0, 1]
    if rho_varpi < 0:
        fc_results['eta'] = -fc_results['eta']
        fc_results['varpi'] = -fc_results['varpi']
        fc_results['rank1_interaction'] = -fc_results['rank1_interaction']
        
        rho_varpi = np.corrcoef(
            fc_results['varpi'][valid_idx], 
            sc_results['varpi'][valid_idx]
        )[0, 1]
        print(f"After adjustment: ρ_ϖ = {rho_varpi:.4f}")
    
    return rho_varpi


def compute_effect_correlations(fc_results, sc_results, adjust_signs=True):
    """Compute correlations between corresponding FC and SC effects."""
    def safe_correlation(x, y, eps=1e-12):
        valid_mask = ~(np.isnan(x) | np.isnan(y))
        if np.sum(valid_mask) < 2:
            return np.nan
        x_valid, y_valid = x[valid_mask], y[valid_mask]
        if np.std(x_valid) < eps or np.std(y_valid) < eps:
            return np.nan
        return np.corrcoef(x_valid, y_valid)[0, 1]
    
    if adjust_signs:
        ensure_positive_correlation_varpi(fc_results, sc_results)
    
    correlations = {
        'rho_alpha': safe_correlation(fc_results['alpha'], sc_results['alpha']),
        'rho_beta': safe_correlation(fc_results['beta'], sc_results['beta']),
        'rho_eta': safe_correlation(fc_results['eta'], sc_results['eta']),
        'rho_varpi': safe_correlation(fc_results['varpi'], sc_results['varpi'])
    }
    
    return correlations