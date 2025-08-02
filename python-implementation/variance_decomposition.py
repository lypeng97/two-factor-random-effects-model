"""
Variance decomposition for random effects models.
"""

import numpy as np


def compute_variance_components(model_results, connectivity_data):
    """Compute variance components from fitted model."""
    mu_hat = model_results['mu']
    alpha_hat = model_results['alpha']
    beta_hat = model_results['beta']
    rank1_interaction = model_results['rank1_interaction']
    epsilon_hat = model_results['epsilon']
    
    # Total variance
    centered_data = connectivity_data - mu_hat
    var_total = np.nanvar(centered_data, ddof=0)
    
    # Component variances
    var_alpha = np.nanvar(alpha_hat, ddof=0)
    var_beta = np.nanvar(beta_hat, ddof=0)
    var_interaction = np.nanvar(rank1_interaction, ddof=0)
    var_epsilon = np.nanvar(epsilon_hat, ddof=0)
    
    return {
        'total': var_total,
        'alpha': var_alpha,
        'beta': var_beta, 
        'interaction': var_interaction,
        'epsilon': var_epsilon
    }


def compute_variance_percentages(variance_components):
    """Convert variance components to percentages."""
    total_var = variance_components['total']
    
    if total_var < 1e-12:
        return {key: 0.0 for key in variance_components.keys()}
    
    return {
        'alpha': (variance_components['alpha'] / total_var) * 100.0,
        'beta': (variance_components['beta'] / total_var) * 100.0,
        'interaction': (variance_components['interaction'] / total_var) * 100.0,
        'epsilon': (variance_components['epsilon'] / total_var) * 100.0,
        'total': 100.0
    }


def print_variance_decomposition(variance_percentages, connectivity_type="Connectivity"):
    """Print formatted variance decomposition table."""
    print(f"\nVariance Decomposition - {connectivity_type}")
    print("=" * 55)
    print("Component                    % of Total Variance")
    print("-" * 55)
    
    component_labels = {
        'alpha': 'Main Edge Effects (α_ij)',
        'beta': 'Main Subject Effects (β^s)', 
        'interaction': 'Interaction Effects (η_ij·ϖ^s)',
        'epsilon': 'Residual Effects (ε_ij^s)'
    }
    
    for component in ['alpha', 'beta', 'interaction', 'epsilon']:
        if component in variance_percentages:
            label = component_labels[component]
            percentage = variance_percentages[component]
            print(f"{label:<28s} {percentage:>15.2f}%")
    
    print("-" * 55)
    print(f"{'Total':<28s} {variance_percentages['total']:>15.2f}%")
    print("")


def analyze_variance_decomposition(fc_results, sc_results, fc_data, sc_data):
    """Complete variance decomposition analysis."""
    
    fc_variance_components = compute_variance_components(fc_results, fc_data)
    sc_variance_components = compute_variance_components(sc_results, sc_data)
    
    fc_variance_percentages = compute_variance_percentages(fc_variance_components)
    sc_variance_percentages = compute_variance_percentages(sc_variance_components)
    
    print_variance_decomposition(sc_variance_percentages, "Structural Connectivity (SC)")
    print_variance_decomposition(fc_variance_percentages, "Functional Connectivity (FC)")
    
    return {
        'fc': {
            'components': fc_variance_components,
            'percentages': fc_variance_percentages
        },
        'sc': {
            'components': sc_variance_components, 
            'percentages': sc_variance_percentages
        }
    }