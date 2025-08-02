"""
Demo script for FC-SC random effects model analysis.

Usage: python demo.py
"""

import numpy as np
from scipy.io import loadmat

from data_loading import (
    load_example_data, 
    load_connectivity_data, 
    flatten_connectivity_arrays, 
    apply_transformations
)
from model_estimation import (
    fit_random_effects_model,
    compute_effect_correlations
)
from variance_decomposition import analyze_variance_decomposition


def print_effect_correlations(correlations):
    """Print FC-SC effect correlations."""
    print("\nFC-SC Effect Correlations")
    print("=" * 60)
    print("Effect Type                           Correlation")
    print("-" * 60)
    
    correlation_info = [
        ('rho_alpha', 'Correlation between main edge effects'),
        ('rho_beta', 'Correlation between main subject effects'),
        ('rho_eta', 'Correlation between interaction edge effects'), 
        ('rho_varpi', 'Correlation between interaction subject effects')
    ]
    
    for corr_key, description in correlation_info:
        if corr_key in correlations:
            correlation_value = correlations[corr_key]
            corr_str = "NaN" if np.isnan(correlation_value) else f"{correlation_value:>8.4f}"
            print(f"{description:<37s} {corr_str:>10s}")
    
    print("-" * 60)
    print("")



def main():
    """Main analysis pipeline."""
    print("FC-SC Random Effects Model Analysis")
    print("=" * 70)
    
    # Choose data loading method
    use_example_data = True  # Use Random Generated Example Data
    
    if use_example_data:
        # Simple example data loading
        print("Loading example data...")
        sc_matrices, fc_matrices, subject_ids = load_example_data()
        
    else:
        # Full data pipeline (update these paths for your actual data)
        print("Loading full dataset...")
        fc_id_mat_path = "Folder/id_list.mat"
        det_mat_path = 'Folder/SC_SD_Stream_Schaefer200_data.mat'
        fc_mat_path = 'Folder/0626_Schaefer_FC.mat'
        
        det_mat = loadmat(det_mat_path)
        fc_mat = loadmat(fc_mat_path)
        
        sc_matrices, _, fc_matrices, subject_ids = load_connectivity_data(
            sc_data=det_mat, fc_data=fc_mat, fc_id_mat_path=fc_id_mat_path
        )
    
    n_subjects = len(subject_ids)
    n_nodes = sc_matrices.shape[1] if len(sc_matrices.shape) == 3 else sc_matrices.shape[2]
    n_edges = n_nodes * (n_nodes - 1) // 2
    print(f"Processed {n_subjects} subjects, {n_nodes} nodes, {n_edges} edges")
    

    sc_flat, _ = flatten_connectivity_arrays(sc_matrices)
    fc_flat, _ = flatten_connectivity_arrays(fc_matrices)
    sc_log, fc_fisherz = apply_transformations(sc_flat, fc_flat)
    
    # Fit models
    print("Fitting models...")
    sc_results = fit_random_effects_model(sc_log)
    fc_results = fit_random_effects_model(fc_fisherz)
    
    # Variance decomposition
    variance_results = analyze_variance_decomposition(
        fc_results, sc_results, fc_fisherz, sc_log
    )
    
    # Effect correlations
    correlations = compute_effect_correlations(fc_results, sc_results)
    print_effect_correlations(correlations)


if __name__ == "__main__":
    main()

