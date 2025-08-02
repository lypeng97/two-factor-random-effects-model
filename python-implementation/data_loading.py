"""
Data loading and preprocessing for random effects models.
"""

import numpy as np
from scipy.io import loadmat

def load_example_data(example_data_path='example_data_connectivity.mat'):
    """Load example connectivity data for demonstration."""
    data = loadmat(example_data_path)
    
    # Extract FC and SC data (both 300×50×50: subjects × nodes × nodes)
    fc_data = data['data_FC']  
    sc_data = data['data_SC']
    
    print(f"Loaded example data:")
    print(f"  FC shape: {fc_data.shape}")
    print(f"  SC shape: {sc_data.shape}")
    
    fc_data = fc_data.astype(np.float32)
    sc_data = sc_data.astype(np.float32)
    
    # Generate simple subject IDs
    n_subjects = fc_data.shape[0]
    subject_ids = [f"subj_{i+1:03d}" for i in range(n_subjects)]
    
    return sc_data, fc_data, subject_ids

def load_fc_ids_from_mat_numeric(fc_id_mat_path, var_name='id_list'):
    """Extract subject IDs"""
    fc_id_data = loadmat(fc_id_mat_path)
    if var_name not in fc_id_data:
        raise KeyError(f"Could not find '{var_name}' in the .mat file.")
    
    id_array = fc_id_data[var_name]
    
    if id_array.ndim == 2 and id_array.shape[0] == 1:
        id_array = id_array[0]
    elif id_array.ndim == 2 and id_array.shape[1] == 1:
        id_array = id_array[:, 0]
    
    fc_ids = [str(int(x)) for x in id_array]
    return fc_ids


def load_connectivity_data(sc_data, fc_data, fc_id_mat_path):
    """Load and align FC and SC connectivity data across subjects."""
    # Extract SC subject IDs
    sc_ids = sc_data['ID']
    if sc_ids.ndim == 2 and sc_ids.shape[0] == 1:
        sc_ids = sc_ids[0]
    sc_ids = [str(x).strip() for x in sc_ids]

    sc_count_all = sc_data['SC_Count']
    sc_fa_all = sc_data['SC_FA']
    print(f"Found {len(sc_ids)} subjects in SC data.")

    # Extract FC data
    fc_array = fc_data['FC']
    n_fc_subjects = fc_array.shape[2]
    print(f"FC array shape: {fc_array.shape} => {n_fc_subjects} FC subjects.")

    fc_ids = load_fc_ids_from_mat_numeric(fc_id_mat_path, var_name='id_list')
    print(f"Found {len(fc_ids)} IDs in FC subject list.")
    
    if len(fc_ids) != n_fc_subjects:
        print("WARNING: Number of FC IDs does not match FC array dimensions.")

    fc_id_to_idx = {subj_id: idx for idx, subj_id in enumerate(fc_ids)}

    sc_count_matrices = []
    sc_fa_matrices = []
    fc_matrices = []
    aligned_subject_ids = []
    missing_fc_subjects = []

    for i, subject_id in enumerate(sc_ids):
        try:
            sc_count_mat = np.array(sc_count_all[i, :, :], dtype=np.float32)
            sc_fa_mat = np.array(sc_fa_all[i, :, :], dtype=np.float32)
            
            if subject_id in fc_id_to_idx:
                fc_idx = fc_id_to_idx[subject_id]
                fc_mat = np.array(fc_array[:, :, fc_idx], dtype=np.float32)
                
                sc_count_matrices.append(sc_count_mat)
                sc_fa_matrices.append(sc_fa_mat)
                fc_matrices.append(fc_mat)
                aligned_subject_ids.append(subject_id)
            else:
                missing_fc_subjects.append(subject_id)
                
        except Exception as e:
            print(f"Error processing subject {subject_id}: {e}")
            continue

    print(f"Successfully aligned {len(aligned_subject_ids)} subjects.")

    return (np.array(sc_count_matrices), 
            np.array(sc_fa_matrices), 
            np.array(fc_matrices), 
            aligned_subject_ids)


def preprocess_connectivity_matrix(matrix):
    """Extract upper triangular values from symmetric connectivity matrix."""
    matrix_copy = np.array(matrix, copy=True)
    np.fill_diagonal(matrix_copy, 0)  # Zero diagonal
    upper_tri_indices = np.triu_indices_from(matrix_copy, k=1)
    return matrix_copy[upper_tri_indices]


def flatten_connectivity_arrays(connectivity_array):
    """Flatten connectivity arrays"""
    n_subjects, n_nodes, _ = connectivity_array.shape
    n_edges = n_nodes * (n_nodes - 1) // 2
    
    flattened_data = np.full((n_subjects, n_edges), np.nan)
    
    for subject_idx in range(n_subjects):
        connectivity_matrix = connectivity_array[subject_idx]
        edge_values = preprocess_connectivity_matrix(connectivity_matrix)
        flattened_data[subject_idx, :] = edge_values
    
    return flattened_data, n_edges


def apply_transformations(sc_data, fc_data, apply_log_transform=True, apply_fisher_z=True):
    """Apply log (for SC) and Fisher-z (for FC) transformations to connectivity data."""
    transformed_sc = sc_data.copy()
    transformed_fc = fc_data.copy()
    
    if apply_log_transform:
        transformed_sc = np.log1p(transformed_sc) 
    
    if apply_fisher_z:
        eps = 1e-8
        transformed_fc = np.clip(transformed_fc, -1 + eps, 1 - eps)
        transformed_fc = np.arctanh(transformed_fc)
    
    return transformed_sc, transformed_fc