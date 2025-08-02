# Two-Factor Random Effects Modeling: Variance Decomposition of Connectivity Data

This project applies a two-factor random effects model with interaction terms to perform variance decomposition of brain connectivity matrices. The total variability in the observed data is decomposed into four components:

- **Main edge effects** ($`\alpha`$)
- **Main subject effects** ($`\beta`$)
- **Interaction effects** ($`\varpi \times \eta`$)
- **Residual error** ($`\epsilon`$)


---

## Implementations

This repository provides two implementations of the random effects model.

### 1. MATLAB Implementation

**Folder:** `matlab-implementation/`
**(requires MATLAB R2018 or later)**

This version uses **Newton–Raphson iterative optimization** to estimate interaction effects and supports:

* Full model fitting and variance decomposition
* Reporting of correlations between Functional Connectivity (FC) and Structural Connectivity (SC) model components.

**To run the analysis:**

1. Open MATLAB.
2. Navigate to the `matlab-implementation` folder.
3. Execute demo.m in the MATLAB command window.

### 2. Python Implementation

**Folder:** `python-implementation/`
**(Requires Python ≥ 3.7)**

This version employs **Singular Value Decomposition (SVD)** for efficient, low-rank estimation of interaction terms. This SVD-based solution is mathematically equivalent to the Newton-Raphson approach and is grounded in the [Eckart–Young–Mirsky theorem](https://en.wikipedia.org/wiki/Low-rank_approximation).

**To run the analysis:**
```bash
cd python-implementation
pip install -r requirements.txt
python demo.py
```

---

## Expected Output

After running either `demo.m` (MATLAB) or `demo.py` (Python), the following results will be generated:

- Percentage of variance explained by edge, subject, interaction, and residual effects
- Correlation between corresponding SC and FC effects, including:
  - $`\alpha`$ (edge effects)
  - $`\beta`$ (subject effects)
  - $`\eta`$ (interaction edge effect)
  - $`\varpi`$ (interaction subject effect)


---


## Files

**MATLAB Files (`matlab-implementation/` folder)**
| File                                | Description                                               |
|-------------------------------------|-----------------------------------------------------------|
| `compute_mu_alpha_beta_eta_varpi.m` | Core function for estimating model parameters.                           |
| `NR_optimization.m`                 | Newton–Raphson optimization for estimating $`\eta`$ and $`\varpi`$ interaction terms.  |
| `print_variance_decomposition.m`    | Outputs of the variance decomposition.                                  |
| `print_effect_correlations.m`       | Computation of correlations between Structural Connectivity (SC) and Functional Connectivity (FC) model components.        |
| `demo.m`                            | MATLAB demo script that runs the full analysis pipeline.                           |

**Python Files (`python-implementation/` folder)**
| File                                | Description                                       |
|-------------------------------------|---------------------------------------------------|
| `data_loading.py`                   | Functions for loading and preprocessing data.        |
| `model_estimation.py`               | Estimates of model parameters using the SVD-based approach.            |
| `variance_decomposition.py`         | Outputs of the variance decomposition.                            |
| `demo.py`                           | Python demo script that runs the full analysis pipeline.         |
| `requirements.txt`                  | Lists of required Python packages for installation.                       |

---



