# Bayesian Nonparametric (BNP) Copula Project

This repository contains R scripts and structured analyses for Bayesian Nonparametric (BNP) copula-based modeling, clustering, and inference on real and simulated datasets using various copula families.

## Directory Structure

- **Data/**: Includes input datasets for real-world and simulated scenarios.
- **Code/**: Contains R scripts for simulation studies, MCMC sampling, and copula density calculations.
- **Result/**: Stores outputs of MCMC samples, estimated parameters, and analysis results.

### Datasets
- **Real Data**: 
  - Wine dataset (`wine_red_rank.rds`)
  - Occupational data (`occupDIFrank.csv`)
- **Simulated Data**: Generated Gaussian Copula data (`mix.rds`, `single.rds`) by using `simulate_gaussian_copula.R`

### Code Folder (`Code/`)

#### Helper functions for running MCMC
- `density.R`: Density and generator functions for Archimedean copulas.
- `post_sampling_helper.R`: Post-processing functions for MCMC results and cluster analysis.
- `sampling_helper.R`: Utility functions for sampling priors and hyperparameter updates.

#### Main Scripts for running MCMC
- `simulation_study1.R`, `simulation_study2.R`, `simulation_study3.R`: Scripts conducting simulation studies evaluating various copula models.
- `real_occup.R`: Copula modeling for occupational dataset.
- `real_wine.R`: Copula modeling for wine dataset.
- `Gaussian_copula.R`: Functions specific to Gaussian copula inference using Metropolis-Hastings.

### Result and Corresponding Analysis Scripts

#### Result
- **Simulation500/**: Simulation study 1 with sample size 500.
- **Simulation200/**: Simulation study 1 with sample size 200.
- **Gaussian_single/**: Simulation study 2 using single-component Gaussian copula simulations.
- **Gaussian_mix/**: Simulation study 2 using mixture Gaussian copula simulations.
- **Clayton_Hdim/**: Simulation study 3 using 4-dimensional Clayton copula simulations.
- **real_wine/**: Analysis and results related to the wine dataset.
- **Real_occup/**: Analysis and results related to occupational data.
- **Figures/**: Visualizations and graphical outputs.

#### Analysis Scripts
- `res_helper.R`: General helper functions for result analysis.
- `simulation_study1_result_analysis.R`: Analysis for simulation study 1.
- `simulation_study2_result_analysis.R`: Analysis for simulation study 2.
- `simulation_study3_result_analysis.R`: Analysis for simulation study 3.
- `real_occupancy_result_analysis.R`: Analysis of occupational data results.
- `red_wine_result_analysis.R`: Analysis of wine data results.

## Workflow

The typical analysis workflow includes:

1. Data loading and preprocessing.
2. Parameter and hyperparameter initialization.
3. MCMC sampling with adaptive tuning.
4. Post-processing of samples and posterior inference.
5. Visualization and result analyses.

## Dependencies
- R packages: `copula`, `foreach`, `doParallel`, `salso`, `tidyverse`

## Running the Analysis

Run scripts via command line:
```bash
Rscript Code/real_wine.R job_name job_num path
```

## Results

Analyses yield structured results:
- Parameter samples and credible intervals.
- Acceptance rates and convergence diagnostics.
- LPML and other model comparison metrics.

---
Refer to the specific scripts in the **Code** directory and respective analysis files for detailed documentation.