# SimulationStudyVPC: Simulation Studies for Variance Partition Coefficient (VPC) in GLMMs

## Author and Contributors

### Author
- **Julius Olaifa** - [julius.olaifa@okstate.edu](mailto:julius.olaifa@okstate.edu)

### Contributor
- **Pratyaydipta Rudra** - Advisor and contributor

## Package Overview

**SimulationStudyVPC** is an R package that facilitates the simulation and comparison of Variance Partition Coefficients (VPC) in Generalized Linear Mixed Models (GLMMs). The package allows users to generate data from multiple model families (e.g., Negative Binomial, Tweedie, Gaussian) and compare the VPC estimates across these families using simulation studies.

The package provides utilities for:
- **Data generation**: Generating data under different GLMM family assumptions.
- **Model fitting**: Fitting GLMMs with user-defined formulas and models.
- **Parallel processing**: Efficiently simulating and estimating VPC across different model families using parallel computing.
- **VPC estimation**: Estimating the VPC for models to understand variance distribution in hierarchical or clustered data.

## Features

- **Multiple GLMM Families Supported**: Compare VPC estimates for models fit with families like Negative Binomial, Tweedie, and Gaussian.
- **Flexible Data Simulation**: Generate datasets under varying assumptions of the underlying model families and simulate different scenarios.
- **Reproducible Simulations**: Set seeds for reproducible simulation studies, ensuring consistent results across different systems and environments.
- **Parallel Computation**: Leverage multiple cores to speed up simulation and model fitting tasks.
- **Customizable Model Formulas**: Define and fit models with flexible formulas to capture relationships between variables in simulated or real-world data.

## Installation

To install the package from GitHub, run:

```r
# Install the development version from GitHub
devtools::install_github("juliusolaifa/SimulationStudyVPC")
