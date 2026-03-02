# Bayesian Homework

This repository contains two Bayesian statistics assignments, each organized as a self-contained homework folder with its own write-up, source data, and generated artifacts.

## Repository Layout

### `hw1/`

Assignment 1 focuses on introductory Bayesian modeling:

- Problem 1 derives Beta posteriors for two movie-review proportions under uniform priors and compares posterior summaries.
- Problem 2 analyzes Wikipedia article lengths, applies a log transform, and performs Bayesian inference for a normal model.

Key files:

- `hw1/README.md`: rendered assignment answers
- `hw1/main.Rmd`: source notebook
- `hw1/Random Wikipedia.txt`: input data
- `hw1/figure/`: generated plots

### `hw2/`

Assignment 2 focuses on hierarchical Bayesian models and JAGS:

- Problem 1 explores hyperpriors for Beta-distributed latent parameters through simulation.
- Problem 2 builds and fits a normal-normal hierarchical model for heart study data using JAGS.

Key files:

- `hw2/README.md`: rendered assignment answers
- `hw2/main.Rmd`: source notebook
- `hw2/Heart Studies.txt`: input data
- `hw2/asgn2template.bug`: JAGS model file
- `hw2/updated_asgn2template.bug`, `hw2/indicator_asgn2template.bug`: model variants
- `hw2/main.pdf`: rendered PDF output
- `hw2/figure/`: generated plots
- `hw2/dag.png`, `hw2/dag2.png`: model diagrams

## How To Read This Repo

Start with the assignment-specific README in each folder:

- `hw1/README.md`
- `hw2/README.md`

If you want the source used to generate the write-ups, open the corresponding `main.Rmd` file in each directory.
