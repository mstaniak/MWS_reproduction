# MSstatsWeightedSummary: reproduction code

This repository contains code required to reproduce the evaluation of weighted protein summarization method described in an article "Relative quantification of proteins and post-translational modifications in proteomic experiments with shared peptides: a weight-based approach".
Due to GitHub's file size limits and data privacy, input data necessary for running the code are only available on request from the owner of the repository.

After obtaining the folder with all input files, it is possible to reproduce the analyses presented in the paper by running the following scripts:

- 00_setup.R to install necessary packages, including two GitHub-based libraries: MSstatsWeightedSummary (which implements the weighted summarization method) and SimulateTMT (which implements model-based protein- and peptide-level data simulation),
- 01_graphs_explanation.R to reproduce peptide-protein graphs and profiles plots of the protein degrader case study,
- brd_explanation_plot.R to reproduce the plot for the toy example showcasing weight-based summarization,
- brd_simulations_main.R to reproduce all feature resampling-based simulations presented in the main article,
- brd_simulations_si.R to reproduce all feature resampling-based simulations presented in the Supplement,
- ptm_all.R to reproduce the results related to the PTM case study,
- simulated_data.R to reproduce the model-based simulations,
- simulated_data_covariance.R to reproduce the plot which describes the relationship between variances of different parameter estimates based on a simple cluster,
- thermal_cluster_plot.R to reproduce analyses of thermal profiling data,
- thermal_sims.R to reproduce results of feature resampling-based simulations based on the thermal profiling case study.
