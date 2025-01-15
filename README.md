# MSstatsWeightedSummary: reproduction code

This repository contains code required to reproduce the evaluation of weighted protein summarization method described in an article "Relative quantification of proteins and post-translational modifications in proteomic experiments with shared peptides: a weight-based approach".
Due to GitHub's file size limits, input data necessary for running the code are available in a separate location: [input\_data folder required for full reproduction](https://drive.google.com/drive/folders/1iJ5M3cEiWrEvgI3SJbLmO7kwSYAOr2w4?usp=sharing).
However, all necessary intermediate results are available in the processed\_data folder, so the raw inputs are not strictly necessary to run the data analyses.

After obtaining the folder with all input files, it is possible to run the following scripts:

- 00_setup.R to install necessary packages, including two currently GitHub-based libraries: MSstatsWeightedSummary (which implements the weighted summarization method) and SimulateTMT (which implements model-based protein- and peptide-level data simulation);
- 01_graphs_explanation.R to reproduce peptide-protein graphs and profiles plots of the protein degrader case study;
- 02_brd_simulations_main.R to reproduce the simulation study based on protein degrader data which evaluates the impact of proposed approach on log-fold change estimation;
- 03_ptm_case_study.R to reproduce the summarization of an example cluster of PTM sites;
- 04_simulated_data.R to reproduce the model-based simulation study which evaluates the impact of proposed approach on log-fold change estimation and statistical inference results; 
- 05_tpp_onepot_data_preparation.R to reproduce the re-processing of protein clusters from onePot and TPP data;
- 06_tpp_onepot_analysis.R to reproduce the re-analysis of onePot and TPP data;
- 07_analytical_cov_matrix_simulated_data.R to reproduce a model-based simulation which evaluates the variance of weights estimation in a simple scenario;
- 08_brd_simulations_si.R to reproduce additional protein degrader-study based simulation describing the influence of loss functions, starting points and lack of unique peptides on the quality of log-fold change estimation.
