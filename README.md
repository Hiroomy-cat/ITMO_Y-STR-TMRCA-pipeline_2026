Y-STR TMRCA Estimation Pipeline
This repository contains a set of tools and scripts developed for accurate Time to the Most Recent Common Ancestor (TMRCA) estimation based on Y-chromosomal Short Tandem Repeats (Y-STRs) from Whole Genome Sequencing (WGS) data.
Overview
The project focuses on validating classical mutation models (IAM, SMM, ASD) and implementing a Bayesian simulation approach (YMrCA) to reconstruct male lineages with high precision using extended 80+ marker panels.
File Descriptions
1. Data Preparation & Pre-processing
prep_for_matlab.py: A preprocessing script that cleans and formats WGS-derived STR genotypes (Platinum Pedigree dataset) for the MATLAB environment. It handles name normalization, converts allele formats, and generates per-locus mutation matrices.
podgotovka_filov.py: (Input Preparation) Automates the filtering of raw HipSTR/GangSTR outputs using locus-specific thresholds for depth of coverage (DP) and genotype quality (Q). It generates clean haplotype tables and identifies kinship pairs from metadata.
2. Core TMRCA Calculation
ymrca_platinum.m: The main Bayesian simulation script (MATLAB). Based on the YMrCA algorithm by Sofie Claerhout et al., modified for cluster performance. It includes an intermediate pruning mechanism to handle high-dimensional data (80+ loci) without memory overflow.
tmrca_from_diffs.py: Implementation of the Walsh (2001) probability grid method. It calculates TMRCA point estimates and confidence intervals for both Infinite Alleles (IAM) and Stepwise Mutation (SMM) models.
3. Simulation & Testing
simulate_ystr_family.py: A simulation tool used to generate synthetic reads with specific coverage levels (1x–20x) and realistic PCR stutter models. Used to benchmark HipSTR vs. STRSensor performance in variable-depth WGS (vdWGS) conditions.
4. Analysis & Visualization
analiz_grishi.py: A comprehensive analysis script that collects results from cluster simulations (coupleX.tsv files), calculates Bayesian consensus between models, and maps genetic TMRCA to autosomal IBD data.
plot_cousins.py: Visualization script focused on distant relatives. It generates correlation plots between whole-genome IBD sharing and Y-STR TMRCA, providing independent biological validation of the paternal molecular clock.
Requirements
Python 3.13+ (Pandas, NumPy, SciPy, Plotly, Seaborn)
MATLAB R2025b (for ymrca_platinum.m)
SLURM Workload Manager (for cluster execution via Job Arrays)
