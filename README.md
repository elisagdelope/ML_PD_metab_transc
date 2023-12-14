# Table of contents
* [Introduction](#introduction)
* [Content](#content)
* [Data](#data)
* [Requirements](#requirements)
* [License](#license)

## Introduction
This repository contains the code for ML analyses performed in Chapter 4 of my PhD thesis "Interpretable Machine Learning on omics data for biomarker discovery in Parkinson's disease". The project consists on performing Parkinson's disease (PD) case-control classification from blood plasma metabolomics measurements at the baseline clinical visit from the LuxPARK cohort, and from whole blood transcriptomics data at baseline as well as dynamic features engineered from a short temporal series of 4 timepoints from the PPMI cohort. The study involves evaluation of different feature selection strategies,  The goal was to build and test a collection of ML models and, most interestingly, identify molecular and higher-level functional representations associated with PD diagnosis.

## Content
The code covers the following main tasks and analyses:

1. Exploration of the data 
2. Preprocessing of the datasets, generation of temporal (dynamic) features
3. Feature selection strategies for classification of PD vs. controls, including unsupervised filters, statistical differential analysis filter, and LASSO regression.
4. Analysis of confounding factors
5. Classification of PD vs. controls using baseline measurements and extracted dynamic features in a nested cross-validation setting, and evaluation on test set for molecular and higher-level functional representations as predictors.
6. Feature importance with SHAP values.
7. Post-hoc analysis with Friedman tests to compare performance results, and plots on the performance.

There is a README.md inside each directory with corresponding explanation for each script:

*ppmi_analyses* contains code related to analysis on transcriptomics data from PPMI.

*luxpark_analyses* contains code related to analysis on metabolomics data from LuxPARK.


## Data
The public transcriptomics data used in this project was derived from the Parkinson’s Progression Markers Initiative (https://www.ppmi-info.org/, RNAseq - IR3).
The metabolomics data from LuxPARK is not publicly available as it is linked to the Luxembourg Parkinson’s Study and its internal regulations. Any requests for accessing the dataset can be directed to request.ncer-pd@uni.lu.

## Requirements
The code for ML modeling and plots was implemented in Python (3.9.13), and that for the pre-processing steps and post-hoc analysis, was implemented in R (R 4.0.3). It has been tested on both current Mac (Ventura) and Linux operating systems (Rocky Linux 8.7 (Green Obsidian)), relying on multiple Python, R, and BioConductor packages correspondingly that are listed at the beginning of each script. The code should be compatible with later versions of Python and R installed on current Mac, Linux or Windows systems. **For the ML modeling, it is necessary to download and install the NestedCV() class implemented in the digipd_ml package and environment from repository [digipd_ml](https://gitlab.lcsb.uni.lu/elisa.gomezdelope/digipd_ml).**

## License
TBD


