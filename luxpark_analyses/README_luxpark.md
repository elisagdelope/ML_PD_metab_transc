Scripts utilized for ML modelling of metabolomics data from the luxPARK cohort for PD case-control sample classification.


### Exploration 

##### lx_descriptive_plots.R, lx_descriptive_plots_V0.R
Generate descriptive figures of the entire dataset and baseline data (V0) respectively.

##### lx_imputation.R
Performs multiple imputation on batch-normalized metabolomics data, to compare with Metabolon imputation method.

##### lx_plot_distributions.R
Generate exploratory plots on the distribution of the data


### Preprocessing

##### lx_data4ML_class.R, lx_data4ML_class_denovo.R
Performs unsupervised filters to generate data for ML modelling of snapshot data (T0) from metabolomics data (all PD/HC and *de novo* PD/HC).

* lx_data4ML_class.R and lx_data4ML_class_denovo.R use as input metabolomics and phenotypical data resulting from previous pre-processing scripts described in repository [statistical_analyses_cross_long_PD](https://gitlab.lcsb.uni.lu/elisa.gomezdelope/statistical_analyses_cross_long_pd) for **parsing data** and **Baseline (T0) PD/HC** (lx_extract_visit.R, lx_denovo_filter.R, lx_generate_pathway_level.R). 


### ML modelling with DA feature selection

##### lx_nCV_DEA_MLclass.R, 
Performs nested crossvalidation including feature selection (through DA) + binary ML classifiers on metabolomics static data (V0). Includes unsupervised feature selection filters prior to the analyses.

##### lx_nCV_rf_MLclass.R
Performs nested crossvalidation in a Random Forest model for binary classification on metabolomics static data (V0). Input taken from *lx_data4ML_class.R*, *lx_data4ML_class_denovo.R* outputs (which include the unsupervised feature selection filters).

##### lx_shapley.R
Calculates shapley values from a given final model and its training set.


### ML modelling with LASSO feature selection

Please, install the environment package *digipd_ml* with the class **NestedCV()** to run these scripts.

##### PD_control_train.py, PD_control_train_PW.py
Performs nested crossvalidation for feature selection (through DA) + binary ML classifiers on metabolomics static data (V0). Input taken from *lx_data4ML_class.R*, *lx_data4ML_class_denovo.R* outputs (which include the unsupervised feature selection filters).

##### PD_control_train_confounders.py
Performs nested crossvalidation for the analysis of clinical and treatment-related confounding factors. 

##### PD_control_finalmodel_plot.py, PD_control_finalmodel_plot_PW.py
Performs PD/HC lassification with a selected model on the test set,  computes SHAP values and generates corresponding plots and tables.

##### PD_control_finalmodel_plot_notest.py
Utilizes a selected model to compute SHAP values and generate corresponding plots and tables in *de novo* PD vs. control classification.

##### plot_PD_control_confounders.py
Genrates figure for the confounder analysis.



### Post-hoc analyses

##### lx_friedman.R
Performs comparisons among model performance, data types performance, and feature selection methods using the friedman test with Bergmann-Hommel corrections for pairwise comparisons, and with Holm corrections for one-to-many comparisons.



