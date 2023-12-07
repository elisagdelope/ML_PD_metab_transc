Scripts utilized for ML modelling of transcritpomics data from the PPMI cohort for PD case-control sample classification.


### Exploration 

##### ppmi_descriptive_plots.R
Generates descriptive figures of the dataset.

##### distributions.R
Tests for different data distributions.


### Preprocessing

##### ppmi_data4ML_class.R
Performs unsupervised filters to generate data for ML modelling of snapshot data (T0) from metabolomics data (all PD/HC and *de novo* PD/HC).

* ppmi_data4ML_class.R employs as input transcriptomics and phenotypical data resulting from previous pre-processing scripts described in repository *statistical_analyses_cross_long_PD* for **parsing data** and **Baseline (T0) PD/HC** (ppmi_filter_gene_expression.R, ppmi_norm_gene_expression.R, ppmi_generate_pathway_level.R). 



### ML modelling with DA feature selection

##### ppmi_nCV_DEA_MLclass_vf.R, 
Performs nested crossvalidation including feature selection (through DA) + binary ML classifiers on transcriptomics static data (BL). Includes unsupervised feature selection filters prior to the analyses.

##### ppmi_nCV_rf_MLclass.R
Performs nested crossvalidation in a Random Forest model for binary classification on transcriptomics static data (BL). Includes unsupervised feature selection filters prior to the analyses.


### ML modelling with LASSO feature selection at baseline and using temporal features

Please, install the environment package *digipd_analysis* to run these scripts.

##### PD_control_train.py, PD_control_train_PW.py
Performs nested crossvalidation for feature selection (through DA) + binary ML classifiers on transcriptomics static data (BL). Input taken from *ppmi_data4ML_class.R* outputs (which include the unsupervised feature selection filters).

##### PD_control_finalmodel_plot.py, PD_control_finalmodel_plot_PW.py
Performs PD/HC lassification with a selected model on the test set,  computes SHAP values and generates corresponding plots and tables.

##### DB_st_plots.ipynb, plots.py
Generate figures from the nested CV results on different dynamic features and baseline.



### Post-hoc analyses

##### ppmi_friedman.R
Performs comparisons among model performance, data types performance, and feature selection methods using the friedman test with Bergmann-Hommel corrections for pairwise comparisons, and with Holm corrections for one-to-many comparisons.



