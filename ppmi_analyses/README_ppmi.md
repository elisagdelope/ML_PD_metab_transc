Scripts utilized for ML modelling of transcritpomics data from the PPMI cohort for PD case-control sample classification.


### Exploration 

##### ppmi_descriptive_plots.R
Generate descriptive figures of the dataset.

##### distributions.R
Test for different data distributions.


### Preprocessing

##### ppmi_data4ML_class.R
Perform unsupervised filters to generate data for ML modelling of snapshot data (T0) from metabolomics data (all PD/HC and *de novo* PD/HC).

* ppmi_data4ML_class.R employs as input transcriptomics and phenotypical data resulting from previous pre-processing scripts described in repository *statistical_analyses_cross_long_PD* for **parsing data** and **Baseline (T0) PD/HC** (ppmi_filter_gene_expression.R, ppmi_norm_gene_expression.R, ppmi_generate_pathway_level.R). 



### ML modelling with DA feature selection

##### ppmi_nCV_DEA_MLclass_vf.R, 
Perform nested crossvalidation including feature selection (through DA) + binary ML classifiers on transcriptomics static data (BL). Includes unsupervised feature selection filters prior to the analyses.

##### ppmi_nCV_rf_MLclass.R
Perform nested crossvalidation in a Random Forest model for binary classification on transcriptomics static data (BL). Includes unsupervised feature selection filters prior to the analyses.


### ML modelling with LASSO feature selection at baseline and using temporal features

Please, install the environment package *digipd_ml* with the class **NestedCV()** to run these scripts.

##### PD_control_train.py, PD_control_train_PW.py
Perform training for PD/HC classification with multiple binary classifiers on transcriptomics data in a nested cross-validation setting, allowing for undersampling, feature scaling, and feature selection with LASSO. Input taken from *ppmi_data4ML_class.R*, *ppmi_data4ML_TS_class.R* outputs (which include the unsupervised feature selection filters).


##### PD_control_finalmodel_plot.py, PD_control_finalmodel_plot_PW.py
Train a selected model on the entire training set for PD/HC classification and validate on the test set, compute SHAP values and generate corresponding plots and tables.

##### DB_st_plots.ipynb, plots.py
Generate figures from the nested CV results on different dynamic features and baseline.



### Post-hoc analyses

##### ppmi_friedman.R
Perform comparisons among model performance, data types performance, and feature selection methods using the friedman test with Bergmann-Hommel corrections for pairwise comparisons, and with Holm corrections for one-to-many comparisons.



