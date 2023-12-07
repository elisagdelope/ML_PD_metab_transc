# import needed packages
import numpy as np
import pandas as pd
from digipd_ml.supervised.classification import NestedCV

# set seed
np.random.seed(111)

# I/O
IN_DIR = "../data/"
OUTPUT_DIR = '../results/'
target = "DIAGNOSIS" # UPDRS__3_binary
input_file = "data_cv_metab_4ML_" + target + ".tsv"
output_file = "results_nestedCV_" + target + ".csv"
output_confounders_file = "results_confounders_nestedCV_" + target + ".csv"
output_confounders_treatment_file = "results_confounders_treatment_nestedCV_" + target + ".csv"
output_met_confounders_file = "results_withconf_nestedCV_" + target + ".csv"
output_met_confounders_treatment_file = "results_withconf_treatment_nestedCV_" + target + ".csv"

# Reading the files
df = pd.read_table(IN_DIR + input_file, index_col=0)
clinical = pd.read_csv(IN_DIR + "pheno_V0.tsv", sep="\t", index_col=2)
clinical = clinical[clinical.index.isin(df.index)]
treatment = pd.read_csv(IN_DIR + "M1342.tsv", sep="\t", index_col=0)
treatment.columns = ["MDOPA3"]
# Filter the indexes of treatment data based on the matching indexes of cv metab data
filtered_indexes = treatment.index.isin(df.index)
treatment = treatment[filtered_indexes]

# Check that patients from preprocessed data ID match diagnosis
if (np.sum(clinical.index == df.index) != df.shape[0]) or (np.sum(treatment.index == df.index) != df.shape[0]) :
    raise ValueError('Patient are not matching')

# save diagnosis
y = df[target]
df = df.drop([target], axis=1)

# save features name and Patient_ID
features_name = [name for name in df.columns]
Patient_ID = list(df.index)

# Check that there is no more missing values
if np.sum(np.isnan(df)).sum() != 0:
    raise ValueError('There is missing data')

# transform data frame to numpy array
y = np.array(y)
X = np.array(df)

# perform nestedCV using class implemented in nestedcv.py file
# the balanced option will perform down sampling if set to False
# meaning that classes are not balanced
#nestedCV = NestedCV()
#nestedCV.fit(X, y, feat_select=True, normalize=True, balanced=True)

# save results into a csv
#nestedCV.get_results().to_csv(OUTPUT_DIR + output_file)


# build confounders design matrix
from sklearn.preprocessing import OneHotEncoder
clinical = clinical.drop(["VISIT", "PATIENT_ID", "DIAGNOSIS"], axis = 1)
enc = OneHotEncoder(handle_unknown='ignore', sparse=False)
encoded_gender = enc.fit_transform(np.array(clinical['GENDER']).reshape(-1, 1))
clinical['MALE'] = encoded_gender[:, 0]
clinical['FEMALE'] = encoded_gender[:, 1]
clinical = clinical.drop('GENDER', axis=1)

from sklearn.impute import KNNImputer
imputer = KNNImputer(n_neighbors=10)
cols = list(clinical.columns)
ix = clinical.index
clinical = imputer.fit_transform(clinical)
clinical = pd.DataFrame(clinical, columns=cols, index=ix)

# create a clinical confounders and a joint clinical and treatment confounders df
confounders = clinical[['AGE', 'BMI', 'MALE']]
confounders_treatment = pd.merge(clinical[['AGE', 'BMI', 'MALE']], treatment, left_index=True, right_index=True)
confounders = np.array(confounders)
confounders_treatment = np.array(confounders_treatment)


# run baseline model with only clinical confounders
# Age, Gender, BMI
nestedCV = NestedCV()
nestedCV.fit(
    confounders, y, feat_select=False, normalize=True, balanced=True)
nestedCV.get_results().to_csv(OUTPUT_DIR + output_confounders_file)

# run baseline model with joint clinical and treatment confounders
# Age, Gender, BMI
nestedCV = NestedCV()
nestedCV.fit(
    confounders_treatment, y, feat_select=False, normalize=True, balanced=True)
nestedCV.get_results().to_csv(OUTPUT_DIR + output_confounders_treatment_file)

# run model using metabolomics features + clinical confounders
X = np.concatenate((X, confounders), axis=1)
features_name.extend(['AGE', 'BMI', 'MALE'])
nestedCV = NestedCV()
nestedCV.fit(X, y, feat_select=True, normalize=True, balanced=False)

# save results into a csv
nestedCV.get_results().to_csv(OUTPUT_DIR + output_met_confounders_file)

# run model using metabolomics features + joint clinical & treatment confounders
X = np.concatenate((X, confounders_treatment), axis=1)
features_name.extend(['AGE', 'BMI', 'MALE', 'MDOPA3'])
nestedCV = NestedCV()
nestedCV.fit(X, y, feat_select=True, normalize=True, balanced=False)

# save results into a csv
nestedCV.get_results().to_csv(OUTPUT_DIR + output_met_confounders_treatment_file)
