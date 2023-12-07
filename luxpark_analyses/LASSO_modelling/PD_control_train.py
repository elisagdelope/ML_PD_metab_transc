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
#input_file = "DENOVO_data_metab_4ML_DIAGNOSIS.tsv"
#output_file = 'DENOVO_results_nestedCV_diagnosis.csv'

# Reading the files
df = pd.read_table(IN_DIR + input_file, index_col=0)

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
nestedCV = NestedCV()
nestedCV.fit(X, y, feat_select=True, normalize=True, balanced=True)

# save results into a csv
nestedCV.get_results().to_csv(OUTPUT_DIR + output_file)

