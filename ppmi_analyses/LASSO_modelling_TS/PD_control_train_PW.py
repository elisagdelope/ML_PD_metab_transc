# import needed packages
import argparse
import numpy as np
import pandas as pd
from digipd_ml.supervised.classification import NestedCV

# parse cmd arguments
parser = argparse.ArgumentParser(prog = 'PD_control_train_PW.py',
                                 description = 'Nested CV at pathway level')
parser.add_argument('db', type=str, help='Aggregation database')
parser.add_argument('st', type=str, help='Aggregation statistic')
parser.add_argument('tempft', type=str, help='Temporal feature to use')
args = parser.parse_args()

# I/O
IN_DIR = "../data/"
OUTPUT_DIR = '../results/'
temp_feature = args.tempft
db = args.db
st = args.st
target = "DIAGNOSIS"
input_file = db + "_" + st + "_" + temp_feature + "_data_cv_expr_4ML_" + target + ".tsv"
output_file = db + "_" + st + "_" + temp_feature + "_results_nestedCV_" + target + ".csv"

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
nestedCV.fit(X, y, feat_select=True, normalize=True, balanced=False)

# save results into a csv
nestedCV.get_results().to_csv(OUTPUT_DIR + output_file)



sbatch pylauncher.sbatch PD_control_train_PW.py GOBP pca sd
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP pca lm-time
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP pca lm-lag
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP mean sd
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP mean lm-time
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP mean lm-lag
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP median sd
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP median lm-time
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP median lm-lag
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP pathifier sd
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP pathifier lm-time
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP pathifier lm-lag
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP sd sd
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP sd lm-time
sbatch pylauncher.sbatch PD_control_train_PW.py GOBP sd lm-lag

sbatch pylauncher.sbatch PD_control_train_PW.py GOCC pca sd
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC pca lm-time
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC pca lm-lag
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC mean sd
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC mean lm-time
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC mean lm-lag
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC median sd
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC median lm-time
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC median lm-lag
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC pathifier sd
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC pathifier lm-time
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC pathifier lm-lag
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC sd sd
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC sd lm-time
sbatch pylauncher.sbatch PD_control_train_PW.py GOCC sd lm-lag
