# import needed packages
import numpy as np
import pandas as pd
import shap
from digipd_ml.supervised.classification import NestedCV
from digipd_ml.plot_utils import plot_shap_values
from digipd_ml.plot_utils import plot_from_shap_values
from digipd_ml.utils import feature_importances_shap_values
from digipd_ml.utils import feature_importances_from_shap_values
from sklearn import metrics
import argparse

# parse cmd arguments
parser = argparse.ArgumentParser(prog='PD_control_finalmodel_plot.py',
                                 description='Nested CV at pathway level')
parser.add_argument('long_feat', type=str, help='Longitudinal feature')
parser.add_argument('model', type=str, help='Name of final model')
args = parser.parse_args()

# I/O
INPUT_DIR = "../data/"
OUTPUT_DIR = '../results/'
target = "DIAGNOSIS" # UPDRS__3_binary
model_name = args.model
temp_feat = args.long_feat
input_cv_file = temp_feat + "_data_cv_expr_4ML_" + target + ".tsv"
input_test_file = temp_feat + "_data_test_expr_4ML_" + target + ".tsv"
out_test_file = temp_feat + "_results_test_" + model_name + "_" + target + ".csv"
out_ft_file = temp_feat + "_top_features_shap_" + target + ".csv"
out_shap_plot = temp_feat + "_shap_plot_" + target + ".pdf"

# Reading the files
df = pd.read_table(INPUT_DIR + input_cv_file, index_col=0)
df_test = pd.read_table(INPUT_DIR + input_test_file, index_col=0)

# save diagnosis
y = df[target]
df = df.drop([target], axis=1)
y_test = df_test[target]
df_test = df_test.drop([target], axis=1)

# save features name and Patient_ID
features_name = [name for name in df.columns]
Patient_ID = list(df.index)

# Check that there is no more missing values
if np.sum(np.isnan(df)).sum() != 0:
    raise ValueError('There is missing data')

# transform data frame to numpy array
y = np.array(y)
X = np.array(df)
y_test = np.array(y_test)
X_test = np.array(df_test)

# fit single best model
nestedCV = NestedCV()
model = nestedCV.fit_single_model(
    X, y, name_model=model_name, normalize=True, feat_select=True, balanced=True)

# apply scaler and selecter before?
# predict & evaluate on test data
y_proba = model.predict_proba(X_test)[:, 1]
y_hat = model.predict(X_test)

auc = metrics.roc_auc_score(y_test, y_proba)
acc = metrics.accuracy_score(y_test, y_hat)
balanced_acc = metrics.balanced_accuracy_score(y_test, y_hat)
confusion = metrics.confusion_matrix(y_test, y_hat)
true_neg = confusion[0][0]
false_neg = confusion[1][0]
true_pos = confusion[1][1]
false_pos = confusion[0][1]
sensitivity = true_pos / (true_pos + false_neg)
specificity = true_neg / (true_neg + false_pos)

test_metrics = pd.DataFrame(data=[auc, acc, balanced_acc, sensitivity, specificity],
                            index=['auc', 'accuracy', 'balanced_accuracy', 'sensitivity', 'specificity'],
                            columns=["test_set"])
# export test_metrics
test_metrics.to_csv(OUTPUT_DIR + out_test_file)

# shap values analysis
X_processed = model['selecter'].transform(model['scaler'].transform(X))  # apply scaler and selecter
feat_names = model['selecter'].get_feature_names_out(features_name)  # features selected (index ordered)
X_processed = pd.DataFrame(data=X_processed, columns=feat_names)

# shap values analysis
if model_name not in ["RBFSVM", "Adaboost"]:
    print("Shaps from explainer")

    # shap plot
    plot_shap_values(model['model'], X_processed, save_fig=True,
                     plot_title="Gene level features",
                     name_file=out_shap_plot, path=OUTPUT_DIR)

    # feature importance
    top_features = feature_importances_shap_values(model['model'], X_processed, n=20)

else:
    print("Shaps from KernelExplainer")
    # generate shap values
    explainer = shap.KernelExplainer(model['model'].predict, shap.kmeans(X_processed, 100))
    shap_values = explainer.shap_values(X_processed)

    # shap plot
    plot_from_shap_values(shap_values, X_processed, save_fig=True,
                          plot_title="Gene level features",
                          name_file=out_shap_plot, path=OUTPUT_DIR)

    # feature importance
    top_features = feature_importances_from_shap_values(shap_values, X_processed, n=20)
top_features.to_csv(OUTPUT_DIR + out_ft_file, index=False)

print("end")



