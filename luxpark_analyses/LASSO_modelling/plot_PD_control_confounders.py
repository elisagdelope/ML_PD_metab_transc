import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# set seed
np.random.seed(111)

dir_results = "../results/"


def _plot_legend_apart(ax, figname, ncol=None):
    """Plot legend apart from figure."""
    # Do all your plots with fig, ax = plt.subplots(),
    # don't call plt.legend() at the end but this instead
    if ncol is None:
        ncol = len(ax.lines)
    fig = plt.figure(figsize=(30, 4), constrained_layout=True)
    fig.legend(ax.lines, [line.get_label() for line in ax.lines], ncol=ncol,
               loc="upper center")
    fig.tight_layout()
    fig.savefig(figname, bbox_inches="tight")
    os.system("pdfcrop %s %s" % (figname, figname))
    return fig


models = [
    "Linear SVM", "RBF SVM", "Random Forest", "Gradient Boosting",
    "Adaboost", "Logistic regression"
]

#dir_results = ""

# load metabolomics results
results_metab = pd.read_csv(dir_results+"results_nestedCV_DIAGNOSIS.csv", index_col=0)
results_metab = results_metab.T
# load clinical confounders only results
results_confounders = pd.read_csv(dir_results+"results_confounders_nestedCV_DIAGNOSIS.csv", index_col=0)
results_confounders = results_confounders.T

# load clinical + treatment confounders results
results_confounders_treatment = pd.read_csv(dir_results+"confounders_iris/results_confounders_treatment_nestedCV_DIAGNOSIS.csv", index_col=0)
results_confounders_treatment = results_confounders_treatment.T

# load metabolomics + clinical confounders results
results_all = pd.read_csv(dir_results+"results_withconf_nestedCV_DIAGNOSIS.csv", index_col=0)
results_all = results_all.T

# load metabolomics + clinical + treatment confounders results
results_all_treatment = pd.read_csv(dir_results+"confounders_iris/results_withconf_treatment_nestedCV_DIAGNOSIS.csv", index_col=0)
results_all_treatment = results_all_treatment.T

x_axis = np.arange(len(results_all.index))

#fig, ax = plt.subplots(figsize=(16, 8))
plt.figure(figsize=(14, 8))
fig, ax = plt.gcf(), plt.gca()
ax.errorbar(
    x_axis - 0.3, results_all_treatment['mean_auc'], yerr=results_all_treatment['std_auc'].values.flatten(),
    fmt='o', color='Black', elinewidth=2, capthick=2, errorevery=1, alpha=1,
    ms=4, capsize=5)
ax.errorbar(
    x_axis - 0.15, results_all['mean_auc'], yerr=results_all['std_auc'].values.flatten(),
    fmt='o', color='Black', elinewidth=2, capthick=2, errorevery=1, alpha=1,
    ms=4, capsize=5)
ax.errorbar(
    x_axis, results_metab['mean_auc'], yerr=results_metab['std_auc'].values.flatten(),
    fmt='o', color='Black', elinewidth=2, capthick=2, errorevery=1, alpha=1,
    ms=4, capsize=5)
ax.errorbar(
    x_axis + 0.15, results_confounders_treatment['mean_auc'],
    yerr=results_confounders_treatment['std_auc'].values.flatten(),
    fmt='o', color='Black', elinewidth=2, capthick=2, errorevery=1, alpha=1,
    ms=4, capsize=5)
ax.errorbar(
    x_axis + 0.3, results_confounders['mean_auc'], yerr=results_confounders['std_auc'].values.flatten(),
    fmt='o', color='Black', elinewidth=2, capthick=2, errorevery=1, alpha=1,
    ms=4, capsize=5)



ax.bar(
    x_axis - 0.3, results_all_treatment['mean_auc'], 0.15, tick_label=results_all_treatment.index,
    label="metabolomics features + clinical + treatment confounders", edgecolor='k', linewidth=2, color="tab:blue")
ax.bar(
    x_axis - 0.15, results_all['mean_auc'], 0.15, tick_label=results_all.index,
    label="metabolomics features + clinical confounders", edgecolor='k', linewidth=2, color="purple")
ax.bar(
    x_axis, results_metab['mean_auc'], 0.15, tick_label=results_metab.index,
    label="metabolomics features", edgecolor='k', linewidth=2, color="tab:green")
ax.bar(
   x_axis + 0.15, results_confounders_treatment['mean_auc'], 0.15, tick_label=results_confounders_treatment.index,
    label="clinical + treatment confounders", edgecolor='k', linewidth=2, color="tab:orange")
ax.bar(
   x_axis + 0.3, results_confounders['mean_auc'], 0.15, tick_label=results_confounders.index,
    label="clinical confounders", edgecolor='k', linewidth=2, color="tab:red")


ax.set_xticks(x_axis)
ax.set_xticklabels(models, fontsize=12)
ax.set_ylabel('Average 10-folds AUC', fontsize=12)       # Label on Y axis
ax.set_ylim((0, 1))
plt.legend(fontsize=12)
#fig.tight_layout()
fig.savefig(dir_results + "PD_vs_control_all_comparisons_treatment.pdf")
plt.show()
