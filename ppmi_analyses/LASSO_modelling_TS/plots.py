import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def collect_cvperf_files(directory, cv_pattern):
    combined_df = pd.DataFrame()

    for filename in os.listdir(directory): 
        if filename.endswith(cv_pattern): 
            file_path = os.path.join(directory, filename)
            # reading content into data frame
            df = pd.read_csv(file_path, index_col=0)
            if not df.empty:
                combined_df = combined_df.append(df.iloc[0], ignore_index=True)
            
    combined_df.index = [filename.replace("_"+cv_pattern, "") for filename in os.listdir(directory) if filename.endswith(cv_pattern)]
    combined_df = combined_df.rename(index={cv_pattern: "GENE"})

    # Specify the desired order of the index
    custom_order = ['GENE', 'GOBP_mean', 'GOBP_median', 'GOBP_pathifier', 'GOBP_pca', 'GOBP_sd', 'GOCC_mean', 'GOCC_median', 'GOCC_pathifier', 'GOCC_pca', 'GOCC_sd', 'CORUM_mean', 'CORUM_median', 'CORUM_pathifier', 'CORUM_pca', 'CORUM_sd']
    combined_df = combined_df.reindex(custom_order)
    return combined_df

def multi_heatmap(list_df, titles, savefig=True):
    # Create subplots with the desired layout
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(11, 6), sharey=True, sharex=True)
    # Initialize minimum and maximum values
    vmin = np.min(np.array(list_df))
    vmax = np.max(np.array(list_df))
    # Loop through the dataframes and corresponding axes to plot the heatmaps
    for i, (df, ax, title) in enumerate(zip(list_df, axes.flatten(), titles)):
        transposed_df = df.transpose()
        sns.heatmap(data=transposed_df, ax=ax, xticklabels=True, yticklabels=True, cbar=False, cmap="viridis_r", annot=False, vmin=vmin, vmax=vmax)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha='right', rotation_mode='anchor', fontsize=9)  # Rotate x-axis labels by 35 degrees
        ax.tick_params(axis='y', labelsize=9)
    # Adjust the position of the subplots
    fig.subplots_adjust(right=0.98, bottom=0.2, wspace=0.05, hspace=0.05)
    # Add a colorbar to the figure
    cbar_ax = fig.add_axes([0.35, 0.05, 0.3, 0.02]) # [left, bottom, width, height]
 #   cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.02])
    cbar = fig.colorbar(axes[0, 0].collections[0], cax=cbar_ax, orientation='horizontal')
 #   cbar.set_clim(vmin, vmax)  # Set colorbar limits based on overall values
    cbar.mappable.set_clim(vmin=vmin,vmax=vmax)
    cbar.set_ticks([0.45, 0.50, 0.55, 0.60, 0.65]) # np.linspace(cbar.vmin, cbar.vmax, 4)
    cbar.set_ticklabels(np.round(cbar.get_ticks(), decimals=2))
    cbar.ax.tick_params(labelsize=9)
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_title('AUC')
    # Adjust the layout
    fig.tight_layout()
    # Set x-label and y-label for the overall figure
    fig.text(0.05, 0.15, 'Data types', ha='center', va='center', weight='bold')
    fig.text(0.01, 0.87, 'Models', ha='center', va='center', rotation='vertical', weight='bold')
    # Display the plot
    plt.show()
    # Save plot
    if savefig:
        fig.savefig('Heatmap_AUC_models_BLTS.png')
        fig.savefig('Heatmap_AUC_models_BLTS.pdf')



def multi_heatmap_vertical(list_df, titles, savefig=True):
    # Create subplots with the desired layout
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(7, 5), sharey=True, sharex=True)
    # Initialize minimum and maximum values
    vmin = np.min(np.array(list_df))
    vmax = np.max(np.array(list_df))
    # Loop through the dataframes and corresponding axes to plot the heatmaps
    for i, (df, ax, title) in enumerate(zip(list_df, axes.flatten(), titles)):
        sns.heatmap(data=df, ax=ax, xticklabels=True, yticklabels=True, cbar=False, cmap="viridis_r", annot=False, vmin=vmin, vmax=vmax)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha='right', rotation_mode='anchor', fontsize=9)  # Rotate x-axis labels by 35 degrees
        if i == 0:
            ax.set_yticks(np.arange(len(ax.get_yticklabels())) + 0.5)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right', fontsize=9)

    # Adjust the position of the subplots
    fig.subplots_adjust(right=0.95, bottom=0.22, wspace=-0.5)
    #fig.subplots_adjust(right=0.92, bottom=0.2, wspace=0.02)
    # Display the colorbar on the right
    cbar_ax = fig.add_axes([0.88, 0.35, 0.02, 0.35])  # Adjust the position and size of the colorbar
    cbar = fig.colorbar(axes[2].collections[0], cax=cbar_ax, orientation='vertical')
    cbar.mappable.set_clim(vmin=vmin, vmax=vmax)
    cbar.set_ticks([0.45, 0.50, 0.55, 0.60, 0.65])
    cbar.set_ticklabels(np.round(cbar.get_ticks(), decimals=2))
    cbar.ax.tick_params(labelsize=9)
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_title('AUC')
    # Set x-label and y-label for the overall figure
    fig.text(0.05, 0.18, 'Models', ha='center', va='center', weight='bold')
    fig.text(0.02, 0.8, 'Data types', ha='center', va='center', rotation='vertical', weight='bold')
    # Save plot
    if savefig:
        fig.savefig('Heatmap_AUC_models_vs_datatypes.png')
        fig.savefig('Heatmap_AUC_models_vs_datatypes.pdf')
    # Display the plot
    plt.show()


list_df = []
#prefix_pattern = ["", "lm-time_", "lm-lag_", "sd_"]
prefix_pattern = ["", "lm-time_", "sd_"]
for cv_pattern in prefix_pattern:
    if cv_pattern == "":
        IN_DIR = "BL-UPDRS3/results" #"BL-PD/results"
    else:
        IN_DIR = "TS-UPDRS3/results" #"TS-PD/results"
    cv_pattern = cv_pattern+"results_nestedCV_UPDRS3_binary.csv"
    df = collect_cvperf_files(IN_DIR, cv_pattern)
    df.columns = ["SVM Linear", "SVM Radial", "Random Forest", "Gradient Boosting", "Adaboost", "Logistic Regression"]
    df = df[["Random Forest", "Logistic Regression", "SVM Linear", "SVM Radial", "Gradient Boosting", "Adaboost"]]
    list_df.append(df)

# Generate figure with 4 heatmaps
#titles = ['Baseline (T0)', 'LM (expr/time)', 'LM (expr/lag1(expr))', 'SD (expr)']  # Replace with desired titles
titles = ['Baseline (T0)', 'LM (expr/time)', 'SD (expr)'] 
multi_heatmap_vertical(list_df, titles, savefig=True)



