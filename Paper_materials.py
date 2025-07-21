#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
from sklearn.metrics import r2_score

def scatter_plot(x, y, title=None, xlab=None, ylab=None, x_log_scale=False, mid_line=False, alpha=0.5, colors=None):
    """
    Set a generic plot with the given parameters.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(x, y, alpha=alpha, c=colors)
    if title is not None: ax.set_title(title, fontsize=18)
    ax.set_xlabel(xlab, fontsize = 16)
    ax.set_ylabel(ylab, fontsize = 16)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.tick_params(axis='both', which='major', labelsize=14)
    if x_log_scale:
        ax.set_xscale('log')
    if mid_line:
        ax.plot([0, 1], [0, 1], color='red', linestyle='--', label='y=x')
        plt.legend()
    plt.show()      

def errorbar_plot(x, y, yerr, title=None, xlab=None, ylab=None, x_1 =None, y_1 = None, yerr_1=None, labels = ["A", "B"],x_log_scale=False):
    """
    Set a generic errorbar plot with the given parameters.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    if x_1 is not None and y_1 is not None and yerr_1 is not None:
        ax.errorbar(x, y, yerr=yerr, fmt='o', label=labels[0])
        ax.errorbar(x_1, y_1, yerr=yerr_1, fmt='o', label=labels[1])
        ax.legend(loc = "lower right")
    else:
        ax.errorbar(x, y, yerr=yerr, fmt='o', label=labels[0], color="black")
    if title is not None: ax.set_title(title, fontsize=18)
    ax.set_xlabel(xlab, fontsize=16)
    ax.set_ylabel(ylab, fontsize=16)
    ax.set_ylim(0, 1)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True, linestyle='--', alpha=0.5)
    if x_log_scale:
        ax.set_xscale('log')
    if any(len(str(x_))>10 for x_ in x):
        plt.xticks(rotation=45)
    plt.show()

#######################   Scatterplot for 5 species   #######################
def scatter_plot_5_sepcies(path):
    df = pd.read_csv(path)
    scatter_plot(df['Test_target'],
                 df['Test_predictions'],
                    xlab='Real proportion of AMR bacteria',
                    ylab='Predicted proportion of AMR bacteria',
                    x_log_scale=False,
                    mid_line=True, alpha=0.5)
    
#############################   Scatterplot multispecie   #######################
def scatter_plot_multispecie(path, xlab="Real proportion of AMR bacteria", ylab="Predicted proportion of AMR bacteria"):
    dati = pd.read_csv(path, index_col=0)
    scatter_plot(dati['Real'],
                 dati['Predictions'],
                 title=None,
                 xlab=xlab,
                 ylab=ylab,
                 x_log_scale=False,
                 mid_line=True, alpha=0.5)

########################   Depth effects andamento per depths   #######################
def R2_plot_from_depth(path_1, path_2=None, title = None):
    """
    Questa si usa sui file del depth effect con depth = 1
    """
    df = pd.read_csv(path_1)
    if path_2:
        df_2 = pd.read_csv(path_2)
        errorbar_plot(df['Subsample'],
                      df['Average R2'],
                      yerr=df['R2 errors'],
                      title=title,
                      xlab='Subsample Size (Gbp)',
                      ylab=r"$R^2$",
                      x_1=df_2['Subsample'],
                      y_1=df_2['Average R2'],
                      yerr_1=df_2['R2 errors'],
                      labels = ["original strain", "5 sub-strain"],
                      x_log_scale=False)
    else:
        errorbar_plot(df['Subsample'],
                      df['Average R2'],
                      yerr=df['R2 errors'],
                      title=title,
                      xlab='Subsample Size (Gbp)',
                      ylab=r"$R^2$",
                      x_log_scale=False)

##########################   Complexity effect andamento per complexity   #######################
def R2_plot_from_complexity(path_1, path_2=None, title = None):
    """
    Questa si usa sui file del complexity effect con complexity = 1
    """
    with open(path_1, 'r') as file:
        data = json.load(file)
    df = {model: reuslts["Test_scores"] for model, reuslts in data.items()}
    df = pd.DataFrame.from_dict(df).T
    if path_2:
        with open(path_2, 'r') as file:
            data2 = json.load(file)
        df2 = {model: reuslts["Test_scores"] for model, reuslts in data2.items()}
        df2 = pd.DataFrame.from_dict(df2).T
    #Plot a scatter plot with error bars
    errorbar_plot(df.columns,
                    df.mean(),
                    yerr=df.std(),
                    title=title,
                    ylab=r"$R^2$",
                    x_1=df2.columns if path_2 else None,
                    y_1=df2.mean() if path_2 else None,
                    yerr_1=df2.std() if path_2 else None,
                    labels=["original strain", "5 sub-strains"] if path_2 else ["R2 score"],
                    x_log_scale=False)

############################   Scatterplot complexity   #######################
def scatter_plot_from_complexity(path, replica_N = 0):
    """
    Questa si usa sui file del complexity effect con complexity = 1
    """
    with open(path, 'r') as file:
        data = json.load(file)
    data = data[list(data.keys())[replica_N]]["Test_predictions"]
    R_tot = None
    for k, v in data.items():
        R  = pd.DataFrame(v).T
        if R_tot is None:
            R_tot = R
        else:
            R_tot = pd.concat([R_tot, R], axis=0)
        t = None
        scatter_plot(R['Test_target'], R['Test_predictions'],
                     title=None,
                     xlab='Real proportion of AMR bacteria',
                     ylab='Predicted proportion of AMR bacteria',
                     x_log_scale=False,
                     mid_line=True, alpha=0.5, colors=R['Color'])
    scatter_plot(R_tot['Test_target'], R_tot['Test_predictions'],
                 title=t,
                 xlab='Real proportion of AMR bacteria',
                 ylab='Predicted proportion of AMR bacteria',
                 x_log_scale=False,
                 mid_line=True, alpha=0.5)

##############################   Feature importance plot   #######################
def feature_importance_plot(path, how_many_cols=30, x_axes_upper_limit=None):
    Features = pd.read_csv(path, index_col=0).T
    Sorted_col = list(Features.mean(axis=0).sort_values(ascending=False).index)
    Just_the_col_i_need = Sorted_col[:how_many_cols]
    Just_the_col_i_need = Just_the_col_i_need[::-1]
    sorted_df = Features[Just_the_col_i_need]
    for name in sorted_df.columns:
        if "ARO" in name.split("|")[0]:
            new_name = "ARO_ID:" + name.split("|")[0].split(":")[1]
        else:
            new_name = "ARO_ID:" + name.split("|")[3].split(":")[1] + "|Wildcard_seq:" + name.split("|")[0].split(":")[1]
        sorted_df.rename(columns={name: new_name}, inplace=True)
    fig, ax = plt.subplots(figsize=(25, 12))
    box = sorted_df.boxplot(ax=ax, vert=False, patch_artist=True)
    if x_axes_upper_limit is not None:
        ax.set_xlim(0, x_axes_upper_limit)
    #Colora le box tutte dello stesso colore
    colors = ["#1f77b4"] * len(sorted_df.columns)  # o qualsiasi altro colore desiderato
    for patch, color in zip(ax.patches, colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.set_xlabel('Gini importance', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=18)
    plt.tight_layout()
    return sorted_df

#############################   Depth multispecie   #######################
def depth_multispecie(path):
    dati = pd.read_csv(path, index_col=0)
    x_ax = [float(x.split("Subsample ")[1]) for x in dati.index]
    if 1 in x_ax:
        x_ax = [x / 10 for x in x_ax]  # Convert to Gbp if needed
    errorbar_plot(x_ax,
                  dati['R2'],
                  yerr=dati['R2_std'],
                  title=None,
                  xlab='Subsample Size (Gbp)',
                  ylab=r"$R^2$",
                  #Punti di colore nero
                  labels=["Multispecie"])

#%%
R2_plot_from_depth("New/Depth_1_strains_10cv/R2s_depth_effects.csv",
                   "New/Depth_5_strains_10cv/R2s_depth_effects.csv")
#%%
#scatter_plot_from_depth("New/Depth_1_strains/Depth_subsample_0.2.csv")
scatter_plot_5_sepcies("New/Depth_1_strains_10cv/Depth_subsample_0.2.csv")
#scatter_plot_from_path("New/Depth_5_strains/Depth_subsample_1.csv")
#%%
R2_plot_from_complexity("New/Complexity_sub0.02_1_strains/Results.json",
                        "New/Complexity_sub0.02_5_strains/Results.json")
#%%
scatter_plot_from_complexity("New/Complexity_sub0.02_1_strains/Results.json", replica_N=1)
scatter_plot_from_complexity("New/Complexity_sub0.02_5_strains/Results.json", replica_N=2)
#%%
scatter_plot_multispecie("New/Depth_multispecie_prev_filter_0.07/scatterplot_data_Pathogens_depth_0.3.csv", xlab="Real proportion of Pathogens", ylab="Predicted proportion of Pathogens")
scatter_plot_multispecie("New/Depth_multispecie_prev_filter_0.07/scatterplot_data_Resistant_depth_0.8.csv")
#%%
depth_multispecie("New/Depth_multispecie_prev_filter_0.07/R2_depth_effects_Pathogens.csv")
depth_multispecie("New/Depth_multispecie_prev_filter_0.07/R2_depth_effects_Resistant.csv")
#%%
feature_importance_plot("New/Feature_imp/Single_specie_PATRIC_5_strains/Feature_importances.csv", x_axes_upper_limit=0.5, how_many_cols=10)
# %%
feature_importance_plot("New/Feature_imp/Multispecie_resistant/Feature_importances.csv", x_axes_upper_limit=0.5, how_many_cols=10)
feature_importance_plot("New/Feature_imp/Multispecie_pathogen/Feature_importances.csv", x_axes_upper_limit=0.5, how_many_cols=10)

# %%
