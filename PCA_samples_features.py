import os
import pandas as pd
import M_learning_images as im
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.cm as cm

JUST_COMPLETE = True
NUMBER_OF_STRAINS = 1
MAPQ_MIN = 0
MAPQ_MAX = 1000
SCORING = "r2" # "r2" or "mean_squared_error"

campioni = pd.read_csv("Database_campioni2.csv", index_col=0)
t_train_PATRIC = campioni[campioni["Size"] == 0.1]
t_train_PATRIC = t_train_PATRIC[t_train_PATRIC["Real/Simul"] == "Simul"]
t_train_PATRIC = t_train_PATRIC[t_train_PATRIC["Aureus_S_N_strains"] == NUMBER_OF_STRAINS]
average_read_length_train = 270
assert all([x==average_read_length_train for x in t_train_PATRIC["reads_length"]])
Train_paths = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_train_PATRIC.index]
Train_data_PATRIC = im.read_many_alignements_json(Train_paths, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, mapq_average_max=MAPQ_MAX, read_length=average_read_length_train)
Train_target_PATRIC = dict((t_train_PATRIC["Aureus_R"] + t_train_PATRIC["Coli"] + t_train_PATRIC["Pseudomonas"]) /
				  (t_train_PATRIC["Arcobacter"] + t_train_PATRIC["Aureus_S"] + t_train_PATRIC["Aureus_R"] + t_train_PATRIC["Coli"] + t_train_PATRIC["Pseudomonas"]))
Train_target_PATRIC_complexity = t_train_PATRIC[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]] != 0
Train_target_PATRIC_complexity = Train_target_PATRIC_complexity.astype(int).sum(axis=1)
print("Ho letto le features")

campioni = pd.read_csv("Database_campioni_assembled.csv", index_col=0)
t_train_Assembled = campioni[campioni["Size"] == 0.1]
t_train_Assembled = t_train_Assembled[t_train_Assembled["Real/Simul"] == "Simul"]
t_train_Assembled = t_train_Assembled[t_train_Assembled["N_strains_for_all"] == NUMBER_OF_STRAINS]
average_read_length_train = 255
assert all([x==average_read_length_train for x in t_train_Assembled["reads_length"]])
Train_paths = ["samples_4/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_train_Assembled.index]
Train_data_Assembled = im.read_many_alignements_json(Train_paths, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, mapq_average_max=MAPQ_MAX, read_length=average_read_length_train)
Train_target_Assembled = dict((t_train_Assembled["Aureus_R"] + t_train_Assembled["Coli"] + t_train_Assembled["Pseudomonas"] +
				   t_train_Assembled["Aureus_R_1"] + t_train_Assembled["Coli_1"] + t_train_Assembled["Pseudomonas_1"]) / 
				   (t_train_Assembled["Arcobacter"] + t_train_Assembled["Aureus_S"] + t_train_Assembled["Aureus_R"] + t_train_Assembled["Coli"] + t_train_Assembled["Pseudomonas"] + 
	 				t_train_Assembled["Arcobacter_1"] + t_train_Assembled["Aureus_S_1"] + t_train_Assembled["Aureus_R_1"] + t_train_Assembled["Coli_1"] + t_train_Assembled["Pseudomonas_1"]))
print("Ho fatto il trian set, features e target")
columns_to_combine = [
    ("Aureus_S", "Aureus_S_1"),
    ("Aureus_R", "Aureus_R_1"),
    ("Coli", "Coli_1"),
    ("Pseudomonas", "Pseudomonas_1"),
    ("Arcobacter", "Arcobacter_1")
]
Train_Assembled_complexity = pd.DataFrame()
t_train_Assembled_ = t_train_Assembled[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter", "Aureus_S_1", "Aureus_R_1", "Coli_1", "Pseudomonas_1", "Arcobacter_1"]] != 0
for col1, col2 in columns_to_combine:
    new_col_name = f"{col1}_combined"
    Train_Assembled_complexity[new_col_name] = (t_train_Assembled_[col1] | t_train_Assembled_[col2]).astype(int)
Train_Assembled_complexity = Train_Assembled_complexity.astype(int).sum(axis=1)

campioni_per_test = pd.read_csv("Database_campioni2.csv", index_col=0)
t_test = campioni_per_test[campioni_per_test["Size"] == 0.1]
t_test = t_test[t_test["Real/Simul"] == "Real"]
average_read_length_test = 250
assert all([x==average_read_length_test for x in t_test["reads_length"]])
t_test = t_test[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
Test_paths = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_test.index]
Test_data = im.read_many_alignements_json(Test_paths, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, read_length=average_read_length_test)
Test_target = dict(t_test[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1) / t_test[["Arcobacter", "Aureus_S", "Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1))	
select_target_complexity = t_test != 0
select_target_complexity = select_target_complexity.astype(int).sum(axis=1)
print("Ho letto il test set")
#seleziona 1 colore per ogni valore dierso di complessitàù
test_complexity_colors = []
for x in select_target_complexity:
    if x == 1:
        test_complexity_colors.append("#000000")
    elif x == 2:
        test_complexity_colors.append("#440154")
    elif x == 3:
        test_complexity_colors.append("#3b528b")
    elif x == 4:
        test_complexity_colors.append("#21918c")
    elif x == 5:
        test_complexity_colors.append("#5ec962")

test_target_values = list(Test_target.values())
norm = plt.Normalize(min(test_target_values), max(test_target_values))
colormap = cm.viridis  # Puoi scegliere altre mappe di colori come 'plasma', 'cividis', ecc.
test_target_colors = [colormap(norm(value)) for value in test_target_values]

Train_data_PATRIC = Train_data_PATRIC.T
Train_data_PATRIC_5_genomes = Train_data_PATRIC[(Train_target_PATRIC_complexity == 5).values]
Train_data_Assembled = Train_data_Assembled.T
Train_Assembled_5_genomes = Train_data_Assembled[(Train_Assembled_complexity == 5).values]
Test_data = Test_data.T
Test_data_complexity_5 = Test_data[select_target_complexity == 5]

Tot_cols = set(Train_data_PATRIC.columns.to_list() + Train_data_Assembled.columns.to_list() + Test_data.columns.to_list())
print(len(Tot_cols))
#Put all the dataframe together adding 0 for the missing columns
Train_data_PATRIC = Train_data_PATRIC.reindex(columns=Tot_cols, fill_value=0)
Train_data_Assembled = Train_data_Assembled.reindex(columns=Tot_cols, fill_value=0)
Test_data = Test_data.reindex(columns=Tot_cols, fill_value=0)

#%%
#Append the dataframes
Tot_data = pd.concat([Train_data_PATRIC, Train_data_Assembled, Test_data], axis=0)
Data_Train_sets = pd.concat([Train_data_PATRIC, Train_data_Assembled], axis=0)
#make a PCA coloring it based on the origin of the data
Data_for_PCA_fitting = Tot_data
pca = PCA(n_components=2)
#Standardizza i dati
PCA_scaler = StandardScaler()
PCA_scaler = PCA_scaler.fit(Data_for_PCA_fitting)
print("Ho fittato lo scaler")
Data_for_PCA_fitting = PCA_scaler.transform(Data_for_PCA_fitting)
pca.fit(Data_for_PCA_fitting)
print("Ho fittato la PCA")
Train_data_PATRIC_PCA = pca.transform(PCA_scaler.transform(Train_data_PATRIC))
print("Ho fatto la PCA di PATRIC")
Train_data_Assembled_PCA = pca.transform(PCA_scaler.transform(Train_data_Assembled))
print("Ho fatto la PCA di Assembled")
Test_data_PCA = pca.transform(PCA_scaler.transform(Test_data))
#print("Ho fatto la PCA di Test")
#%%
fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # Crea una figura con due subplot affiancati

# Primo scatterplot
axes[0].scatter(Train_data_PATRIC_PCA[:, 0], Train_data_PATRIC_PCA[:, 1], c="red", alpha=0.3, label="PATRIC")
axes[0].scatter(Train_data_Assembled_PCA[:, 0], Train_data_Assembled_PCA[:, 1], c="orange", alpha=0.3, label="Assembled")
axes[0].scatter(Test_data_PCA[:, 0], Test_data_PCA[:, 1], c=test_complexity_colors, label="Test (Complexity)")
axes[0].set_title("Scatterplot con test_complexity_colors")
axes[0].set_xlabel("PC1")
axes[0].set_ylabel("PC2")

# Secondo scatterplot
axes[1].scatter(Train_data_PATRIC_PCA[:, 0], Train_data_PATRIC_PCA[:, 1], c="red", alpha=0.3, label="PATRIC")
axes[1].scatter(Train_data_Assembled_PCA[:, 0], Train_data_Assembled_PCA[:, 1], c="orange", alpha=0.3, label="Assembled")
axes[1].scatter(Test_data_PCA[:, 0], Test_data_PCA[:, 1], c=test_target_colors, label="Test (Target)")
axes[1].set_title("Scatterplot con test_target_colors")
axes[1].set_xlabel("PC1")
axes[1].set_ylabel("PC2")

# Aggiungi una legenda per i colori dei dati di test
custom_legend_complexity = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#000000", markersize=10, label="Complexity 1"),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#440154", markersize=10, label="Complexity 2"),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#3b528b", markersize=10, label="Complexity 3"),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#21918c", markersize=10, label="Complexity 4"),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#5ec962", markersize=10, label="Complexity 5"),
    plt.Line2D([0], [0], marker='o', color='w', alpha=0.3, markerfacecolor="orange", markersize=10, label="Assembled"),
    plt.Line2D([0], [0], marker='o', color='w', alpha=0.3, markerfacecolor="red", markersize=10, label="PATRIC"),
]
axes[0].legend(handles=custom_legend_complexity, title="Test Complexity Colors", loc="upper right")

custom_legend_resistance = [
    plt.Line2D([0], [0], marker='o', color='w', alpha=0.3, markerfacecolor="orange", markersize=10, label="Assembled"),
    plt.Line2D([0], [0], marker='o', color='w', alpha=0.3, markerfacecolor="red", markersize=10, label="PATRIC"),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="black", markersize=10, label="Dal viola (Suscettibile) al giallo (Resistente)"),
]
axes[1].legend(handles=custom_legend_resistance, title="Test Resistance Colors", loc="upper right")
# Mostra il grafico
plt.tight_layout()
plt.savefig("PCA_samples_features.png")


#%%
#PCA con i campioni che hanno 5 speci

