#%%
import json
import sys
import pandas as pd
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import r2_score
import numpy as np
import matplotlib.pyplot as plt
import M_learning_images as im
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, make_scorer
import os

def rmse(y_true, y_pred):
    return np.sqrt(mean_squared_error(y_true, y_pred))
rmse_scorer = make_scorer(rmse, greater_is_better=False)

print("Si parteeeee")

JUST_COMPLETE = False
NUMBER_OF_STRAINS = 1
MAPQ_MIN = 0
MAPQ_MAX = 1000
SCORING = "r2" # "r2" or "mean_squared_error"

#%%
campioni_PATRIC = pd.read_csv("Database_campioni2.csv", index_col=0)
campioni_PATRIC = campioni_PATRIC[campioni_PATRIC["Size"] == 0.1]
campioni_PATRIC = campioni_PATRIC[campioni_PATRIC["Real/Simul"] == "Simul"]
campioni_PATRIC = campioni_PATRIC[campioni_PATRIC["Aureus_S_N_strains"] == NUMBER_OF_STRAINS]
average_read_length_train = 270
assert all([x==average_read_length_train for x in campioni_PATRIC["reads_length"]])
Train_paths = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in campioni_PATRIC.index]
Train_data_PATRIC = im.read_many_alignements_json(Train_paths, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, mapq_average_max=MAPQ_MAX, read_length=average_read_length_train).T
print("Ho letto le features")
Train_target_PATRIC = dict((campioni_PATRIC["Aureus_R"] + campioni_PATRIC["Coli"] + campioni_PATRIC["Pseudomonas"]) /
				  (campioni_PATRIC["Arcobacter"] + campioni_PATRIC["Aureus_S"] + campioni_PATRIC["Aureus_R"] + campioni_PATRIC["Coli"] + campioni_PATRIC["Pseudomonas"]))

campioni_Assembled = pd.read_csv("Database_campioni_assembled.csv", index_col=0)
campioni_Assembled = campioni_Assembled[campioni_Assembled["Size"] == 0.1]
campioni_Assembled = campioni_Assembled[campioni_Assembled["Real/Simul"] == "Simul"]
campioni_Assembled = campioni_Assembled[campioni_Assembled["N_strains_for_all"] == NUMBER_OF_STRAINS]
average_read_length_train = 255
assert all([x==average_read_length_train for x in campioni_Assembled["reads_length"]])
Train_paths = ["samples_4/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in campioni_Assembled.index]
Train_data_Assembled = im.read_many_alignements_json(Train_paths, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, mapq_average_max=MAPQ_MAX, read_length=average_read_length_train).T
Train_target_Assembled = dict((campioni_Assembled["Aureus_R"] + campioni_Assembled["Coli"] + campioni_Assembled["Pseudomonas"] +
					   campioni_Assembled["Aureus_R_1"] + campioni_Assembled["Coli_1"] + campioni_Assembled["Pseudomonas_1"]) / 
					   (campioni_Assembled["Arcobacter"] + campioni_Assembled["Aureus_S"] + campioni_Assembled["Aureus_R"] + campioni_Assembled["Coli"] + campioni_Assembled["Pseudomonas"] + 
		 				campioni_Assembled["Arcobacter_1"] + campioni_Assembled["Aureus_S_1"] + campioni_Assembled["Aureus_R_1"] + campioni_Assembled["Coli_1"] + campioni_Assembled["Pseudomonas_1"]))
print("Ho fatto il trian set, features e target")

#%%
campioni_per_test = pd.read_csv("Database_campioni2.csv", index_col=0)
t_test = campioni_per_test[campioni_per_test["Size"] == 0.1]
t_test = t_test[t_test["Real/Simul"] == "Real"]
average_read_length_test = 250
assert all([x==average_read_length_test for x in t_test["reads_length"]])
t_test = t_test[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
Test_paths = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_test.index]
Test_data = im.read_many_alignements_json(Test_paths, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, read_length=average_read_length_test)
Test_data = Test_data.T
Test_target = dict(t_test[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1) / t_test[["Arcobacter", "Aureus_S", "Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1))	
select_target_complexity = t_test != 0
select_target_complexity = select_target_complexity.astype(int)
complexity_score = select_target_complexity.sum(axis=1)

#Fa i colori
complexity = 3
complexity_combinations = pd.DataFrame(select_target_complexity[select_target_complexity.sum(axis=1) == complexity])
complexity_combinations["combination"] = complexity_combinations.apply(lambda row: '_'.join(row.astype(str)), axis=1)
unique_combinations = complexity_combinations['combination'].unique()
colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_combinations)))
#make colors in a readable format


#%%
################
#	 Print sum of all resistance genes in the training data
################
somma_PARTIC = Train_data_PATRIC.sum().sum()
somma_Assembled = Train_data_Assembled.sum().sum()
print("Train data PATRIC sum: ", somma_PARTIC)
print("Train data Assembled sum: ", somma_Assembled)
print("Difference: ", somma_PARTIC - somma_Assembled)
################
#     Plot target distribution in histograms
################
plt.figure()
plt.hist(list(Train_target_PATRIC.values()), bins=20, alpha=0.5, label="PATRIC")
plt.hist(list(Train_target_Assembled.values()), bins=20, alpha=0.5, label="Assembled")
plt.legend()
plt.title("Train target distribution")
plt.show()
################
print("Proportion of samples with 0 resistance in the training data, expected 0.1 with 5 genomes")
print("PATRIC: ", len([x for x in Train_target_PATRIC.values() if x == 0]) / len(Train_target_PATRIC))
print("Assembled: ", len([x for x in Train_target_Assembled.values() if x == 0]) / len(Train_target_Assembled))
print("Proporzion of samples with 1 resistance in the training data, expected 0.2 with 5 genomes")
print("PATRIC: ", len([x for x in Train_target_PATRIC.values() if x == 1]) / len(Train_target_PATRIC))
print("Assembled: ", len([x for x in Train_target_Assembled.values() if x == 1]) / len(Train_target_Assembled))
print("Nel caso dei campiolni assemblati ho in realtà 10 genomi e quindi cambi al proporzione, probabilemente di un fattore 2")
#%%
################ RESISTANGE GENES IN THE TRAINING AND TEST SAMPLES
#	 Plot an histogram with the row sums of Train_data_PATRIC and Train_data_Assembled and Test_data
################
bins = np.histogram_bin_edges(np.concatenate((Train_data_PATRIC.sum(axis=1), Train_data_Assembled.sum(axis=1), Test_data.sum(axis=1))), bins=100)
plt.figure()
plt.hist(Train_data_PATRIC.sum(axis=1), bins=bins, alpha=0.5, label="PATRIC")
plt.hist(Train_data_Assembled.sum(axis=1), bins=bins, alpha=0.5, label="Assembled")
plt.hist(Test_data.sum(axis=1), bins=bins, alpha=0.8, label="Test", color='green')
plt.legend()
plt.title("Distribuzione della somma delle righe nei dati di addestramento")
plt.show()
################
# 	Vedo qualche outlier nei dati di test
print("Outliers in test data:")
print(Test_data.sum(axis=1)[Test_data.sum(axis=1) > 10000])
print("Gli outliers sono i campioni con solo escherigha coli")
# %%
################ GENES DISTRIBUTION IN THE SAMPLES (Normalized) -> SI NOTA UN OUTLIER
#     Plot sum of each column in Train_data, sorting them from the most to the least
#     Plot sum of each column in Assembled_data, sorting them from the most to the least
#     Plot sum of each column in Test_data, sorting them from the most to the least
#     Use log scale on y axis
################
Gene_abu_P = Train_data_PATRIC.sum(axis=0) / Train_data_PATRIC.shape[0]
Gene_abu_ass = Train_data_Assembled.sum(axis=0) / Train_data_Assembled.shape[0]
Gene_abu_Test = Test_data.sum(axis=0) / Test_data.shape[0]
bins = np.histogram_bin_edges(np.concatenate((Gene_abu_P, Gene_abu_ass, Gene_abu_Test)), bins=100)
plt.figure()
plt.hist(Gene_abu_P, bins=bins, alpha=0.5, label="PATRIC")
plt.hist(Gene_abu_ass, bins=bins, alpha=0.5, label="Assembled")
plt.hist(Gene_abu_Test, bins=bins, alpha=0.5, label="Test")
plt.yscale("log")
plt.legend()
plt.title("Abbondanza diei geni nei database")
plt.show()
################
#    Ho trovato un outlier nei dati di patric e 3 outliers nei dati di test
#	 voglio vedere come è distribuito nei campioni di Train (Assembled e PATRIC) e test set
################
#outlier = Gene_abu_P[Gene_abu_P > 80000].index
outlier_PATRIC = [x for x in Gene_abu_P[Gene_abu_P > 200].index]
print("Outliers in patric:")
print(Gene_abu_P[Gene_abu_P > 200])
outliers_Test = [x for x in Gene_abu_Test[Gene_abu_Test > 200].index]
print("Outliers in Test:")
print(Gene_abu_Test[Gene_abu_Test > 200])
for outlier in outliers_Test:
	vai = True
	all_data = Test_data[outlier]
	if outlier in Train_data_PATRIC.columns: all_data = np.concatenate([all_data, Train_data_PATRIC[outlier]])
	if outlier in Train_data_Assembled.columns: all_data = np.concatenate([all_data, Train_data_Assembled[outlier]])
	fig, ax = plt.subplots()
	bins = np.histogram_bin_edges(all_data, bins=100)
	ax.hist(Test_data[outlier], bins=bins, alpha=0.9, label="Test", color="green")
	if outlier in Train_data_Assembled.columns: ax.hist(Train_data_Assembled[outlier], bins=bins, alpha=0.5, label="Assembled")
	if outlier in Train_data_PATRIC.columns: ax.hist(Train_data_PATRIC[outlier], bins=bins, alpha=0.5, label="PATRIC")
	ax.legend()
	ax.set_title("Distribuzione di"+outlier+" nei dati di addestramento")
	plt.show()
	#######################################
	# Plot the outlier distribution in the test set.
	# Create a barplot with a bar for each sample of the test set
	outlier_test = Test_data[outlier][:23]
	fig, ax = plt.subplots()
	outlier_test.plot(kind='bar', ax=ax)
	ax.set_title("Distribuzione degli outlier nel set di test")
	ax.set_xlabel("Campioni del set di test")
	ax.set_ylabel("Somma degli outlier")
	plt.show()

#%%
############### 	SOMMA DEI GENI RESISTENZA NEI CAMPIONI
#	 Plot the sum of the resistance genes in the samples
###############
bins = np.histogram_bin_edges(np.concatenate((Train_data_PATRIC.sum(axis=1), Train_data_Assembled.sum(axis=1), Test_data.sum(axis=1))), bins=100)
fig, ax = plt.subplots()
Train_data_PATRIC.sum(axis=1).plot(kind='hist',bins=bins, ax=ax, alpha=0.5, label="PATRIC")
Train_data_Assembled.sum(axis=1).plot(kind='hist',bins=bins, ax=ax, alpha=0.5, label="Assembled")
Test_data.sum(axis=1).plot(kind='hist',bins=bins, ax=ax, alpha=0.5, label="Test")
ax.legend()
ax.set_title("Somma dei geni di resistenza nei campioni")
plt.show()

# %%
############ Vediamo i campioni islati
TDP = campioni_PATRIC.drop(["Size", "Real/Simul", "Original sample", "Seed", "Aureus_S_N_strains", "Aureus_R_N_strains", "Coli_N_strains", "Pseudomonas_N_strains", "Arcobacter_N_strains", "reads_length"], axis=1)
TDA = campioni_Assembled.drop(["Size", "Real/Simul", "Original sample", "Seed", "N_strains_for_all", "reads_length"], axis=1)
for isolate in TDP.columns:
	TDP_ = TDP[TDP[isolate] == 100]
	TDA_ = TDA[[isolate, isolate+"_1"]].sum(axis=1) == 100
	TDA_ = TDA[TDA_]
	Train_data_PATRIC_ = Train_data_PATRIC.loc[TDP_.index]
	Train_data_Assembled_ = Train_data_Assembled.loc[TDA_.index]
	Gene_abu_P_ = Train_data_PATRIC_.sum(axis=0) / Train_data_PATRIC_.shape[0]
	Gene_abu_A_ = Train_data_Assembled_.sum(axis=0) / Train_data_Assembled_.shape[0]
	bins = np.histogram_bin_edges(np.concatenate((Gene_abu_P_, Gene_abu_A_)), bins=100)
	fig, ax = plt.subplots()
	ax.hist(Gene_abu_P_, bins=bins, alpha=0.5, label="PATRIC")
	ax.hist(Gene_abu_A_, bins=bins, alpha=0.5, label="Assembled")
	ax.set_yscale("log")
	ax.legend()
	ax.set_title("Abbondanza dei geni nei campioni isolati da "+isolate)
	plt.show()
	########################
	# print TDP_ and TDA_ most abundant genes
	########################
	print("Most abundant genes in "+isolate+" in PATRIC:")
	print(Gene_abu_P_.sort_values(ascending=False).head(10))
	print("Most abundant genes in "+isolate+" in Assembled:")
	print(Gene_abu_A_.sort_values(ascending=False).head(10))
# %%
