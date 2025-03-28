import json
import sys
import pandas as pd
from sklearn.ensemble import AdaBoostRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.metrics import r2_score
import numpy as np
import matplotlib.pyplot as plt
import M_learning_images as im
from sklearn.model_selection import train_test_split
import seaborn as sns
import random
import os


#Lunghezza media delle reads simulate e reali, comunque un fattore di normalizzazione che dovrebbe valere per entrambi
#potrebbe essere diverso tra simulate e reali
average_read_length = 270

Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
Ada =  AdaBoostRegressor(estimator=Base_estimator_for_adaboost,n_estimators=500, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1, 0.5, 1, 2]}
#SELECTING THE MODEL TO USE

model = Ada
params = params_Ada

#####################################
#		Multispecie - Pathogen	   	Fatto
#		Multispecie - Resistant	   	Sto facendo
#		Single specie - 1 strain	Da fare
#		Single specie - 5 strain	Da fare

STRAINS = 1
JUST_COMPLETE = True
MAPQ_MIN = 0
MAPQ_MAX = 1000
Che_training_set = "Single specie Assembled"
TARGET_MULTISPECIE = "Resistant"
save_path = "Feature_imp"
if os.path.exists(save_path) == False:
	os.makedirs(save_path)

if Che_training_set == "Multispecie":
	save_path = save_path + "/Multispecie/"
	if os.path.exists(save_path) == False: os.makedirs(save_path)
	campioni = pd.read_csv("Database_campioni_multispecie.csv", index_col=0)
	t = campioni[campioni["Size"] == 0.1]
	assert STRAINS == 1, "Strains must be 1 for multispecie"
	t = t[t["N_strains_for_all"] == STRAINS]
	Train = ["samples_3/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t.index]
	Target = t.drop(columns=["Size", "N_strains_for_all", "Real/Simul", "reads_length", "Seed", "Original sample"])
	Target = Target.div(Target.sum(axis=1), axis=0)
	if TARGET_MULTISPECIE == "Pathogen":
		save_path = save_path + "_pathogen/"
		Target = Target[["Aureus_R", "Aureus_S", "Coli", "Pseudomonas", "Arcobacter"]].sum(axis=1)
	elif TARGET_MULTISPECIE == "Resistant":
		save_path = save_path + "_resistant/"
		Target = Target[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
	else: raise ValueError("Target must be Pathogen or Resistant")
	Train_target = dict(zip(Target.index, Target.values))
elif Che_training_set == "Single specie PATRIC":
	save_path = save_path + "/Single_specie_PATRIC/"
	if os.path.exists(save_path) == False: os.makedirs(save_path)
	campioni = pd.read_csv("Database_campioni2.csv", index_col=0)
	t = campioni[campioni["Size"] == 0.1]
	t = t[t["Aureus_S_N_strains"] == STRAINS]
	t = t[t["Real/Simul"] == "Simul"]
	Train = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t.index]
	Target = t[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
	Target = Target.div(Target.sum(axis=1), axis=0) #non si sa mai
	Target = Target[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
	Train_target = dict(zip(Target.index, Target.values))
elif Che_training_set == "Single specie Assembled":
	save_path = save_path + "/Single_specie_Assembled/"
	if os.path.exists(save_path) == False: os.makedirs(save_path)
	campioni = pd.read_csv("Database_campioni_assembled.csv", index_col=0)
	t = campioni[campioni["Size"] == 0.1]
	t = t[t["N_strains_for_all"] == STRAINS]
	t = t[t["Real/Simul"] == "Simul"]
	Train = ["samples_4/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t.index]
	Target = t[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter", "Aureus_S_1", "Aureus_R_1", "Coli_1", "Pseudomonas_1", "Arcobacter_1"]]
	Target = Target.div(Target.sum(axis=1), axis=0) #non si sa mai
	Target = Target[["Aureus_R", "Coli", "Pseudomonas", "Aureus_R_1", "Coli_1", "Pseudomonas_1"]].sum(axis=1)
	Train_target = dict(zip(Target.index, Target.values))

print(save_path)

Train_data = im.read_many_alignements_json(Train,
											 just_completely_mapped=JUST_COMPLETE,
											 mapq_average_min=MAPQ_MIN,
											 mapq_average_max=MAPQ_MAX,
											 read_length=average_read_length,
											 subsample=1).T
print("Dati letti.\nSto usando "+str(len(Train_data.columns))+" colonne nel fit")
Feature_importances_total = []
print("Letto tutto, ora faccio il fit")
###################################################
#  		ALTRO MODO IN CUI POTREI FARE
#		Divido i dati in 10 fold all'inizio
#		Per ogni 9-fold set, faccio il fit.
#		    In questo caso, comunque ottimizzo i parametri nei 9 fold con una ulteriore cross validation
#		Osservo la distribuzione delle feature importance di questi 10 modselli
###################################################
for i in range(5):
	model.random_state = i
	print(i)
	cross_val_folds = KFold(n_splits=10, shuffle=True, random_state=i)
	grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=cross_val_folds, verbose=0, scoring="r2")
	print("Start looking for best hyperparams")
	grid.fit(Train_data, list(Train_target.values()))
	F_imp = grid.best_estimator_.feature_importances_
	F = Train_data.columns
	F_imp = {F[i]:F_imp[i] for i in range(len(F))}
	Feature_importances_total.append(F_imp)

Features = pd.DataFrame(Feature_importances_total)

WRITING_OR_READING = "w"
if WRITING_OR_READING == "w":
	Features.T.to_csv(save_path + "/Feature_importances.csv")
elif WRITING_OR_READING == "r":
	path = save_path+ "/Feature_importances.csv"
	Features = pd.read_csv(path, index_col=0).T
How_many_cols = 30
Sorted_col = list(Features.mean().sort_values(ascending=False).index)
Just_the_col_i_need = Sorted_col[:How_many_cols]
Just_the_col_i_need = Just_the_col_i_need[::-1]
sorted_df = Features[Just_the_col_i_need]
for name in sorted_df.columns:
	if "ARO" in name.split("|")[0]:
		new_name = "ARO_ID:" + name.split("|")[0].split(":")[1]
	else:
		new_name = "ARO_ID:" + name.split("|")[3].split(":")[1] + "|Wildcard_seq:" + name.split("|")[0].split(":")[1]
	sorted_df.rename(columns={name: new_name}, inplace=True)
fig, ax = plt.subplots(figsize=(20, 12))
sorted_df.boxplot(ax=ax, vert=False)
ax.set_title('Feature importance assigned to CARD genes', fontsize=15)
ax.set_xlabel('Feature Importance', fontsize=18)
ax.set_ylabel('CARD gene identifier', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.savefig(save_path + "/Top_30_features.png")