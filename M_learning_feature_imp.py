#%%
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
import random
import os
import sys


#Lunghezza media delle reads simulate e reali, comunque un fattore di normalizzazione che dovrebbe valere per entrambi
#potrebbe essere diverso tra simulate e reali
average_read_length = 270

Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
Ada =  AdaBoostRegressor(estimator=Base_estimator_for_adaboost,n_estimators=100, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[1]}
#SELECTING THE MODEL TO USE

model = Ada
params = params_Ada

#####################################
#		Multispecie - Pathogen	   	Fatto
#		Multispecie - Resistant	   	Sto facendo
#		Single specie - 1 strain	Da fare
#		Single specie - 5 strain	Da fare

if len(sys.argv) == 1:
	Che_training_set = "Single specie PATRIC"
elif len(sys.argv) > 2:
	assert sys.argv[1] in ["Multispecie", "Single specie PATRIC"], "First argument must be Multispecie or Single specie PATRIC"
	Che_training_set = sys.argv[1]
	TARGET_MULTISPECIE = sys.argv[2]
else:
	raise ValueError("You must provide no arguments (For single specie PATRIC) or two arguments (For multispecie).")

JUST_COMPLETE = True
MAPQ_MIN = 0
MAPQ_MAX = 1000

FILTRO_ABBONDANZA = 6
FILTRO_PREVALENZA = 0.1
save_path = "New/Feature_imp"
if os.path.exists(save_path) == False:
	os.makedirs(save_path)

if Che_training_set == "Multispecie":
	FILTRO_PREVALENZA = 0.07
	campioni = pd.read_csv("Database_campioni_multispecie.csv", index_col=0)
	t = campioni[campioni["Size"] == 0.1]
	t = t[t["N_strains_for_all"] == 1]
	Train = ["samples_3/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t.index]
	Target = t.drop(columns=["Size", "N_strains_for_all", "Real/Simul", "reads_length", "Seed", "Original sample"])
	Target = Target.div(Target.sum(axis=1), axis=0)
	save_path = save_path + "/Multispecie"
	if TARGET_MULTISPECIE == "Pathogen":
		save_path = save_path + "_pathogen/"
		Target = Target[["Aureus_R", "Aureus_S", "Coli", "Pseudomonas", "Arcobacter"]].sum(axis=1)
	elif TARGET_MULTISPECIE == "Resistant":
		save_path = save_path + "_resistant/"
		Target = Target[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
	else: raise ValueError("Target must be Pathogen or Resistant")
	if os.path.exists(save_path) == False: os.makedirs(save_path)
	Train_target = dict(zip(Target.index, Target.values))
elif Che_training_set == "Single specie PATRIC":
	save_path = save_path + "/Single_specie_PATRIC/"
	if os.path.exists(save_path) == False: os.makedirs(save_path)
	campioni = pd.read_csv("Database_campioni2.csv", index_col=0)
	t = campioni[campioni["Size"] == 0.1]
	t = t[t["Aureus_S_N_strains"] == 1]
	t = t[t["Real/Simul"] == "Simul"]
	Train = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t.index]
	Target = t[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
	Target = Target.div(Target.sum(axis=1), axis=0) #non si sa mai
	Target = Target[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
	Train_target = dict(zip(Target.index, Target.values))
	t_test = pd.read_csv("Database_campioni2.csv", index_col=0)
	t_test = t_test[t_test["Size"] == 0.1]
	t_test = t_test[t_test["Real/Simul"] == "Real"]
	campioni_veri_che_non_usiamo = ["DCWGS"+i+"_" for i in ["1", "2", "3", "4", "5", "6", "7", "AR", "PA", "EC", "SAS", "SAR"]]
	for i in campioni_veri_che_non_usiamo:
		t_test = t_test[~t_test.index.str.contains(i)]
	Test_paths = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_test.index]


print(save_path)

Train_data = im.read_many_alignements_json_1(Train,
											 just_completely_mapped=JUST_COMPLETE,
											 mapq_average_min=MAPQ_MIN,
											 mapq_average_max=MAPQ_MAX,
											 read_length=average_read_length,
											 subsample=1,
											 filtro_abbondanza=FILTRO_ABBONDANZA).T

col_prevalence = (Train_data > 0).sum(axis=0)/ len(Train_data)
Train_data = Train_data.loc[:, col_prevalence > FILTRO_PREVALENZA]

if Che_training_set == "Single specie PATRIC":
	Test_data = im.read_many_alignements_json_1(Test_paths,
												just_completely_mapped=JUST_COMPLETE,
												mapq_average_min=MAPQ_MIN,
												mapq_average_max=MAPQ_MAX,
												read_length=average_read_length,
												subsample=1,
												filtro_abbondanza=FILTRO_ABBONDANZA).T
	column_intersection = Train_data.columns.intersection(Test_data.columns)
	Train_data = Train_data[column_intersection]

print("Dati letti.\nSto usando "+str(len(Train_data.columns))+" colonne nel fit")
Feature_importances_total = []
print("Letto tutto, ora faccio il fit")

#%%
cross_val_folds = KFold(n_splits=10, shuffle=True, random_state=0)
grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=cross_val_folds, verbose=0, scoring="r2")
print("Start looking for best hyperparams")
grid.fit(Train_data, list(Train_target.values()))
best_params = grid.best_params_
print("Best params: ", best_params)
folds_spliter = KFold(n_splits=5, shuffle=True, random_state=0)
new_folds = folds_spliter.split(Train_data)
for train, test in new_folds:
	print("Train")
	data_folds = Train_data.iloc[train]
	target_folds = [Train_target[idx] for idx in data_folds.index]
	data_test = Train_data.iloc[test]
	target_test = [Train_target[idx] for idx in data_test.index]
	model.set_params(**best_params)
	model.fit(data_folds, target_folds)
	F_imp = model.feature_importances_
	F = Train_data.columns
	F_imp = {F[i]:F_imp[i] for i in range(len(F))}
	Feature_importances_total.append(F_imp)
Features = pd.DataFrame(Feature_importances_total)

Features.T.to_csv(save_path + "/Feature_importances.csv")
#path = save_path+ "/Feature_importances.csv"
#Features = pd.read_csv(path, index_col=0).T
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