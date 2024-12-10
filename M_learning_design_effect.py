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
import random

no_bias = {'Aureus_S': 1,
	'Aureus_R': 1,
 	'Coli': 1,
 	'Pseudomonas': 1,
 	'Arcobacter': 1}

print("Si parteeeee")
#Lunghezza media delle reads simulate e reali, comunque un fattore di normalizzazione che dovrebbe valere per entrambi
#potrebbe essere diverso tra simulate e reali
images_dir = "M_learning_images_4/"
campioni = pd.read_csv("Database_campioni_multispecie.csv", index_col=0)

Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
Ada =  AdaBoostRegressor(estimator=Base_estimator_for_adaboost,n_estimators=100)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1, 0.5, 1, 2]}#SELECTING THE MODEL TO USE

model = Ada
params = params_Ada
subsample_counts = 1

WHICH = "Perfected"
JUST_COMPLETE = True
NUMBER_OF_STRAINS = 1
MAPQ_MIN = -1
genes_to_take_out = []
print(WHICH)
#Deciding which samples to peack from the database and defining the train and test sets
t = campioni[campioni["Size"] == 0.1]
#Selecting train samples
#t_train = t[t["Real/Simul"] == "Simul"]
#t_train = t[t["Real/Simul"] == "De_Caesar_simul"]
t_train = t[t["Real/Simul"] == "Simul_perfected"]
t_train = t_train[t_train["Aureus_S_N_strains"] == NUMBER_OF_STRAINS]
Train = ["samples_3/" + x for x in t_train.index]
Train_data = im.read_bwt_2_aro(Train, subsample_counts, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, av_read_length=255)
Train_data = Train_data.drop(columns=[x for x in genes_to_take_out if x in Train_data.columns])
print("Ho letto le features")
Train_target = im.read_metadata(Train)["Abundances_r"]
print("Ho letto tuttooooooo")

#Preparing the test samples for selection
#t_test = t[t["Real/Simul"] == "Real"]
t_test = t[t["Real/Simul"] == "Real_new_trim"]
t_test = t_test[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
select_target_complexity = t_test != 0
select_target_complexity = select_target_complexity.astype(int)

print("Entriamo nel looooop")
Results = {}
for i in range(2,6):
	ttt = t_test[select_target_complexity.sum(axis=1) == i]
	Test = ["samples_2/" + x for x in ttt.index]
	Test_data = im.read_bwt_2_aro(Test, subsample_counts, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, av_read_length=245)
	Test_data = Test_data.drop(columns=[x for x in genes_to_take_out if x in Test_data.columns])
	if True:#adding colors
		ttc = pd.DataFrame(select_target_complexity[select_target_complexity.sum(axis=1) == i])
		tot_row_sums = Test_data.sum(axis=1)
		Test_data_ = Test_data.copy()
		Test_data_.index = [x.split("/")[1] for x in Test_data.index]
		row_sums = Test_data_.loc[ttc.index].sum(axis=1)
		norm = plt.Normalize(tot_row_sums.min(), tot_row_sums.max())
		normalized_row_sums = norm(row_sums)
		colors = plt.cm.rainbow(normalized_row_sums)
		colors_mapped = pd.Series(dict(zip(t_test.index, colors)))
	if False:#adding colors
		ttc = pd.DataFrame(select_target_complexity[select_target_complexity.sum(axis=1) == i])
		ttc["combination"] = ttc.apply(lambda row: '_'.join(row.astype(str)), axis=1)
		unique_combinations = ttc['combination'].unique()
		colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_combinations)))
		color_map = dict(zip(unique_combinations, colors))
		colors_mapped = ttc['combination'].map(color_map)
		#print(color_map)
	print("Sto facendo il fit per "+str(i)+" specie")
	#Sto prendendo il target dai metadata
	Test_target_nominale = im.read_metadata(Test)["Abundances_r"]
	ttt = ttt[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
	Test_target_aggiustato = dict(zip(Test, ttt))
	#Dealing with columns absent in one of the two datasets 
	if False: #Drop columns which are not common to both datasets	
		missing_in_Train = [c for c in Test_data.columns if c not in Train_data.columns]
		missing_in_Test = [c for c in Train_data.columns if c not in Test_data.columns]
		Test_data_clean = Test_data.drop(columns=missing_in_Train).sort_index(axis=1)
		Train_data_clean = Train_data.drop(columns=missing_in_Test).sort_index(axis=1)
	else: #Add columns which are missing in one of the two datasets filling with 0
		all_columns = list(set(Test_data.columns).union(set(Train_data.columns)))
		Test_data_clean = Test_data.reindex(columns = all_columns, fill_value=0).sort_index(axis=1)
		Train_data_clean = Train_data.reindex(columns = all_columns, fill_value=0).sort_index(axis=1)
	print("Sto usando "+str(len(Train_data_clean.columns))+" colonne nel fit")
	R2s_for_species_number = []
	R2s_train_for_species_number_adj = []
	#Per ogni numero di specie, faccio 5 fit diversi per ottenere dei valori di R2 piu' robusti
	for iii in range(0,1):
		grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=10, verbose=0, scoring="r2", )
		grid.fit(Train_data_clean, list(Train_target.values()))
		train_average_score = np.mean([np.mean(v) for k,v in grid.cv_results_.items() if "split" in k])
		print("Train_set_av_score   "+str(train_average_score)+"   .")
		best = grid.best_estimator_
		if False:
			model.fit(Train_data_clean, list(Train_target.values()))
			Test_predictions = model.predict(Test_data_clean)
		elif True:
			Test_predictions = grid.predict(Test_data_clean)
		Test_predictions = dict(zip(Test_data_clean.index, Test_predictions))
		R2 = r2_score(list(Test_target_aggiustato.values()), list(Test_predictions.values())) #il R2 ridotto non mi sembra adatto, non ho un numero preciso di parametri che fitto
		print("R2   "+str(R2)+"   .")
		im.plot_results_2(list(Test_target_aggiustato.values()), list(Test_predictions.values()), colors=colors_mapped, legend = None, filename = images_dir+WHICH+"_complexity_" + str(i)+ "_" + str(iii) + ".png", R2 = R2, y_low_lim=0, y_up_lim=1, x_low_lim=0, x_up_lim=1)
		R2s_for_species_number.append(R2)
		R2s_train_for_species_number_adj.append(train_average_score)
	Results["Number of species: " + str(i)] = {"R2" : R2s_for_species_number, "R2_train" : R2s_train_for_species_number_adj}

R2s = [[int(k.split("Number of species: ")[1]), v["R2"], v["R2_train"]] for k,v in Results.items()]
R2s = pd.DataFrame(R2s, columns = ["Number of species", "R2", "R2_train"])
R2s["Average R2"] = R2s["R2"].apply(np.mean)
R2s["R2 errors"] = R2s["R2"].apply(np.std)
R2s["Average R2 train"] = R2s["R2_train"].apply(np.mean)
R2s["R2 errors train"] = R2s["R2_train"].apply(np.std)
R2s.to_csv(images_dir+WHICH+"_R2s_complexity_effects.csv")

fig, ax = plt.subplots()
ax.errorbar(R2s['Number of species'], R2s['Average R2'], yerr=R2s['R2 errors'], fmt='o')
ax.set_xlabel('Sample complexity')
ax.set_xticks([2,3,4,5])
ax.axhline(y=1, color='r', linestyle='--')
ax.set_ylabel('R2')
ax.set_title("Effect of sample complexity")
#ax.set_ylim(0.5,1.1)


# Show the plot
plt.savefig(images_dir+"/Complexity_effect_"+WHICH+"_synth.png")