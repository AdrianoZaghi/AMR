import json
import sys
import pandas as pd
import os
from sklearn.ensemble import AdaBoostRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import GridSearchCV, KFold, cross_val_predict
from sklearn.metrics import r2_score
import numpy as np
import matplotlib.pyplot as plt
import M_learning_images as im
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, make_scorer

def rmse(y_true, y_pred):
    return np.sqrt(mean_squared_error(y_true, y_pred))
rmse_scorer = make_scorer(rmse, greater_is_better=False)

#Lunghezza media delle reads simulate e reali, comunque un fattore di normalizzazione che dovrebbe valere per entrambi
#potrebbe essere diverso tra simulate e reali
campioni = pd.read_csv("Database_campioni_multispecie.csv", index_col=0)

Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
Ada =  AdaBoostRegressor(estimator=Base_estimator_for_adaboost,n_estimators=100, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1]}
#SELECTING THE MODEL TO USE

model = Ada
params = params_Ada

t = campioni[campioni["Size"] == 0.1]
Train = ["samples_3/" + x for x in t[t["Real/Simul"] == "Simul_perfected"].index]
Train_features_path= [t+"/bwt_out/bwt.allele_mapping_data.json" for t in Train]

Train_target = im.read_metadata_multispecie(Train)
Train_target_pathogen_abundance = [x["Pathogen_abundance"] for x in Train_target.values()]
Train_target_resistant_abundance = [x["Resistant_abundance"] for x in Train_target.values()]
Train_target_resistant_in_pathogens = [x["Resistant_in_pathogens"] for x in Train_target.values()]

JUST_COMPLETE = True
MAPQ_MIN = 0
MAPQ_MAX = 1000
average_read_length = 270
SCORING = ["R2"]
TARGET = ["Pathogens", "Resistant"]
depths = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1] # for testing purposes
FILTRO_ABBONDANZA = 6
FILTRO_PREVALENZA = 0.07
saving_path = "New/Depth_multispecie_prev_filter_"+str(FILTRO_PREVALENZA)+"/"
if not os.path.exists(saving_path):
	os.makedirs(saving_path)
Fold_Number = 10
for scoring in SCORING:
	for target in TARGET:
		print("Target: "+target + " Scoring: "+scoring)
		if target == "Pathogens":
			Train_target = Train_target_pathogen_abundance
		elif target == "Resistant":
			Train_target = Train_target_resistant_abundance
		elif target == "Resistant_in_pathogens":
			Train_target = Train_target_resistant_in_pathogens
		Results = {}
		for sub in depths:
			print("Depth: "+str(sub))
			#Train_data = im.read_bwt_2_aro(Train, sub)
			Train_data = im.read_many_alignements_json_1(Train_features_path, 
													just_completely_mapped=JUST_COMPLETE,
													mapq_average_min=MAPQ_MIN,
													mapq_average_max=MAPQ_MAX,
													read_length=average_read_length,
													subsample=sub,
													filtro_abbondanza=FILTRO_ABBONDANZA).T
			col_prevalence = (Train_data > 0).sum(axis=0)/ len(Train_data)
			Train_data = Train_data.loc[:, col_prevalence > FILTRO_PREVALENZA]
			if scoring == "R2":
				grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=Fold_Number, verbose=0, scoring="r2")
				score_sign = 1
			elif scoring == "RMSE":
				grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=Fold_Number, verbose=0, scoring=rmse_scorer)
				score_sign = -1
			print("Start looking for best hyperparams")
			grid.fit(Train_data, list(Train_target))
			best_model = grid.best_estimator_
			cross_val_del_best = KFold(n_splits=Fold_Number, shuffle=True, random_state=0)
			print("Start cross validation with best model")
			predictions = cross_val_predict(best_model, Train_data, list(Train_target), cv=cross_val_del_best, n_jobs=-1)
			#Per ottenere la media dello score, facio una media dei valori ottenuti con la cross validation
			print("Average_score",r2_score(list(Train_target), predictions))
			R2_per_fold = []
			for _, test_id in cross_val_del_best.split(Train_data):
				R2_per_fold.append(round(r2_score(np.array(Train_target)[test_id], np.array(predictions)[test_id]), 3))
			#Per ottenere la deviazione standard, faccio la deviazione standard dei valori ottenuti con la cross validation
			#La deviazione standard della media è la deviazione standard divisa per la radice del numero di fold???
			train_std_score = np.std(R2_per_fold)/np.sqrt(Fold_Number)
			Results["Subsample "+str(sub)] = {scoring:np.mean(R2_per_fold), scoring + "_std":train_std_score, "scores" : R2_per_fold}
			predictions = np.array(predictions)
			Train_target = np.array(Train_target)
			#Salvo le previsioni e i valori reali
			dati_scatterplot = pd.DataFrame({"Real":Train_target, "Predictions":predictions})
			dati_scatterplot.to_csv(saving_path + "scatterplot_data_" + target + "_depth_" + str(sub) + ".csv")
			fig, ax = plt.subplots(figsize=(8, 6))
			ax.scatter(Train_target, predictions, color='blue')
			ax.plot([0, 1], [0, 1], color='red', linestyle='--')  # Linea di riferimento y=x
			ax.set_xlabel('Valori Reali')
			ax.set_ylabel('Predizioni')
			ax.set_title('Scatterplot delle Predizioni vs Valori Reali')
			ax.grid(True)
			plt.savefig(saving_path + scoring + "_scatterplot_" + target + "_depth_" + str(sub) + ".png")
		save = saving_path+scoring+"_depth_effects_"+target
		pd.DataFrame(Results).T.to_csv(save+".csv")
		Scores = [v[scoring] for _,v in Results.items()]
		Scores_std = [v[scoring+"_std"] for _,v in Results.items()]
		fig, ax = plt.subplots(figsize=(10, 6))
		ax.errorbar(depths, Scores, yerr=Scores_std, fmt='o')
		ax.set_xscale('log')
		ax.axhline(y=1, color='r', linestyle='--')
		if scoring == "R2":
			ax.set_ylim(0, 1.1)
		ax.set_title("Effect of sample depth")
		ax.set_xlabel('Subsample size (Gbp)')
		ax.set_ylabel(scoring)
		plt.savefig(save + ".png")