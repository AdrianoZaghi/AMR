import json
import sys
import pandas as pd
#from sklearn.linear_model import ElasticNet
#from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import GridSearchCV
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
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1, 0.5, 1, 2]}
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
MAPQ_MIN = -1
MAPQ_MAX = 1000
average_read_length = 270
#SCORING = ["RMSE", "R2"]
SCORING = ["R2"]
TARGET = ["Pathogens", "Resistant", "Resistant_in_pathogens"]
depths = [0.001, 0.003, 0.007, 0.01, 0.05, 0.1, 0.5, 1]
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
			for i in range(0,1):
				#Train_data = im.read_bwt_2_aro(Train, sub)
				Train_data = im.read_many_alignements_json(Train_features_path, 
														just_completely_mapped=JUST_COMPLETE,
														mapq_average_min=MAPQ_MIN,
														mapq_average_max=MAPQ_MAX,
														read_length=average_read_length,
														subsample=sub).T
				if scoring == "R2":
					grid = GridSearchCV(model, params, n_jobs=-1, refit=False, cv=10, verbose=0, scoring="r2")
					score_sign = 1
				elif scoring == "RMSE":
					grid = GridSearchCV(model, params, n_jobs=-1, refit=False, cv=10, verbose=0, scoring=rmse_scorer)
					score_sign = -1
				grid.fit(Train_data, list(Train_target))
				train_average_score = score_sign * np.mean(grid.cv_results_["mean_test_score"])
				print("Score: "+str(train_average_score))
				train_std_score = np.std(grid.cv_results_["mean_test_score"])
				Results["Subsample "+str(sub)] = {scoring:train_average_score, scoring + "_std":train_std_score}

		save_path = "Sembra_che_funzioni_multispecie/"+scoring+"_depth_effects_"+target
		pd.DataFrame(Results).T.to_csv(save_path+".csv")
		Scores = [v[scoring] for _,v in Results.items()]
		Scores_std = [v[scoring+"_std"] for _,v in Results.items()]

		fig, ax = plt.subplots(figsize=(10, 6))
		# Crea il grafico a barre con le barre di errore
		ax.errorbar(depths, Scores, yerr=Scores_std, fmt='o')
		ax.set_xscale('log')
		ax.axhline(y=1, color='r', linestyle='--')
		if scoring == "R2":
			ax.set_ylim(0, 1.1)
		ax.set_title("Effect of sample depth")
		ax.set_xlabel('Subsample size (Gbp)')
		ax.set_ylabel(scoring)

		# Salva il grafico
		plt.savefig(save_path + ".png")