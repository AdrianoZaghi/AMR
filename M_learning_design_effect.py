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
import seaborn as sns
import random

#Lunghezza media delle reads simulate e reali, comunque un fattore di normalizzazione che dovrebbe valere per entrambi
#potrebbe essere diverso tra simulate e reali
average_read_length = 270
images_dir = "M_learning_images/"
campioni = pd.read_csv("Database_campioni2.csv", index_col=0)

#%%
Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
Ada =  AdaBoostRegressor(estimator=Base_estimator_for_adaboost,n_estimators=500)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1]}
#SELECTING THE MODEL TO USE

model = Ada
params = params_Ada
subsample_counts = 1
USE_BIAS = True
WHICH = "_"

#Deciding which samples to peack from the database and defining the train and test sets
t = campioni[campioni["Size"] == 0.1]
#Selecting train samples
t_train = t[t["Real/Simul"] == "De_Caesar_simul"]
Train = ["samples_2/" + x for x in t_train.index]
Train_data = im.read_bwt_2_aro(Train, subsample_counts)
Train_target = im.read_metadata(Train)["Abundances_r"]

#Preparing the test samples for selection
t_test = t[t["Real/Simul"] == "Real"]
t_test = t_test[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
select_target_complexity = t_test != 0
select_target_complexity = select_target_complexity.astype(int)

#bias = {'Aureus_S': 6.37698087,
# 'Aureus_R': 8.46693667,
# 'Coli': 1,
# 'Pseudomonas': 1.83036887,
# 'Arcobacter': 1.8665934}
bias = {'Aureus_S': 4.93,
	'Aureus_R': 4.93,
 	'Coli': 1,
 	'Pseudomonas': 1,
 	'Arcobacter': 1}

bias_uncertanty = {'Aureus_S': 0.99561185,
'Aureus_R': 1.16214196,
'Coli': 0,
'Pseudomonas': 0.2528959,
'Arcobacter': 0.48897503}

Results = {}
for i in range(2,6):
	ttt = t_test[select_target_complexity.sum(axis=1) == i]
	Test = ["samples_2/" + x for x in ttt.index]
	Test_data = im.read_bwt_2_aro(Test, subsample_counts)
	#Sto prendendo il target dai metadata
	Test_target_nominale = im.read_metadata(Test)["Abundances_r"]
	if USE_BIAS:
		#Moltiplico tutte le abbondanze per il loro bias
		ttt = ttt.mul(bias, axis=1)
	#Normalizzo le abbondanze aggiustate e sommo le abbondanze di Aureus_R, Coli e Pseudomonas per ottenere la resistenza
	ttt = ttt.div(ttt.sum(axis=1), axis=0)[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
	Test_target_aggiustato = dict(zip(Test, ttt))
	Te_t = {k:v.set_index("Genome_ID")["Abundance"] for k,v in im.read_metadata(Test)["Metadata"].items()}
	#DEVO AGGIUSTARE STO BIASSSSS
	#Dealing with columns absent in one of the two datasets 
	missing_in_Train = [c for c in Test_data.columns if c not in Train_data.columns]
	missing_in_Test = [c for c in Train_data.columns if c not in Test_data.columns]
	Test_data_clean = Test_data.drop(columns=missing_in_Train).sort_index(axis=1)
	Train_data_clean = Train_data.drop(columns=missing_in_Test).sort_index(axis=1)
	print("Sto usando "+str(len(Train_data_clean.columns))+" colonne nel fit")
	R2s_for_species_number = []
	R2s_for_species_number_adj = []
	for iii in range(0,1):
		grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=10, verbose=0, scoring="r2", )
		grid.fit(Train_data_clean, list(Train_target.values()), )
		train_scores = {k:v for k,v in grid.cv_results_.items() if "split" in k}
		best = grid.best_estimator_
		Test_predictions = grid.predict(Test_data_clean)
		Test_predictions = dict(zip(Test_data_clean.index, Test_predictions))
		R2 = r2_score(list(Test_target_nominale.values()), list(Test_predictions.values())) #il R2 ridotto non mi sembra adatto, non ho un numero preciso di parametri che fitto
		#im.plot_results(list(Test_target_nominale.values()), list(Test_predictions.values()), "M_learning_images/" + WHICH + ".png", R2 = R2)
		#F_imp = im.feature_importance(Train_data, best, "M_learning_images/"+WHICH+"_feature_importance.png")
		#plt.clf()
		R2_adj = r2_score(list(Test_target_aggiustato.values()), list(Test_predictions.values())) #il R2 ridotto non mi sembra adatto, non ho un numero preciso di parametri che fitto
		im.plot_results(list(Test_target_aggiustato.values()), list(Test_predictions.values()), "M_learning_images/" + WHICH + ".png", R2 = R2_adj)
		#F_imp = im.feature_importance(Train_data, best, "M_learning_images/"+WHICH+"_feature_importance.png")
		plt.clf()
		R2s_for_species_number.append(R2)
		R2s_for_species_number_adj.append(R2_adj)
	Results["Number of species: " + str(i)] = {"R2" : R2s_for_species_number, "R2_adj" : R2s_for_species_number_adj}
	#Results["Number of species: " + str(i)]["Feature_importance"] = F_imp
	#Results["Number of species: " + str(i)]["Best_estimator"] = best
#%%
R2s = [[int(k.split("Number of species: ")[1]), v["R2"]] for k,v in Results.items()]
R2s = pd.DataFrame(R2s, columns = ["Number of species", "R2"])
R2s["Average R2"] = R2s["R2"].apply(np.mean)
R2s["R2 errors"] = R2s["R2"].apply(np.std)
R2s.to_csv("M_learning_images/"+WHICH+"_R2s_depth_effects.csv")
R2s_adj = [[int(k.split("Number of species: ")[1]), v["R2_adj"]] for k,v in Results.items()]
R2s_adj = pd.DataFrame(R2s_adj, columns = ["Number of species", "R2_adj"])
R2s_adj["Average R2"] = R2s_adj["R2_adj"].apply(np.mean)
R2s_adj["R2 errors"] = R2s_adj["R2_adj"].apply(np.std)
R2s_adj.to_csv("M_learning_images/"+WHICH+"_R2s_adj_depth_effects.csv")

#USIAMO R2_adj PER FAVORE
fig, ax = plt.subplots()
ax.errorbar(R2s_adj['Number of species'], R2s_adj['Average R2'], yerr=R2s_adj['R2 errors'], fmt='o')
ax.set_xlabel('Sample complexity')
ax.set_xticks([2,3,4,5])
ax.axhline(y=1, color='r', linestyle='--')
ax.set_ylabel('R2')

# Show the plot
plt.show()