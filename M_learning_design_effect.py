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

print("Si parteeeee")
#Lunghezza media delle reads simulate e reali, comunque un fattore di normalizzazione che dovrebbe valere per entrambi
#potrebbe essere diverso tra simulate e reali

campioni = pd.read_csv("Database_campioni2.csv", index_col=0)

def rmse(y_true, y_pred):
    return np.sqrt(mean_squared_error(y_true, y_pred))
rmse_scorer = make_scorer(rmse, greater_is_better=False)

Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
Ada =  AdaBoostRegressor(estimator=Base_estimator_for_adaboost,n_estimators=100)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1, 0.5, 1, 2]}#SELECTING THE MODEL TO USE

model = Ada
params = params_Ada

WHICH = "De_Caesar_simul"#"Simul"#
JUST_COMPLETE = True
NUMBER_OF_STRAINS = 5
MAPQ_MIN = -1
MAPQ_MAX = 1000
SCORING = "r2" # "r2" or "mean_squared_error"
images_dir = f"Sembra_che_funzionii_{WHICH}_strains_{NUMBER_OF_STRAINS}_scoring_{SCORING}_mapping_{MAPQ_MIN}_{JUST_COMPLETE}/"
print(images_dir)
if not os.path.exists(images_dir):
	os.makedirs(images_dir)
print(WHICH)
#Deciding which samples to peack from the database and defining the train and test sets
t = campioni[campioni["Size"] == 0.1]
#Selecting train samples
#t_train = t[t["Real/Simul"] == "Simul"]
t_train = t[t["Real/Simul"] == WHICH]
t_train = t_train[t_train["Aureus_S_N_strains"] == NUMBER_OF_STRAINS]
Train_ = ["samples_2/" + x for x in t_train.index]
#Train_data = im.read_bwt_2_aro(Train, subsample_counts, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, av_read_length=255)
Train = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_train.index]
Train_data = im.read_many_alignements_json(Train, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, mapq_average_max=MAPQ_MAX, read_length=255)
Train_data = Train_data.T
print("Ho letto le features")
Train_target = dict((t_train["Aureus_R"] + t_train["Coli"] + t_train["Pseudomonas"]) / (t_train["Arcobacter"] + t_train["Aureus_S"] + t_train["Aureus_R"] + t_train["Coli"] + t_train["Pseudomonas"]))
print("Ho letto tuttooooooo")
#Preparing the test samples for selection
t_test = t[t["Real/Simul"] == "Real"]
#t_test = t[t["Real/Simul"] == "Real_new_trim"]
t_test = t_test[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
select_target_complexity = t_test != 0
select_target_complexity = select_target_complexity.astype(int)

print("Entriamo nel looooop")
Results = {}
for i in range(2,6):
	t_test_ = t_test[select_target_complexity.sum(axis=1) == i]
	Test = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_test_.index]
	Test_data = im.read_many_alignements_json(Test, just_completely_mapped=JUST_COMPLETE, mapq_average_min=MAPQ_MIN, read_length=255)
	Test_data = Test_data.T
	if True:#adding colors for each bacteria combination
		ttc = pd.DataFrame(select_target_complexity[select_target_complexity.sum(axis=1) == i])
		ttc["combination"] = ttc.apply(lambda row: '_'.join(row.astype(str)), axis=1)
		unique_combinations = ttc['combination'].unique()
		colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_combinations)))
		color_map = dict(zip(unique_combinations, colors))
		colors_mapped = ttc['combination'].map(color_map)
		#print(color_map)
	print("Sto facendo il fit per "+str(i)+" specie")
	#Sto prendendo il target dai metadata
	Test_target = dict(t_test_[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1) / t_test_[["Arcobacter", "Aureus_S", "Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1))
	#Dealing with columns absent in one of the two datasets
	if True: #Drop columns which are not common to both datasets	
		missing_in_Train = [c for c in Test_data.columns if c not in Train_data.columns]
		missing_in_Test = [c for c in Train_data.columns if c not in Test_data.columns]
		Test_data_clean = Test_data.drop(columns=missing_in_Train).sort_index(axis=1)
		Train_data_clean = Train_data.drop(columns=missing_in_Test).sort_index(axis=1)
	else: #Add columns which are missing in one of the two datasets filling with 0
		all_columns = list(set(Test_data.columns).union(set(Train_data.columns)))
		Test_data_clean = Test_data.reindex(columns = all_columns, fill_value=0).sort_index(axis=1)
		Train_data_clean = Train_data.reindex(columns = all_columns, fill_value=0).sort_index(axis=1)
	print("Sto usando "+str(len(Train_data_clean.columns))+" colonne nel fit")
	Train_scores_for_species_number = []
	Test_scores_for_species_number = []
	#Per ogni numero di specie, faccio 5 fit diversi per ottenere dei valori di R2 piu' robusti
	for iii in range(0,5):
		if SCORING == "r2":
			grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=10, verbose=0, scoring="r2")
		elif SCORING == "mean_squared_error":
			grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=10, verbose=0, scoring=rmse_scorer)
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
		if SCORING == "r2":
			Score = r2_score(list(Test_target.values()), list(Test_predictions.values())) #il R2 ridotto non mi sembra adatto, non ho un numero preciso di parametri che fitto
		elif SCORING == "mean_squared_error":
			Score = rmse(list(Test_target.values()), list(Test_predictions.values()))
		print("Score   "+str(Score)+"   .")
		im.plot_results_2(list(Test_target.values()), list(Test_predictions.values()), colors=colors_mapped, legend = None, filename = images_dir+WHICH+"_complexity_" + str(i)+ "_" + str(iii) + ".png", Score = Score, Score_name=SCORING, y_low_lim=0, y_up_lim=1, x_low_lim=0, x_up_lim=1)
		Test_scores_for_species_number.append(Score)
		Train_scores_for_species_number.append(train_average_score)
	Results["Number of species: " + str(i)] = {"Test_scores" : Test_scores_for_species_number, "Train_scores" : Train_scores_for_species_number}

Scores = [[int(k.split("Number of species: ")[1]), v["Train_scores"], v["Test_scores"]] for k,v in Results.items()]
Scores = pd.DataFrame(Scores, columns = ["Number of species", "Train_scores", "Test_scores"])
Scores["Average_test_score"] = Scores["Test_scores"].apply(np.mean)
Scores["test_std"] = Scores["Test_scores"].apply(np.std)
Scores["Average_train_score"] = Scores["Train_scores"].apply(np.mean)
Scores["train_std"] = Scores["Train_scores"].apply(np.std)
Scores.to_csv(images_dir+"_"+SCORING+"_complexity_effects.csv")

fig, ax = plt.subplots()
ax.errorbar(Scores['Number of species'], Scores['Average_test_score'], yerr=Scores['test_std'], fmt='o')
ax.set_xlabel('Sample complexity')
ax.set_xticks([2,3,4,5])
if SCORING == "r2":
	ax.axhline(y=1, color='r', linestyle='--')
elif SCORING == "mean_squared_error":
	ax.axhline(y=0, color='r', linestyle='--')
ax.set_ylabel(SCORING)
ax.set_title("Effect of sample complexity_"+SCORING)
#ax.set_ylim(0.5,1.1)


# Show the plot
plt.savefig(images_dir+"/Complexity_effect_synth.png")