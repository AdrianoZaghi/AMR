#%%
import json
import pandas as pd
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
import dask.dataframe as dd
import sys

def rmse(y_true, y_pred):
    return np.sqrt(mean_squared_error(y_true, y_pred))
rmse_scorer = make_scorer(rmse, greater_is_better=False)

def plot_x_y(x,y):
	plt.figure(figsize=(10, 6))
	plt.plot(x, y, marker='o')
	plt.xlabel('Columns')
	plt.ylabel('Sum of values')
	plt.yscale("log")
	plt.tight_layout()
	plt.show()
	plt.close()

Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
#params = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1, 0.5, 1, 2]}#SELECTING THE MODEL TO USE
params = {"loss":["linear", "square", "exponential"], "learning_rate":[1]}#SELECTING THE MODEL TO USE

print("Si parteeeee")

JUST_COMPLETE = True
fold_number = 5
MAPQ_MIN = 0
MAPQ_MAX = 1000
FILTRO_ABBONDANZA = 6
FILTRO_PREVALENZA = 0.1
SUBSAMPLE = 0.2

if len(sys.argv) > 2:
	NUMBER_OF_STRAINS = int(sys.argv[1])
else:
	NUMBER_OF_STRAINS = 1

images_dir = "New/Complexity_sub0.02__" + str(NUMBER_OF_STRAINS) + "_strains/"
print(images_dir)
if not os.path.exists(images_dir):
	os.makedirs(images_dir)

SCORING = "r2" # "r2" or "mean_squared_error"
#%%
campioni = pd.read_csv("Database_campioni2.csv", index_col=0)
t = campioni[campioni["Size"] == 0.1]
t_train = t[t["Real/Simul"] == "Simul"]
t_train = t_train[t_train["Aureus_S_N_strains"] == NUMBER_OF_STRAINS]
average_read_length_train = 270
Train_paths = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_train.index]
Train_data = im.read_many_alignements_json_1(Train_paths,
										   just_completely_mapped=JUST_COMPLETE,
										   mapq_average_min=MAPQ_MIN,
										   mapq_average_max=MAPQ_MAX,
										   read_length=average_read_length_train,
										   subsample=SUBSAMPLE,
										   filtro_abbondanza=FILTRO_ABBONDANZA)
print("Ho letto le train features")
Train_target = dict((t_train["Aureus_R"] + t_train["Coli"] + t_train["Pseudomonas"]) /
				  (t_train["Arcobacter"] + t_train["Aureus_S"] + t_train["Aureus_R"] + t_train["Coli"] + t_train["Pseudomonas"]))

print("Ho fatto il trian set, features e target")

Train_data = Train_data.T
col_prevalence = (Train_data > 0).sum(axis=0)/ len(Train_data)
Train_data = Train_data.loc[:, col_prevalence > FILTRO_PREVALENZA]

#%%
campioni_per_test = pd.read_csv("Database_campioni2.csv", index_col=0)
t_test = campioni_per_test[campioni_per_test["Size"] == 0.1]
t_test = t_test[t_test["Real/Simul"] == "Real"]
campioni_veri_che_non_usiamo = ["DCWGS"+i+"_" for i in ["1", "2", "3", "4", "5", "6", "7", "AR", "PA", "EC", "SAS", "SAR"]]
for i in campioni_veri_che_non_usiamo:
	t_test = t_test[~t_test.index.str.contains(i)]
average_read_length_test = 250
#assert all([x==average_read_length_test for x in t_test["reads_length"]])
t_test = t_test[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
Test_paths = ["samples_2/" + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_test.index]
Test_data = im.read_many_alignements_json_1(Test_paths,
											just_completely_mapped=JUST_COMPLETE,
											mapq_average_min=MAPQ_MIN,
											read_length=average_read_length_test,
											mapq_average_max=MAPQ_MAX,
											subsample=SUBSAMPLE,
											filtro_abbondanza=FILTRO_ABBONDANZA)
Test_data = Test_data.T
Test_target = dict(t_test[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1) / t_test[["Arcobacter", "Aureus_S", "Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1))	
select_target_complexity = t_test != 0
select_target_complexity = select_target_complexity.astype(int)
print("Ho fatto il test set, features e target")
COLUMNS = list(set(Train_data.columns).intersection(set(Test_data.columns)))
Test_data_all_col = Test_data.reindex(columns = COLUMNS, fill_value=0).sort_index(axis=1)
Train_data_all_col = Train_data.reindex(columns = COLUMNS, fill_value=0).sort_index(axis=1)

#%%
Results = {}
for iii in range (0, 5): #5 replicas
	print("Replica N " + str(iii))
	model =  AdaBoostRegressor(estimator=Base_estimator_for_adaboost,n_estimators=100, random_state=iii)
	grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=fold_number, verbose=0, scoring="r2")
	print("Start fit")
	grid.fit(Train_data_all_col, list(Train_target.values()))
	train_average_score = np.mean([np.mean(v) for k,v in grid.cv_results_.items() if "split" in k])
	print("Train_set_av_score   "+str(train_average_score)+"   .")
	best = grid.best_estimator_
	Test_scores_for_a_model = {}
	Test_predictions_for_a_model = {}
	for i in range(2,6):
		Test_data_all_col_ = Test_data_all_col[select_target_complexity.sum(axis=1) == i]
		Test_target_ = {k : v for k,v in Test_target.items() if k in Test_data_all_col_.index}
		if True:#adding colors for each bacteria combination
			ttc = pd.DataFrame(select_target_complexity[select_target_complexity.sum(axis=1) == i])
			ttc["combination"] = ttc.apply(lambda row: '_'.join(row.astype(str)), axis=1)
			unique_combinations = ttc['combination'].unique()
			colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_combinations)))
			color_map = dict(zip(unique_combinations, colors))
			colors_mapped = ttc['combination'].map(color_map)
		Test_predictions = grid.predict(Test_data_all_col_).tolist()
		if SCORING == "r2":
			Score = round(r2_score(list(Test_target_.values()), list(Test_predictions)), 3)
		elif SCORING == "mean_squared_error":
			Score = rmse(list(Test_target_.values()), list(Test_predictions))
		print("Score   "+str(Score)+"   .")
		Test_predictions_for_a_model["Number of species: " + str(i)] = {Test_data_all_col_.index[j]:
																  			{"Test_target": Test_target_[Test_data_all_col_.index[j]],
					  														"Test_predictions": Test_predictions[j],
																			"Color": colors_mapped[Test_data_all_col_.index[j]].tolist()} 
																		for j in range(len(Test_data_all_col_.index))}
		#Test_predictions_for_a_model.update(dict(zip(Test_data_all_col_.index, Test_predictions)))
		Test_scores_for_a_model["Number of species: " + str(i)] = Score
	Results["Model replica N " + str(iii)] = {"Test_scores" : Test_scores_for_a_model, "Test_predictions": Test_predictions_for_a_model, "Train_scores" : float(train_average_score)}
#Salva results in json
with open(images_dir+"Results.json", 'w') as f:
	json.dump(Results, f)

#%%
if False:#Serve per fare gli scatter plot direttamente da i risultati una volta che lo script è finito
	with open(images_dir+"Results.json", 'r') as f:
		Results = json.load(f)
		M = "Model replica N 2"
		for num_spec in Results[M]["Test_predictions"].keys():
			print(num_spec)
			test_target = [x["Test_target"] for x in Results[M]["Test_predictions"][num_spec].values()]
			test_predictions = [x["Test_predictions"] for x in Results[M]["Test_predictions"][num_spec].values()]
			colors = [x["Color"] for x in Results[M]["Test_predictions"][num_spec].values()]
			im.plot_results_2(test_target,
						test_predictions,
						colors = colors)
#%%
Test_scores = pd.DataFrame({k : v["Test_scores"] for k,v in Results.items()}).T
Test_scores.to_csv(images_dir+"Test_scores.csv")
Average_score_for_specie = Test_scores.mean(axis=0).values
Standars_dev_for_specie = Test_scores.std(axis=0).values
Test_predictions = pd.DataFrame({k : v["Test_predictions"] for k,v in Results.items()}).T
Test_predictions.to_csv(images_dir+"Test_predictions.csv")
Train_scores = pd.DataFrame([[k ,v["Train_scores"]] for k,v in Results.items()], columns = ["Model", "Train_scores"])
Train_scores.to_csv(images_dir+"Train_scores.csv")

fig, ax = plt.subplots()
ax.errorbar(range(2,6), Average_score_for_specie, yerr=Standars_dev_for_specie, fmt='o')
ax.set_xlabel('Sample complexity')
ax.set_xticks([2,3,4,5])
if SCORING == "r2":
	ax.axhline(y=1, color='r', linestyle='--')
elif SCORING == "mean_squared_error":
	ax.axhline(y=0, color='r', linestyle='--')
ax.set_ylabel(SCORING)
ax.set_ylim(0, 1.1)
ax.set_title("Effect of sample complexity_"+SCORING)
#ax.set_ylim(0.5,1.1)

# Show the plot
plt.savefig(images_dir+"/Complexity_effect_synth.png")
# %%
