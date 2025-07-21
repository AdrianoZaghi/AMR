#%%
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
#import seaborn as sns
import random
import os
import dask.dataframe as dd
import json

#Lunghezza media delle reads simulate e reali, comunque un fattore di normalizzazione che dovrebbe valere per entrambi
#potrebbe essere diverso tra simulate e reali

Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
params = {"loss":["linear", "square", "exponential"], "learning_rate":[1]}

species_per_genome = 1
PREVALENCE_FILTER = 0

campioni = pd.read_csv("Database_campioni2.csv", index_col=0)
samples_dir = "samples_2/"
t = campioni[campioni["Size"] == 0.1]
t_train = t[t["Aureus_R_N_strains"] == species_per_genome]
t_train = t_train[t_train["Real/Simul"] == "Simul"]
average_read_length_train = 270
assert all([x==average_read_length_train for x in t_train["reads_length"]])
t_train = t_train[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
t_train = t_train.div(t_train.sum(axis=1), axis=0) #non si sa mai
t_train = t_train[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
#t_train = t_train[["Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]].sum(axis=1)

Train = [samples_dir + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_train.index]
Train_target = dict(zip(t_train.index, t_train.values))

save_dir = "DIOO_Filtered_depth_"+str(species_per_genome) +"_strains_prev"+str(PREVALENCE_FILTER)+"_train_col/"
print(save_dir)
print("Train data: "+str(len(Train)))
if os.path.exists(save_dir) == False:
	os.makedirs(save_dir)

campioni = pd.read_csv("Database_campioni2.csv", index_col=0)
samples_dir = "samples_2/"
t_test = campioni[campioni["Size"] == 0.1]
t_test = t_test[t_test["Real/Simul"] == "Real"]
#elimina i campioni veri che non usiamo
campioni_veri_che_non_usiamo = ["DCWGS"+i+"_" for i in ["1", "2", "3", "4", "5", "6", "7", "AR", "PA", "EC", "SAS", "SAR"]]
for i in campioni_veri_che_non_usiamo:
	t_test = t_test[~t_test.index.str.contains(i)]
average_read_length_test = 250
assert all([x==average_read_length_test for x in t_test["reads_length"]])
Test = [samples_dir + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_test.index]
t_test = t_test[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
t_test = t_test.div(t_test.sum(axis=1), axis=0) #non si sa mai
t_test = t_test[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
#t_test = t_test[["Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]].sum(axis=1)
Test_target = t_test.to_dict()
five_species_samples = ["DCWGS41", "DCWGS42", "DCWGS43", "DCWGS40"]
Test_target_five_species = {x: Test_target[x] for x in Test_target if any([y in x for y in five_species_samples])}

import matplotlib.pyplot as plt
#%%
plt.figure(figsize=(8, 5))
plt.hist(t_test.values, bins=20, color='skyblue', edgecolor='black')
plt.xlabel('Valori di t_test')
plt.ylabel('Frequenza')
plt.title('Istogramma dei valori di t_test')
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(t_train.values, bins=20, color='skyblue', edgecolor='black')
plt.xlabel('Valori di t_train')
plt.ylabel('Frequenza')
plt.ylim(0, 200)
plt.title('Istogramma dei valori di t_train')
plt.tight_layout()
plt.show()
assert 1==2

#%%
col_numbers_collect = []
MAPQ_MIN = 0
MAPQ_MAX = 1000
JUST_COMPLETE = True
Results_collect = []
for FILTRO_ABBONDANZA in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
	print("##########################################")
	print("FILTRO ABBONDANZA: "+str(FILTRO_ABBONDANZA))
	depths = [0.01, 0.05, 0.1, 0.15, 0.30, 0.5, 0.70, 0.85, 0.95, 1]
	Results = {}
	col_numbers = {}
	for sub in depths:
		print("Depth: "+str(sub)+" "+" Filtro_abbondanza: "+str(FILTRO_ABBONDANZA), end = '')
		Train_data = im.read_many_alignements_json_1(Train,
												 just_completely_mapped=JUST_COMPLETE,
												 mapq_average_min=MAPQ_MIN,
												 mapq_average_max=MAPQ_MAX,
												 read_length=average_read_length_train,
												 subsample=sub,
												 filtro_abbondanza=FILTRO_ABBONDANZA).T
		Test_data = im.read_many_alignements_json_1(Test,
												just_completely_mapped=JUST_COMPLETE,
												mapq_average_min=MAPQ_MIN,
												mapq_average_max=MAPQ_MAX,
												read_length=average_read_length_test,
												subsample=sub,
												filtro_abbondanza=FILTRO_ABBONDANZA).T
		C_1 = set(Train_data.columns).intersection(set(Test_data.columns))
		print(len(C_1), "Colonne in comune before prevalence filtering")
		col_prevalence = (Train_data > 0).sum(axis=0)/ len(Train_data)
		Train_data = Train_data.loc[:, col_prevalence > 0.2]
		C_2 = set(Train_data.columns).intersection(set(Test_data.columns))
		print(len(C_2), "Colonne in comune after prevalence filtering")
		#Test_data, Train_data = Test_data.sort_index(axis=1), Train_data.sort_index(axis=1)
		#DI NUOVO, TENGO SOLO LE COLONNE CHE SONO COMUNI AI DUE DATAFRAME
		#COLUMNS = list(set(Train_data.columns).intersection(set(Test_data.columns)))
		COLUMNS = list(Train_data.columns)
		if len(COLUMNS) == 0:
			print("Nessuna colonna in comune, salto il subsample "+str(sub))
			Results["Subsample " + str(sub)] = {"R2" : [0,0], "R2_five_species": [0,0]}
			continue
		Test_data =  Test_data.reindex(columns = COLUMNS, fill_value=0).sort_index(axis=1)
		Train_data = Train_data.reindex(columns = COLUMNS, fill_value=0).sort_index(axis=1)
		print("Dati letti.\nSto usando "+str(len(Train_data.columns))+" colonne nel fit")
		R2s_for_a_dept = []
		R2s_for_a_dept_five = []
		for i in range(0,2):
			print(".", end = '')
			#Test_data, Train_data = im.read_bwt_2_aro(Test, sub), im.read_bwt_2_aro(Train, sub)
			model =  AdaBoostRegressor(estimator=Base_estimator_for_adaboost,n_estimators=100, random_state=i)
			grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=5, verbose=0, scoring="r2")
			#Fitting the model
			grid.fit(Train_data, list(Train_target.values()))
			#Retreaving best estimator
			best = grid.best_estimator_
			#Making predicitons and evaluate performances
			Test_predictions = grid.predict(Test_data)
			Test_predictions = pd.Series(Test_predictions, index=Test_data.index)
			R2 = r2_score(list(Test_target.values()), list(Test_predictions)) #il R2 ridotto non mi sembra adatto, non ho un numero preciso di parametri che fitto
			R2_five_species = r2_score(list(Test_target_five_species.values()),
									 [Test_predictions[x] for x in Test_target_five_species if x in Test_predictions.index])
			R2s_for_a_dept.append(R2)
			R2s_for_a_dept_five.append(R2_five_species)
		Results["Subsample " + str(sub)] = {"R2" : R2s_for_a_dept, "R2_five_species": R2s_for_a_dept_five}
		json.dump(Results, open(save_dir+"/Temporary_results_depth_effects.json", "w"), indent=4)	
		predictions = pd.DataFrame({
			"ID": Test_data.index,
			"Test_target": list(Test_target.values()),
			"Test_predictions": list(Test_predictions)
			}).to_csv(save_dir+"/Depth_subsample_" + str(sub) + ".csv")
		im.plot_results(list(Test_target.values()), Test_predictions.tolist(), save_dir+"/Depth_subsample_" + str(sub) + ".png", R2 = R2)
	col_numbers["Depth " + str(sub)] = {"Pre_prevalence_filtering": len(C_1),
						"Post_prevalence_filtering": len(C_2),
						"FILTRO_ABBONDANZA": FILTRO_ABBONDANZA,
						"Subsample": depths}
	Results_collect.append(Results)
	col_numbers_collect.append(col_numbers)
json.dump(Results_collect, open(save_dir+"/Results_depth_effects.json", "w"), indent=4)
json.dump(col_numbers_collect, open(save_dir+"/col_numbers_depth_effects.json", "w"), indent=4)
assert 1==2

R2s = [[float(k.split("Subsample ")[1]), v["R2"]] for k,v in Results.items()]
R2s = pd.DataFrame(R2s, columns = ["Subsample", "R2"])
R2s["Subsample"] = R2s["Subsample"]*0.1
R2s["Average R2"] = R2s["R2"].apply(np.mean)
R2s["R2 errors"] = R2s["R2"].apply(np.std)
R2s.to_csv(save_dir+"/R2s_depth_effects.csv")
fig, ax = plt.subplots(figsize=(10, 6))
# Plot the data with error bars
ax.errorbar(R2s['Subsample'], R2s['Average R2'], yerr=R2s['R2 errors'], fmt='o')
# Set the x-axis scale to logarithmic
ax.set_xscale('log')
# Add a horizontal line at y=1
ax.axhline(y=1, color='r', linestyle='--')
# Set y-axis limits
ax.set_ylim(0, 1.1)
# Set the title and labels
ax.set_title("Effect of Sample Depth")
ax.set_xlabel('Subsample Size (Gbp)')
ax.set_ylabel('R2')

# Show the plot
#plt.show()
plt.savefig(save_dir+"/R2s_depth_effects.png")
# %%
