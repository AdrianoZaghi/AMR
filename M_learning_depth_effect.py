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

#Lunghezza media delle reads simulate e reali, comunque un fattore di normalizzazione che dovrebbe valere per entrambi
#potrebbe essere diverso tra simulate e reali
average_read_length = 270
images_dir = "M_learning_images/"
campioni = pd.read_csv("Database_campioni2.csv", index_col=0)

Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
Ada =  AdaBoostRegressor(estimator=Base_estimator_for_adaboost,n_estimators=100, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1, 0.5, 1, 2]}
#SELECTING THE MODEL TO USE

model = Ada
params = params_Ada
Transform_the_targets = None # "Log_shift" "Logistic" None
Images = True
samples_dir = "samples_2/"
species_per_genome = 5
save_dir = "Sembra_che_funzioni_depth_"+ str(species_per_genome) +"_strains/"
if os.path.exists(save_dir) == False:
	os.makedirs(save_dir)

t = campioni[campioni["Size"] == 0.1]

t_train = t[t["Aureus_R_N_strains"] == species_per_genome]
t_train = t_train[t_train["Real/Simul"] == "Simul"]
Train = [samples_dir + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_train.index]
t_train = t_train[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
t_train = t_train.div(t_train.sum(axis=1), axis=0) #non si sa mai
t_train = t_train[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
Train_target = dict(zip(t_train.index, t_train.values))

t_test = t[t["Aureus_R_N_strains"] == 0]
t_test = t_test[t_test["Real/Simul"] == "Real"]
Test = [samples_dir + x + "/bwt_out/bwt.allele_mapping_data.json" for x in t_test.index]
t_test = t_test[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]
t_test = t_test.div(t_test.sum(axis=1), axis=0) #non si sa mai
t_test = t_test[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
Test_target = dict(zip(t_test.index, t_test.values))

MAPQ_MIN = 0
MAPQ_MAX = 1000
JUST_COMPLETE = True
average_read_length = 270
depths = [0.001, 0.003, 0.007, 0.01, 0.05, 0.1, 0.5, 1]
Results = {}
#%%
for sub in depths:
	print("Depth: "+str(sub)+" ", end = '')
	R2s_for_a_dept = []
	for i in range(0,6):
		print(".", end = '')
		#Test_data, Train_data = im.read_bwt_2_aro(Test, sub), im.read_bwt_2_aro(Train, sub)
		Train_data = im.read_many_alignements_json(Train,
											 just_completely_mapped=JUST_COMPLETE,
											 mapq_average_min=MAPQ_MIN,
											 mapq_average_max=MAPQ_MAX,
											 read_length=average_read_length,
											 subsample=sub).T
		Test_data = im.read_many_alignements_json(Test,
											just_completely_mapped=JUST_COMPLETE,
											mapq_average_min=MAPQ_MIN,
											mapq_average_max=MAPQ_MAX,
											read_length=average_read_length,
											subsample=sub).T											 	
		missing_in_Test, missing_in_Train = [c for c in Train_data.columns if c not in Test_data.columns], [c for c in Test_data.columns if c not in Train_data.columns]
		Test_data, Train_data = Test_data.drop(columns=missing_in_Train), Train_data.drop(columns=missing_in_Test)
		Test_data, Train_data = Test_data.sort_index(axis=1), Train_data.sort_index(axis=1)
		print("Dati letti.\nSto usando "+str(len(Train_data.columns))+" colonne nel fit")
		grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=5, verbose=0, scoring="r2")
		#Fitting the model
		grid.fit(Train_data, list(Train_target.values()))
		#Retreaving best estimator
		best = grid.best_estimator_
		#Making predicitons and evaluate performances
		Test_predictions = grid.predict(Test_data)
		Test_predictions = dict(zip(Test_data.index, Test_predictions))
		R2 = r2_score(list(Test_target.values()), list(Test_predictions.values())) #il R2 ridotto non mi sembra adatto, non ho un numero preciso di parametri che fitto
		R2s_for_a_dept.append(R2)
	print(R2s_for_a_dept)
	Results["Subsample " + str(sub)] = {"R2" : R2s_for_a_dept, "Best_estimator" : best}
	if Images:
		#Plotting results["Best_estimator"]
		im.plot_results(list(Test_target.values()), list(Test_predictions.values()), save_dir+"/Depth_subsample_" + str(sub) + ".png", R2 = R2)
		#Feature importance
		#F_imp = im.feature_importance(Train_data, best, "M_learning_images/"+WHICH+"_feature_importance.png")
		#plt.clf()
		#Results["Subsample " + str(sub)]["Feature_importance"] = F_imp

R2s = [[float(k.split("Subsample ")[1]), v["R2"]] for k,v in Results.items()]
R2s = pd.DataFrame(R2s, columns = ["Subsample", "R2"])
R2s["Subsample"] = R2s["Subsample"]*0.1
R2s["Average R2"] = R2s["R2"].apply(np.mean)
R2s["R2 errors"] = R2s["R2"].apply(np.std)
R2s.to_csv("M_learning_images/R2s_depth_effects_patric.csv")
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
