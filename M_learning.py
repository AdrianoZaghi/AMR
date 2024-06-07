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
#Declaring the models and initializing their grid of hyperparameters
El = ElasticNet(max_iter=100000, tol=0.0005)
params_El = {"alpha":[0.2, 0.5, 0.75, 0.85,1], "l1_ratio":[0,0.5,1]}
params_El = {"alpha":[0.2, 0.5, 0.75, 0.85,1, 2], "l1_ratio":[0, 0.3, 0.5, 0.7 ,1]}
#Con questa non serve levare i campinoi estremi che anzi hanno un'effetto stabilizzante
#È evidente la presenza di un bias

Forest = RandomForestRegressor(random_state = 0, n_estimators=1000, n_jobs=-1, oob_score=True)
params_Forest = {"criterion":["squared_error", "absolute_error", "friedman_mse"]}
#Gli ensamble non lineari, tendono a classificare gli 0, gli uno e qulli intermedi. Performano male su quelli intermedi

Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=5)
Ada =  AdaBoostRegressor(base_estimator=Base_estimator_for_adaboost,n_estimators=100, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1]}
#SELECTING THE MODEL TO USE

model = Ada
params = params_Ada
Transform_the_targets = None # "Log_shift" "Logistic" None
subsample_counts = 1
Eliminate_extreme_samples = False
Subsample_train_samples = 0
WHICH = "Synth_on_Real_assembled_genomes_0.1" #"Synth_on_Real_patric_genomes" #"Campioni_per_Ettore" #"Just_1_specie" #"Synth_on_Real_assembled_genomes_0.1"
MANUAL_FEATURE_SELECTION = False
Plotting_train_target = False
Zeros_analysis = False
Add_second_degree_features = False
Compare_test_targets = False
regression = False
READ_TARGET_FROM =  "Database_campioni" # "Tax_Ettore" #

#Deciding which samples to peack from the database and defining the train and test sets
print(WHICH)
if(WHICH == "Synth_on_Synth_assembled_genomes"):
	t = campioni[campioni["Size"] == 0.01]
	Simul = ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Synth_on_Synth_assembled_genomes")
elif(WHICH == "Synth_on_Synth_patric_genomes"):
	t = campioni[campioni["Size"] == 0.01]
	Simul = ["samples_2/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Synth_on_Synth_patric_genomes")
elif(WHICH == "Synth_on_Real_assembled_genomes_0.1"):
	t = campioni[campioni["Size"] == 0.1]
	Train = ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
	Test = ["samples_2/" + x for x in t[t["Real/Simul"] == "Real"].index]
	print("Synth_on_Real_assembled_genomes_0.1")
elif(WHICH == "Synth_on_Real_assembled_genomes"):
	t = campioni[campioni["Size"] == 0.01]
	Train = ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
	if False:
		Train = random.sample(Train, 200)
	Test = ["samples_2/" + x for x in t[t["Real/Simul"] == "Real"].index]
	print("Synth_on_Real_assembled_genomes")
elif(WHICH == "Synth_on_Real_patric_genomes"):
	t = campioni[campioni["Size"] == 0.01]
	Train = ["samples_2/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Test = ["samples_2/" + x for x in t[t["Real/Simul"] == "Real"].index]
	print("Synth_on_Real_patric_genomes")
elif(WHICH == "Just_1_specie"):
	t = campioni[campioni["Size"] == 0.1]
	t = t[t["Arcobacter"] > 0.5] # ettore non trova arcobacter
	t = t[t["Aureus_S"] > 0.5] # ettore lo sovrastima
	t = t[t["Aureus_R"] > 0.5] # ettore da dei conteggi ad aureus_s che secondo De Cesare sono di aureus_r
	t = t[t["Pseudomonas"] > 0.5] # ettore da dei conteggi aureus_s che secondo De Cesare sono di pseudomonas
	t = t[t["Coli"] > 0.5] # in accordo
	Train = ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
	Test = ["samples_2/" + x for x in t[t["Real/Simul"] == "Real"].index]
	print("Just_1_specie")
elif(WHICH == "Aureus_S_e__R_Synth_on_synth"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[t["Coli"] == 0]
	t = t[t["Arcobacter"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = ["samples_2/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Train, Test = train_test_split(t, test_size=0.1, random_state=0)
	print("Synth_on_Real_patric_genomes")
elif(WHICH == "Campioni_per_Ettore"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[t["Real/Simul"] == "Simul"]
	t = t[t.max(axis=1) == 100]


if Subsample_train_samples > 0 and Subsample_train_samples < len(Train):
		Train = random.sample(Train, Subsample_train_samples)
print("Quanti trains: " + str(len(Train)))
print("Quanti tests: " + str(len(Test)))

#################READING DATA
#Features (this funciton is also scaling the features to make them relatable between samples)
Train_data = im.read_bwt_2_aro(Train, subsample_counts)
Test_data = im.read_bwt_2_aro(Test, subsample_counts)

#Dealing with columns absent in one of the two datasets 
missing_in_Train = [c for c in Test_data.columns if c not in Train_data.columns]
missing_in_Test = [c for c in Train_data.columns if c not in Test_data.columns]
Test_data = Test_data.drop(columns=missing_in_Train)
Train_data = Train_data.drop(columns=missing_in_Test)

Test_data = Test_data.sort_index(axis=1)
Train_data = Train_data.sort_index(axis=1)
print("Sto usando "+str(len(Train_data.columns))+" colonne nel fit")

#Adding second degree features
if Add_second_degree_features:
	for col in Test_data.columns:
		if col in ['3000178', '3000621', '3000676', '3000565', '3003905', '3000839', '3000815', '3001209', '3005046', "3007041"]:
			Test_data[col + '_squared'] = Test_data[col].apply(np.square)
			Train_data[col + '_squared'] = Train_data[col].apply(np.square)
			Test_data[col + '_cube'] = Test_data[col].apply(lambda x: x**3)
			Train_data[col + '_cube'] = Train_data[col].apply(lambda x: x**3)
#Manual feature selection
			
if MANUAL_FEATURE_SELECTION:
	selected = ["3005046", "3000621", "3004661", "3005097", "3001209", "3000347"]
	Train_data = Train_data[selected]
	Test_data = Test_data[selected]
#	de_selected = ["3005046", "3000621", "3004661", "3005097", "3001209", "3000347"]
#	Train_data = Train_data.drop(columns=de_selected)
#	Test_data = Test_data.drop(columns=de_selected)

All_zero_data = pd.DataFrame(0, index=range(10), columns=Train_data.columns)

######################### GUARDANDO LE HEATMAP CON LE ABBONDANZE DEI CMAPIONI DI ARCOBACTER
# 	WHICH == "Just_1_specie"
#TOT = pd.concat([Train_data, Test_data])
#im.heatmap_abundance(TOT, Normalize = False, figsize = (15,20))

#reading targets
Train_target= im.read_metadata(Train)["Abundances_r"]
if READ_TARGET_FROM == "Tax_Ettore":
	Test_target = im.read_abundances_from_tax_ettore(Test)
	Test_target = {k:round(v["Aureus_R"]+v["Pseudomonas"]+v["Coli"],5) for k,v in Test_target.iterrows()}
	Test_data = Test_data[Test_data.index.isin(Test_target.keys())]
elif READ_TARGET_FROM == "Database_campioni":
	Test_target = im.read_metadata(Test)["Abundances_r"]

if Eliminate_extreme_samples:
	Ltr, Lte = len(Train_target), len(Test_target)
	Train_target = {k:v for k,v in Train_target.items() if v != 1 and v != 0}
	Test_target = {k:v for k,v in Test_target.items() if v != 1 and v != 0}
	Train_data = Train_data.loc[Train_target.keys()]
	Test_data = Test_data.loc[Test_target.keys()]
	print("Eliminate "+str(Ltr - len(Train_target))+" extreme_samples from train")
	print("Eliminate "+str(Lte - len(Test_target))+" extreme_samples from test")

if Compare_test_targets:
	quanti = 30
	Test_target_ettore = im.read_abundances_from_tax_ettore(Test).head(quanti)
	Test_target_decesare = [ind for ind, row in campioni.iterrows() if "samples_2/"+ind in Test_target_ettore.index]
	Test_target_decesare = campioni.loc[Test_target_decesare]
	Test_target_decesare = Test_target_decesare.sort_index().reset_index(drop=True)
	Test_target_ettore = Test_target_ettore.sort_index().reset_index(drop=True)
	x = np.arange(Test_target_decesare.shape[0])
	for c in Test_target_ettore.columns:
		print(c)
		to_plot=pd.DataFrame([Test_target_decesare[c] , Test_target_ettore[c]])
		to_plot.index =["De Cesare", "Ettore"]
		fig, ax = plt.subplots(layout='constrained')
		i = 0
		for attribute, measurement in to_plot.iterrows():
			rects = ax.bar(x + 0.25*i, measurement, 0.25, label=attribute)
			i += 1
		ax.legend(loc='upper left', ncols=2)
		ax.set_title(c)
		plt.show()

if Transform_the_targets == "Log_shift":	
	shift = 0.01
	Train_target_trans = im.log_shift(Train_target, shift = shift)
	Test_target_trans = im.log_shift(Test_target, shift = shift)
	print("Log_shift transformation del target")
elif Transform_the_targets == "Logistic":
	border = 0.001
	Target_transformation = "Logistic"
	Train_target_trans = im.logistic(Train_target, min = border, max = 1-border)
	Test_target_trans = im.logistic(Test_target, min = border, max = 1-border)
	print("Logistic transformation del target")
else:
	Train_target_trans = Train_target
	Test_target_trans = Test_target

if regression:
	#########################REGRESSION
	#Optimizing hyperparameters
	grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=5, verbose=0, scoring="r2")
	#Fitting the model
	grid.fit(Train_data, list(Train_target_trans.values()))
	#Retreaving best estimator
	best = grid.best_estimator_
	#Making predicitons and evaluate performances
	Test_predictions = grid.predict(Test_data)
	Test_predictions = dict(zip(Test_data.index, Test_predictions))
	R2 = r2_score(list(Test_target_trans.values()), list(Test_predictions.values())) #il R2 ridotto non mi sembra adatto, non ho un numero preciso di parametri che fitto
	im.plot_results(list(Test_target_trans.values()), list(Test_predictions.values()), "M_learning_images/" + WHICH + ".png", R2 = R2)
	#Feature importance
	F_imp = im.feature_importance(Train_data, best, "M_learning_images/"+WHICH+"_feature_importance.png")
	##Ordered heatmap
	#im.heatmap_abundance(Test_data, Target = Test_target, F_imp = F_imp, Normalize = True, filename = "M_learning_images/"+WHICH+"_Test_heatmap.png")
	#im.heatmap_abundance(Train_data, Target = Train_target_trans, F_imp = F_imp, Normalize = True, filename = "M_learning_images/"+WHICH+"_Train_heatmap.png")

if Plotting_train_target:
	plt.clf()
	values = list(Train_target.values())
	plt.hist(values, bins='auto')  # Choose an appropriate number of bins or 'auto'
	plt.title('Histogram of Train_target values')
	plt.xlabel('Value')
	plt.ylabel('Frequency')
	plt.show()

if Zeros_analysis:
	#DENOISING USING PREDICTIONS OF ZERO RESISTANCE TRAIN SAMPLES
	if Transform_the_targets == "Log_shift":
		Z_train = {k : v for k , v in Train_target_trans.items() if v == np.log(shift)}
		Z_test = {k : v for k , v in Test_target_trans.items() if v == np.log(shift)}
	elif Transform_the_targets == "Logistic":
		Z_train = {k : v for k , v in Train_target_trans.items() if v == -np.log((1/border)-1)}
		Z_test = {k : v for k , v in Test_target_trans.items() if v == -np.log((1/border)-1)}
	else:
		Z_train = {k : v for k , v in Train_target_trans.items() if v == 0}
		Z_test = {k : v for k , v in Test_target_trans.items() if v == 0}
	Z_train_data = Train_data.loc[Z_train.keys()]
	Z_train_pred = grid.predict(Z_train_data)
	Z_train_pred = dict(zip(Z_train_data.index, Z_train_pred))
	
	#Noto che ci sono sue cluster di predizioni, separati dal valore -6
	#In quelli < -6 c'è sono Aureus_S, in quelli > -6 ci c'è Arcobacter e Arcobacter + Aureus_S
	Z_train_pred_low = {k:v for k,v in Z_train_pred.items() if v < -6}
	Z_train_data_low = im.read_metadata([k for k,v in Z_train_pred_low.items()])
	Batteri_in_low = {k: v["Species_Name"] for k,v in Z_train_data_low["Metadata"].items()}
	Z_train_pred_high = {k:v for k,v in Z_train_pred.items() if v >= -6}
	Z_train_data_high = im.read_metadata([k for k,v in Z_train_pred_high.items()])
	Batteri_in_high = {k: v["Species_Name"] for k,v in Z_train_data_high["Metadata"].items()}

	Z_test_data = Test_data.loc[Z_test.keys()]
	Z_test_pred = [v for k,v in Test_predictions.items() if k in Z_test.keys()]
	Z_test_pred = dict(zip(Z_test.keys(), Z_test_pred))

	plt.clf()
	plt.hist(Z_train_pred.values(), bins='auto', alpha=0.5, label='Train_set_predictions', color='blue')
	plt.hist(Z_test_pred.values(), bins='auto', alpha=0.5, label='Test_set_predictions', color='red')
	if Transform_the_targets == "Log_shift":
		plt.axvline(np.log(shift), color='blue', linestyle='dashed', linewidth=2)
	elif Transform_the_targets == "Logistic":
		plt.axvline(-np.log((1/border)-1), color='blue', linestyle='dashed', linewidth=2)
	plt.legend(loc='upper right')
	plt.show()

if False:

	A = im.heatmap_abundance(Train_data, Target=Train_target, F_imp = F_imp, Normalize = False, rempove_empthy_col=True, filename = None)
	B = im.heatmap_abundance(Test_data, Target=Test_target, F_imp = F_imp, Normalize = False, rempove_empthy_col=True, filename = None)
	A.max().max()
	B.max().max()
