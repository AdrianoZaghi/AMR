import json
import sys
import pandas as pd
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
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

#Declaring the models and initializing their grid of hyperparameters
El = ElasticNet(max_iter=100000, tol=0.05)
params_El = {"alpha":[0, 0.001, 0.01]+list(range(1, 100, 10)), "l1_ratio":[0, .1, .5, .7, .9, .95, .99, 1]}

Forest = RandomForestRegressor(random_state = 0)
params_Forest = {"criterion":["squared_error", "absolute_error", "friedman_mse", "poisson"]}

Ada =  AdaBoostRegressor(n_estimators=1000, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.5, 1, 1.5, 2 ,3, 5,10]}


#SELECTING THE MODEL TO USE
model = Ada
params = params_Ada

#Deciding which samples to peack from the database and defining the train and test sets
WHICH = "Synth_on_Real_patric_genomes"
print(WHICH)
if(WHICH == "Synth_on_Synth_assembled_genomes"):
	t = campioni[campioni["Size"] == 0.01]
	Simul = ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Synth_on_Synth_assembled_genomes")
elif(WHICH == "Synth_on_Real_assembled_genomes"):
	t = campioni[campioni["Size"] == 0.01]
	Train = ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
	Test = ["samples_2/" + x for x in t[t["Real/Simul"] == "Real"].index]
	print("Synth_on_Real_assembled_genomes")
elif(WHICH == "Synth_on_Synth_patric_genomes"):
	t = campioni[campioni["Size"] == 0.01]
	Simul = ["samples_2/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Synth_on_Synth_patric_genomes")
elif(WHICH == "Synth_on_Real_patric_genomes"):
	t = campioni[campioni["Size"] == 0.01]
	Train = ["samples_2/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Test = ["samples_2/" + x for x in t[t["Real/Simul"] == "Real"].index]
	print("Synth_on_Real_patric_genomes")
print("Quanti trains: " + str(len(Train)))
print("Quanti tests: " + str(len(Test)))

#################READING DATA
#Features (this funciton is also scaling the features to make them relatable between samples)
Train_data = im.read_bwt_2_aro(Train)
Test_data = im.read_bwt_2_aro(Test)
#Dealing with columns absent in one of the two datasets 
missing_in_Train = [c for c in Test_data.columns if c not in Train_data.columns]
missing_in_Test = [c for c in Train_data.columns if c not in Test_data.columns]
if False:
	Test_data[missing_in_Test] = 0
	Train_data[missing_in_Train] = 0
else:
	Test_data = Test_data.drop(columns=missing_in_Train)
	Train_data = Train_data.drop(columns=missing_in_Test)
Test_data = Test_data.sort_index(axis=1)
Train_data = Train_data.sort_index(axis=1)
print("Sto usando le seguenti colonne nel fit")
print(Train_data.columns)
#Manual feature selection
MANUAL_FEATURE_SELECTION = False
if MANUAL_FEATURE_SELECTION:
	selected = ["3000216"]#, "3004517"]
	Train_data = Train_data[selected]
	Test_data = Test_data[selected]

#%%
#reading targets
Train_target= im.read_metadata(Train)["Abundances_r"]
Test_target= im.read_metadata(Test)["Abundances_r"]

Target_transformation = "Niente"
if Target_transformation == "Log_shift":
	Train_target = im.log_shift(Train_target, shift = 1)
	Test_target = im.log_shift(Test_target, shift = 1)
	print("Log_shift del target")
elif Target_transformation == "Log_odds":
	Train_target = im.log_odds(Train_target, upper_saturation=0.9999999999, lower_saturation=0.0000000001)
	Test_target = im.log_odds(Test_target, upper_saturation=0.9999999999, lower_saturation=0.0000000001)
	print("Log_odds del target")

#########################REGRESSION
#Optimizing hyperparameters
grid = GridSearchCV(model, params, n_jobs=-1, refit=True)
#Fitting the model
grid.fit(Train_data, list(Train_target.values()))
#Retreaving best estimator
best = grid.best_estimator_
#Making predicitons and evaluate performances
Test_predictions = grid.predict(Test_data)
R2 = r2_score(list(Test_target.values()), Test_predictions) #il R2 ridotto non mi sembra adatto, non ho un numero preciso di parametri che fitto
im.plot_results(list(Test_target.values()), Test_predictions, "M_learning_images/" + WHICH + ".png", R2 = R2)
#Feature importance
F_imp = im.feature_importance(Train_data, best, "M_learning_images/"+WHICH+"_feature_importance.png")
#Ordered heatmap
im.heatmap_abundance(Test_data, Target = Test_target, F_imp = F_imp, Normalize = True, filename = "M_learning_images/"+WHICH+"_Test_heatmap.png")
im.heatmap_abundance(Train_data, Target = Train_target, F_imp = F_imp, Normalize = True, filename = "M_learning_images/"+WHICH+"_Train_heatmap.png")

#VALUTARE LA BONTÀ DEI NOSTRI RISULTATI:
# - R2 PROBABILMENTE NON È IL MIGLIORE PERCHÈ NON DA TROPPA INMPORTANZA A VARIAZIONI PICCOLE SU VALORI PICCOLI.
# - PROVATE A FARE I LOGARITMI DEI VALORI DI ABUNDANCE E VEDERE SE IL R2 MIGLORA