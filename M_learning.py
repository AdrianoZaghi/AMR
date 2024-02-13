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
campioni = pd.read_csv("Database_campioni.csv", index_col=0)

#Declaring the models and initializing their grid of hyperparameters
El = ElasticNet(max_iter=100000, tol=0.05)
params_El = {"alpha":[0, 0.001, 0.01]+list(range(1, 100, 10)), "l1_ratio":[0, .1, .5, .7, .9, .95, .99, 1]}

Forest = RandomForestRegressor(random_state = 0)
params_Forest = {"criterion":["squared_error", "absolute_error", "friedman_mse", "poisson"]}

Ada =  AdaBoostRegressor(n_estimators=1000, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.5, 1, 1.5, 2 ,3, 5,10]}

#%%
#SELECTING THE MODEL TO USE
model = Ada
params = params_Ada

#Deciding which samples to peack from the database and defining the train and test sets
WHICH = "AureusS_Arcobacter_AureusR_5_strain"
print(WHICH)
if(WHICH == "Aureus_1_strain"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] == 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus 1 strain train")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Aureus_5_strain"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [5,0] for x in t["N strains"]]]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] == 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus 5 strain train")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "AureusS_Arcobacter_AureusR_1_strain"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] != 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Aureus_R 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "AureusS_Arcobacter_Coli_1_strain"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] == 0]
	t = t[t["Coli"] != 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Coli 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "AureusS_Arcobacter_Pseudomonas_1_strain"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] == 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] != 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "AureusS_Arcobacter_AureusR_5_strain"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [5,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] != 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Aureus_R 5 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "AureusS_Arcobacter_AureusR_5_strain_0.1Gbp"):
	t = campioni[campioni["Size"] == 0.1]
	t = t[[x in [5,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] != 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Aureus_R 5 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "AureusS_Arcobacter_Coli_5_strain"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [5,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] == 0]
	t = t[t["Coli"] != 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Coli 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "AureusS_Arcobacter_e_Pseudomonas_5_strain"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [5,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] == 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] != 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Pseudomonas 5 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Tutti_tranne_AureusR"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] == 0]
	t = t[t["Coli"] != 0]
	t = t[t["Pseudomonas"] != 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Tutti_tranne_Coli"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] != 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] != 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Tutti_tranne_Pseudomonas"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] != 0]
	t = t[t["Coli"] != 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Tutti"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [5,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] != 0]
	t = t[t["Coli"] != 0]
	t = t[t["Pseudomonas"] != 0]
	t = t[t["Arcobacter"] != 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Synth_on_Synth_AureusS_Arcobacter_AureusR"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] != 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	Simul = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Synth_on_Synth_AureusS_Arcobacter_Coli"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] == 0]
	t = t[t["Coli"] != 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	Simul = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Synth_on_Synth_AureusS_Arcobacter_Pseudomonas"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] == 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] != 0]
	t = t[t["Arcobacter"] != 0]
	Simul = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Synth_on_Synth_No_Aureus_R"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] == 0]
	t = t[t["Coli"] != 0]
	t = t[t["Pseudomonas"] != 0]
	t = t[t["Arcobacter"] != 0]
	Simul = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Synth_on_Synth_No_Coli"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] != 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] != 0]
	t = t[t["Arcobacter"] != 0]
	Simul = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Synth_on_Synth_No_Pseudomonas"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] != 0]
	t = t[t["Coli"] != 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	Simul = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Sinth_on_synth_AureusS_Arcobacter_AureusR_5_strain_0.1Gbp"):
	t = campioni[campioni["Size"] == 0.1]
	t = t[[x in [5,0] for x in t["N strains"]]]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] != 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	Simul = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0) 
	print("Aureus_S e Arcobacter e Pseudomonas 1 strain")
	print("Quanti trains: " + str(len(Train)))
	print("Quanti tests: " + str(len(Test)))
elif(WHICH == "Confornto base"):
	Train = ["samples/Patric_" + x for x in ["Aureus_R", "Aureus_S", "Coli", "Pseudomonas", "Arcobacter"]]
	Test = ["samples/Decesare_" + x for x in ["Aureus_R", "Aureus_S", "Coli", "Pseudomonas", "Arcobacter"]]

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

#reading targets
Train_target= im.read_metadata(Train)["Abundances_r"]
Test_target= im.read_metadata(Test)["Abundances_r"]

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