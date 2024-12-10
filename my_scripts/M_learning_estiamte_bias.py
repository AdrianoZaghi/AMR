import json
import sys
import pandas as pd
from sklearn.linear_model import LinearRegression
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
from scipy.optimize import curve_fit

#Lunghezza media delle reads simulate e reali, comunque un fattore di normalizzazione che dovrebbe valere per entrambi
#potrebbe essere diverso tra simulate e reali
average_read_length = 270
images_dir = "M_learning_images/"
campioni = pd.read_csv("Database_campioni2.csv", index_col=0)

#%%
Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
Ada =  AdaBoostRegressor(base_estimator=Base_estimator_for_adaboost,n_estimators=100, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1, 0.5, 1, 2]}
#SELECTING THE MODEL TO USE
model = Ada
params = params_Ada
Transform_the_targets = None # "Log_shift" "Logistic" None
subsample_counts = 1
Eliminate_extreme_samples = False
Subsample_train_samples = 0
MANUAL_FEATURE_SELECTION = False
Plotting_train_target = False
Zeros_analysis = False
Add_second_degree_features = False
Compare_test_targets = False
regression = False
READ_TARGET_FROM =  "Database_campioni" # "Tax_Ettore" #


t = campioni[campioni["Size"] == 0.1]
Train = t[t["Aureus_S_N_strains"] == 1]
#Selecting train and test samples: data and target
#Train = ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
Train = ["samples_2/" + x for x in Train[Train["Real/Simul"] == "Simul"].index]
#Train = Train + ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
Test = ["samples_2/" + x for x in t[t["Real/Simul"] == "Real"].index]
Train_data = im.read_bwt_2_aro(Train, subsample_counts)
Train_target = im.read_metadata(Train)["Abundances_r"]
Test_data = im.read_bwt_2_aro(Test, subsample_counts)
Test_metadata = im.read_metadata(Test)
Test_target = Test_metadata["Abundances_r"]

#Dealing with columns absent in one of the two datasets 
missing_in_Train = [c for c in Test_data.columns if c not in Train_data.columns]
missing_in_Test = [c for c in Train_data.columns if c not in Test_data.columns]
Test_data = Test_data.drop(columns=missing_in_Train)
Train_data = Train_data.drop(columns=missing_in_Test)
Test_data = Test_data.sort_index(axis=1)
Train_data = Train_data.sort_index(axis=1)
print("Sto usando "+str(len(Train_data.columns))+" colonne nel fit")
grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=10, verbose=0, scoring="r2")
grid.fit(Train_data, list(Train_target.values()))
train_scores = {k:v for k,v in grid.cv_results_.items() if "split" in k}
best = grid.best_estimator_
Test_predictions = grid.predict(Test_data)
Test_predictions = dict(zip(Test_data.index, Test_predictions))
assert Test_predictions.keys() == Test_target.keys()

print("Retriving nominal abundances")
#Questi sono i valori delle predizioni dei miei modelli che, in questo caso, consideriamo come vere
nominal_abundances = [Test_metadata["Metadata"][k] for k in Test_target.keys()]
nominal = []
for x in nominal_abundances:
	x = x.set_index("Genome_ID")
	Ab = []
	Ab.append(x.loc["Aureus_S", "Abundance"]) if "Aureus_S" in x.index else Ab.append(0)
	Ab.append(x.loc["Aureus_R", "Abundance"]) if "Aureus_R" in x.index else Ab.append(0)
	Ab.append(x.loc["Coli", "Abundance"]) if "Coli" in x.index else Ab.append(0)
	Ab.append(x.loc["Pseudomonas", "Abundance"]) if "Pseudomonas" in x.index else Ab.append(0)
	Ab.append(x.loc["Arcobacter", "Abundance"]) if "Arcobacter" in x.index else Ab.append(0)
	nominal.append(Ab)
nominal = pd.DataFrame(nominal, columns=["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"])

print("Estimating the bias")
#FUNZIONE PER CORREGGERE IL BIAS AGGIUSTANDO PER UNA DIVERSA VITALITÀ DELLE COLTURE
#Su consiglio di Giampieri, fisso a 1 il parametro di Coli altrimenti il problema è sovradeterminato
#ho scelto Coli perchè sembra essere quello più sottostimato
def Vitality_adjustement(Nominal_abundances, x, y, w, z):
	Aureus_S, Aureus_R, Coli, Pseudomonas, Arcobacter = Nominal_abundances
	Real_total = Aureus_S*x + Aureus_R*y + Coli + Pseudomonas*w + Arcobacter*z
	Resistance = (Aureus_R*y + Coli + Pseudomonas*w)/Real_total
	return Resistance

#fUNZIONE PER CORREGGERE IL BIAS AGGIUSTANDO PER LA DIVERSA FACILITÀ DI ESTAZIONE TRA GRAM POSITIVI E GRAM NEGATIVI
def Gram_positiv_or_negativ_adjustement(Nominal_abundances, x):
	Aureus_S, Aureus_R, Coli, Pseudomonas, Arcobacter = Nominal_abundances
	Real_total = Aureus_S*x + Aureus_R*x + Coli + Pseudomonas + Arcobacter
	Resistance = (Aureus_R*x + Coli + Pseudomonas)/Real_total
	return Resistance

#Ottimizziamo sta funzione
#Valore iniziale dei parametri è 1,1,1,1
X =  [nominal["Aureus_S"], nominal[ "Aureus_R"], nominal["Coli"], nominal["Pseudomonas"], nominal["Arcobacter"]]
Resistance_ML = [Test_predictions[k] for k in Test_target.keys()]
fit = curve_fit(Gram_positiv_or_negativ_adjustement, X, Resistance_ML, p0=[1])
Adjustment = fit[0][0]
cov_matrix = fit[1]
covv = np.sqrt(cov_matrix.diagonal())
COR_matrix = cov_matrix/covv.reshape(-1,1) * covv.reshape(1,-1)

print("Adjustment factor: "+str(Adjustment))
print("Covariance matrix:")
print(COR_matrix)
print("Correlation matrix:")
print(COR_matrix/np.sqrt(COR_matrix.diagonal().reshape(-1,1) * COR_matrix.diagonal().reshape(1,-1)))