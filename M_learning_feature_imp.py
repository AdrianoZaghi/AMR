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

Base_estimator_for_adaboost = DecisionTreeRegressor(max_depth=10)
Ada =  AdaBoostRegressor(estimator=Base_estimator_for_adaboost,n_estimators=500, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.1, 0.5, 1, 2]}
#SELECTING THE MODEL TO USE

model = Ada
params = params_Ada

WHICH = "De_caesar_big"

t = campioni[campioni["Size"] == 0.1]
Train = t[t["Aureus_S_N_strains"] == 1]
Train = ["samples_2/" + x for x in Train[Train["Real/Simul"] == "De_Caesar_simul"].index]
#Train = Train + ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
Test = ["samples_2/" + x for x in t[t["Real/Simul"] == "Real"].index]


Train_target = im.read_metadata(Train)["Abundances_r"]
Test_target = im.read_metadata(Test)["Abundances_r"]

t_test = t[t["Real/Simul"] == "Real"]
t_test = t_test[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]

bias_De_caesar = {'Aureus_S': 4.96,
 'Aureus_R': 4.96,
 'Coli': 1,
 'Pseudomonas': 1,
 'Arcobacter': 1}

bias_simul = {'Aureus_S': 5.03,
'Aureus_R': 5.03,
'Coli': 1,
'Pseudomonas': 1,
'Arcobacter': 1}

t_test = t_test.mul(bias_simul, axis=1)
#Normalizzo le abbondanze aggiustate e sommo le abbondanze di Aureus_R, Coli e Pseudomonas per ottenere la resistenza
t_test = t_test.div(t_test.sum(axis=1), axis=0)
t_test = t_test[["Aureus_R", "Coli", "Pseudomonas"]].sum(axis=1)
Test_target_aggiustato = dict(zip(t_test.index, t_test.values))

Feature_importances_total = []
for i in range(5):
	model.random_state = i
	print(".", end = '')
	Test_data, Train_data = im.read_bwt_2_aro(Test,1), im.read_bwt_2_aro(Train,1)
	missing_in_Test, missing_in_Train = [c for c in Train_data.columns if c not in Test_data.columns], [c for c in Test_data.columns if c not in Train_data.columns]
	Test_data, Train_data = Test_data.drop(columns=missing_in_Train), Train_data.drop(columns=missing_in_Test)
	Test_data, Train_data = Test_data.sort_index(axis=1), Train_data.sort_index(axis=1)
	print("Dati letti.\nSto usando "+str(len(Train_data.columns))+" colonne nel fit")
	
	grid = GridSearchCV(model, params, n_jobs=-1, refit=True, cv=5, verbose=0, scoring="r2")
	#Fitting the model
	grid.fit(Train_data, list(Train_target.values()))
	#Retreaving best estimator
	#Making predicitons and evaluate performances

	F_imp = grid.best_estimator_.feature_importances_
	F = Train_data.columns
	F_imp = {F[i]:F_imp[i] for i in range(len(F))}
	Feature_importances_total.append(F_imp)

Features = pd.DataFrame(Feature_importances_total)
How_many_cols = 30
Sorted_col = list(Features.mean().sort_values(ascending=False).index)
Just_the_col_i_need = Sorted_col[:How_many_cols]
Just_the_col_i_need = Just_the_col_i_need[::-1]
sorted_df = Features[Just_the_col_i_need]

fig, ax = plt.subplots(figsize=(20, 12))
sorted_df.boxplot(ax=ax, vert=False)
ax.set_title('Feature importance assigned to CARD genes', fontsize=15)
ax.set_xlabel('Feature Importance', fontsize=18)
ax.set_ylabel('CARD gene identifier', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()

plt.savefig("M_learning_images/Feature_imp_"+WHICH+".png")