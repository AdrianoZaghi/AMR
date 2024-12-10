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

average_read_length = 270
images_dir = "M_learning_images/"

def read_bwt_2_aro(popolazioni):
	#Legge le abbondanze delle ARO entries rileate dall'alineamento di rgi bwt
	#Le aro entries sono identificativi dei resistomi presenti in rgi
	#Questa funzione restituisce un database con le abbondanze stimate di ogni reistoma per ogni campione
	resistomes = {}
	for P in popolazioni:
		with open(P+"/bwt_out/bwt.allele_mapping_data.json", mode = "r") as f:
			dati_ = json.load(f)
			dati = [di for di in dati_ if isinstance(di, dict)]
			R={}
			for aro in dati:
				if str(aro["aro_accession"]) in R.keys():
					R[str(aro["aro_accession"])] += int(aro["reads"]["all"]) * average_read_length / int(aro["reference"]["sequence_length"])
				else:
					R[str(aro["aro_accession"])] = int(aro["reads"]["all"]) * average_read_length / int(aro["reference"]["sequence_length"])
			resistomes[P] = R
	Return = pd.DataFrame.from_dict(resistomes, orient="index").fillna(0)
	return Return

def read_bwt_2_various(popolazioni):
	#Questa funzione restituisce dei dataframe che contengono features associate ai resistomi
	#Le "abbondanze" di qeuste feature sono la somma delle abondanze dei resistomi aventi quella feature nel campione
	#Le features sono:
	#	Patogeno a cui è associato il resistoma
	#	Antibiotico per cui viene sviluppata la resistenza
	#	Farmaco per cui viene sviluppata la resistenza
	Path = {}
	Anti = {}
	Drug = {}
	for P in popolazioni:
		with open(P+"/bwt_out/bwt.allele_mapping_data.json", mode = "r") as f:
			dati_ = json.load(f)
			dati = [di for di in dati_ if isinstance(di, dict)]
			Pathogens = {}
			Antibiotics = {}
			Drugs_class = {}
			for aro in dati:
				No_antibiotics_found_for_aros = 0
				No_drug_class_found_for_aros = 0
				if str(aro["observed_in_pathogens"]) in Pathogens.keys():
					Pathogens[str(aro["observed_in_pathogens"])] += int(aro["reads"]["all"]) * average_read_length / int(aro["reference"]["sequence_length"])
				else:
					Pathogens[str(aro["observed_in_pathogens"])] = int(aro["reads"]["all"]) * average_read_length / int(aro["reference"]["sequence_length"])
				if "Antibiotic" in aro["resistomes"].keys():
					for anti in aro["resistomes"]["Antibiotic"]:
						if str(anti) in Antibiotics.keys():
							Antibiotics[anti] += int(aro["reads"]["all"]) * average_read_length / int(aro["reference"]["sequence_length"])
						else:
							Antibiotics[anti] = int(aro["reads"]["all"]) * average_read_length / int(aro["reference"]["sequence_length"])
				else:
					No_antibiotics_found_for_aros += 1
				if "Drug Class" in aro["resistomes"].keys():
					for d_class in aro["resistomes"]["Drug Class"]:
						if str(d_class) in Drugs_class.keys():
							Drugs_class[d_class] += int(aro["reads"]["all"]) * average_read_length / int(aro["reference"]["sequence_length"])
						else:
							Drugs_class[d_class] = int(aro["reads"]["all"]) * average_read_length / int(aro["reference"]["sequence_length"])
				else:
					No_drug_class_found_for_aros += 1
		Path[P], Anti[P], Drug[P] = Pathogens, Antibiotics, Drugs_class
	print(No_drug_class_found_for_aros, "ARO non avevano drug class associate")
	print(No_antibiotics_found_for_aros, "ARO non avevano antibiotici associati")
	return({"Pathogens":pd.DataFrame.from_dict(Path, orient="index").fillna(0),\
		"Antibiotics": pd.DataFrame.from_dict(Anti, orient="index").fillna(0),\
		"Drugs_class": pd.DataFrame.from_dict(Drug, orient="index").fillna(0)})

def read_metadata(popolazioni):
	#Legge i metadati del campione, in particolare abbondanza di strains resistenti
	#DA MODIFICARE IN CASO
	metadata = {}
	for P in popolazioni:
		M_pop = pd.read_csv(P+"/input.tsv", sep="\t")
		metadata[P] = M_pop
	[metadata[m]["Phenotype"].replace(["susceptible", "resistant"], [0,1], inplace=True) for m in metadata]
	r_abundances = []
	Aureus_r_abundance = []
	for pop_name in metadata:
		M = metadata[pop_name]
		tot = sum(M["Abundance"])
		r_abundances.append(M[["Abundance", "Phenotype"]].prod(axis=1).sum(axis=0)/tot)
		if("Aureus_R" in M.values):
			Aureus_r_abundance.append(M[M["Genome_ID"]=="Aureus_R"].iloc[0]["Abundance"]/tot)
		else: Aureus_r_abundance.append(0)
	return {"Metadata":metadata, "Abundances_r":r_abundances, "Aureus_r": Aureus_r_abundance}

def regression(Train, Test, Train_target, classifier, parameter_space, Keep_columns_that_appear_just_in_one_set = False):
	#GESTISCE LA GRIDSEARCH PER ALGORITMI DI REGRESSIONE 
	missing_in_Train = [c for c in Test.columns if c not in Train.columns]
	missing_in_Test = [c for c in Train.columns if c not in Test.columns]
	if Keep_columns_that_appear_just_in_one_set:
		Test[missing_in_Test] = 0
		Train[missing_in_Train] = 0
	else:
		Test = Test.drop(columns=missing_in_Train)
		Train = Train.drop(columns=missing_in_Test)
	print("Sto usando le seguenti colonne nel fit")
	print(Train.columns)
	Test=Test.sort_index(axis=1)
	Train=Train.sort_index(axis=1)
	grid = GridSearchCV(classifier, parameter_space, n_jobs=-1, refit=True)
	grid.fit(Train, Train_target)
	best = grid.best_estimator_
	Test_predictions = best.predict(Test)
	#print("Target \t Predictions")
	#print(Test_target, Test_predictions)
	return({"Best estimator":best, "Predictions":Test_predictions})

campioni = pd.read_csv("Database_campioni.csv", index_col=0)

WHICH = "Synth on Synth Aureus_S, Arcobacter e Pseudomonas"

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
elif(WHICH == "Synth_on_Synth_from_real_genome_assembly_Aureus_S_Pseudomonas_Arcobacter"):
	t = campioni_2[campioni_2["Size"] == 0.01]
	t = t[t["Aureus_S"] != 0]
	t = t[t["Aureus_R"] == 0]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] != 0]
	t = t[t["Arcobacter"] != 0]
	Simul = ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
	Train, Test = train_test_split(Simul, test_size=0.1, random_state=0)
elif(WHICH == "Confornto base"):
	Train = ["samples/Patric_" + x for x in ["Aureus_R", "Aureus_S", "Coli", "Pseudomonas", "Arcobacter"]]
	Test = ["samples/Decesare_" + x for x in ["Aureus_R", "Aureus_S", "Coli", "Pseudomonas", "Arcobacter"]]

El = ElasticNet(max_iter=100000, tol=0.05)
params_El = {"alpha":[0, 0.001, 0.01]+list(range(1, 100, 10)), "l1_ratio":[0, .1, .5, .7, .9, .95, .99, 1]}

Forest = RandomForestRegressor(random_state = 0)
params_Forest = {"criterion":["squared_error", "absolute_error", "friedman_mse", "poisson"]}

Ada =  AdaBoostRegressor(n_estimators=1000, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.5, 1, 1.5, 2 ,3, 5,10]}

def log_ods(vector):
	#	Ci metto un -50 per evitare che il log di 0 dia errore
	ods = [vv/(1 - vv) for vv in vector]
	log_ods = [np.log(vv) if vv < 20 else np.log(20,) for vv in ods]
	return([vv if vv > -4 else -4 for vv in log_ods])

def log_shift(vector):
	return(np.log([v+1 for v in vector]))

def normalize_counts(data):
	r = data.apply(lambda x: x / x.max())

Train_data = read_bwt_2_aro(Train)
Test_data = read_bwt_2_aro(Test)
MANUAL_FEATURE_SELECTION = True
if MANUAL_FEATURE_SELECTION:
	selected = ["3000216"]#, "3004517"]
	Train_data = Train_data[selected]
	Test_data = Test_data[selected]
Train_target= dict(zip(Train, read_metadata(Train)["Abundances_r"]))
Test_target= dict(zip(Test, read_metadata(Test)["Abundances_r"]))

#Results = regression(Train_data, Test_data, list(Train_target.values()), El, params_El)
#R2 = r2_score(list(Test_target.values()), Results["Predictions"]) #il R2 ridotto non mi sembra adatto, non ho un numero preciso di parametri che fitto
#im.plot_results(list(Test_target.values()), Results["Predictions"], "M_learning_images/prova", R2 = R2)

missing_in_Train = [c for c in Test_data.columns if c not in Train_data.columns]
missing_in_Test = [c for c in Train_data.columns if c not in Test_data.columns]
if False:
	Test_data[missing_in_Test] = 0
	Train_data[missing_in_Train] = 0
else:
	Test_data = Test_data.drop(columns=missing_in_Train)
	Train_data = Train_data.drop(columns=missing_in_Test)
print("Sto usando le seguenti colonne nel fit")
print(Train_data.columns)
grid = GridSearchCV(El, params_El, n_jobs=-1, refit=True)
grid.fit(Train_data, list(Train_target.values()))
best = grid.best_estimator_
best.fit(Train_data, list(Train_target.values()))
Test_predictions = best.predict(Test_data)

R2 = r2_score(list(Test_target.values()), Test_predictions) #il R2 ridotto non mi sembra adatto, non ho un numero preciso di parametri che fitto
im.plot_results(list(Test_target.values()), Test_predictions, "M_learning_images/prova", R2 = R2)


#PLOT FEATURE IMPORTANCE
Feature_im = False
if(Feature_im):
	if(isinstance(Results["Best estimator"],ElasticNet)):
		F_imp = Results["Best estimator"].coef_
	else:
		F_imp = Results["Best estimator"].feature_importances_
	F = Train_data.columns
	F_imp = {F[i]:F_imp[i] for i in range(len(F))}
	F_imp = dict(sorted(F_imp.items(), key=lambda item: item[1]))
	F_imp_head = {i : F_imp[i] for i in list(F_imp.keys())[-15:]}
	plt.figure(figsize=(15, 6))
	plt.bar(F_imp_head.keys(), F_imp_head.values())

#HEATMAP
if (True):
	if False:
		tt = {"data": Test_data, "target": Test_target}
	else:
		tt = {"data": Train_data, "target": Train_target}
	print(tt)
	sorted_keys = sorted(tt["target"], key=tt["target"].get, reverse=True)
	print({x : tt["target"][x] for x in sorted_keys})
	d = tt["data"].reindex(sorted_keys)
	if False:
		d = d[F_imp.keys()]
	d_norm = d.apply(lambda x: x / x.max())
	plt.figure(figsize=(15, 6))
	sns.heatmap(d, cmap='viridis')
	plt.figure(figsize=(15, 6))
	sns.heatmap(d_norm, cmap='viridis')

if False:
	Info_controlli = read_bwt_2_various(Train)
	im.plot_controls(Info_controlli["Antibiotics"], images_dir+"Control_antibiotics.pdf")
	im.plot_controls(Info_controlli["Drugs_class"], images_dir+"Control_drugs.pdf")
	im.plot_controls(Info_controlli["Pathogens"], images_dir+"Control_pathogens.pdf")
	[Info_controlli[c].to_csv("Controllo_"+c+".tsv", sep="\t") for c in Info_controlli]


#ci sono dei geni che non compaiono mai in uno dei due stes.
#prova a normalizzare i conteggi, però probabilmente andrebbe fatto considerando la massima abbondanza di train e test
#poi fai selezione manuale delle features