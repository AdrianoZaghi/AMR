import json
import sys
import pandas as pd
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt

import M_learning_images as im

average_read_length = 270

def read_bwt_2_aro(popolazioni):
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
	return({"Pathognes":pd.DataFrame.from_dict(Path, orient="index").fillna(0),\
		"Antibiotics": pd.DataFrame.from_dict(Anti, orient="index").fillna(0),\
		"Drugs_class": pd.DataFrame.from_dict(Drug, orient="index").fillna(0)})

def read_metadata(popolazioni):
	metadata = {}
	for P in popolazioni:
		M_pop = pd.read_csv(P+"/input.tsv", sep="\t")
		metadata[P] = M_pop
	[metadata[m]["Phenotype"].replace(["susceptible", "resistant"], [0,1], inplace=True) for m in metadata]
	r_abundances = [float(metadata[pop_name][["Abundance", "Phenotype"]].prod(axis=1).sum(axis=0)) for pop_name in metadata]
	return {"Metadata":metadata, "Abundances":r_abundances}

def regression(Train, Test, Train_target, Test_target, classifier, parameter_space):
	missing_in_Train = [c for c in Test.columns if c not in Train.columns]
	missing_in_Test = [c for c in Train.columns if c not in Test.columns]
	Test[missing_in_Test] = 0
	Train[missing_in_Train] = 0
	Test=Test.reindex(sorted(Test.columns), axis=1)
	Train=Train.reindex(sorted(Train.columns), axis=1)
	grid = GridSearchCV(classifier, parameter_space, n_jobs=-1, refit=True).fit(Train, Train_target)
	best = grid.best_estimator_
	Test_predictions = best.predict(Test)
	print("Target \t Predictions")
	print(Test_target, Test_predictions)
	return({"Best estimator":best, "Predictions":Test_predictions})



Train_0 = ["Train_000", "Train_001", "Train_002", "Train_003", "Train_004", "Train_005", "Train_006", "Train_007", "Train_008", "Train_009"]
Test_0 = ["Test_000", "Test_001", "Test_002", "Test_003", "Test_004", "Test_005", "Test_006", "Test_007", "Test_008", "Test_009"]
Train_1 = ["Train_100", "Train_101", "Train_102", "Train_103", "Train_104", "Train_105", "Train_106", "Train_107", "Train_108", "Train_109"]
Test_1 = ["Test_100", "Test_101", "Test_102", "Test_103", "Test_104", "Test_105", "Test_106", "Test_107", "Test_108", "Test_109"]
Train_2 = ["Train_200", "Train_201", "Train_202", "Train_203", "Train_204", "Train_205", "Train_206", "Train_207", "Train_208", "Train_209"]
Test_2 = ["Test_200", "Test_201", "Test_202", "Test_203", "Test_204", "Test_205", "Test_206", "Test_207", "Test_208", "Test_209", "Test_210", "Test_211", "Test_212", "Test_213", "Test_214"]
Controlli = ["Controllo_Aureus_S" ,"Controllo_Aureus_R" ,"Controllo_Coli" ,"Controllo_Arcobacter" ,"Controllo_Pseudomonas"]

El = ElasticNet(max_iter=100000, tol=0.05)
params_El = {"alpha":[0.01]+list(range(1, 100, 10)), "l1_ratio":[0, .1, .5, .7, .9, .95, .99, 1]}

Forest = RandomForestRegressor(random_state = 0)
params_Forest = {"criterion":["squared_error", "absolute_error", "friedman_mse", "poisson"]}

Ada =  AdaBoostRegressor(n_estimators=100, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.5, 1, 1.5, 3, 5]}

Train = Train_0
Test = Test_0
Info_controlli = read_bwt_2_various(Controlli)
Train_data = read_bwt_2_aro(Train)
Test_data = read_bwt_2_aro(Test)
Train_target= read_metadata(Train)["Abundances"]
Test_target= read_metadata(Test)["Abundances"]
#[Info_controlli[c].to_csv("Controllo_"+c+".tsv", sep="\t") for c in Info_controlli]
El_results = regression(Train_data, Test_data, Train_target, Test_target, El, params_El)
im.plot_results(Test_target, El_results["Predictions"], "M_learning_images/prova")
Forest_results = regression(Train_data, Test_data, Train_target, Test_target, Forest, params_Forest)
Ada_results = regression(Train_data, Test_data, Train_target, Test_target, Ada, params_Ada)