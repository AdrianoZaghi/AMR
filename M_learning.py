import json
import sys
import pandas as pd
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.model_selection import GridSearchCV

average_read_length = 270

def read_bwt_2(popolazioni):
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
	#		print(P + " ha saltato " + str(len(dati_)-len(dati)) + " entrate su " + str(len(dati_)))
	Return = pd.DataFrame.from_dict(resistomes, orient="index").fillna(0)
	return Return

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
	print(Train)
	print(Test)
	grid = GridSearchCV(classifier, parameter_space, n_jobs=-1, refit=True).fit(Train, Test_target)
	best = grid.best_estimator_
	Test_predictions = best.predict(Test)
	print("Target \t Predictions")
	print(Test_target)
	print(Test_predictions)
	return(best)

Te=["Test_000", "Test_003", "Test_005", "Test_100","Train_103"]
Tr=["Train_000","Train_001","Train_005", "Train_102","Train_104"]
El = ElasticNet(max_iter=100000, tol=0.05)
params = {"alpha":[0.01]+list(range(1, 100, 10)), "l1_ratio":[0, .1, .5, .7, .9, .95, .99, 1]}

Train= read_bwt_2(Tr)
Test= read_bwt_2(Te)
Train_target= read_metadata(Tr)["Abundances"]
Test_target= read_metadata(Te)["Abundances"]
best = regression(Train, Test, Train_target, Test_target, El, params)