#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 14:57:20 2022

@author: utente
"""

import os.path
import json
import pandas as pd
import modin.pandas as PD
import numpy as np
import matplotlib.pyplot as plt
import statistics
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import AdaBoostRegressor
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_validate
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import GridSearchCV
from multiprocessing import Pool

from sklearn.model_selection import train_test_split

average_read_length = 170

def filter_kmer(dataframe, limit):
	R = dataframe[(dataframe.abundance > limit)]
	#print("ho droppato "+str(dataframe.shape[0]-R.shape[0])+" kmers delle "+str(dataframe.shape[0])+" che erano; limite = "+str(limit))
	return R

def read_single_kmer_file(P, limit = 10):
	if os.path.isfile(P+"/kmc_out/kmer_count.txt"):
		with open(P+"/kmc_out/kmer_count.txt") as d:
			D = pd.read_csv(d)
			D.columns = ["abundance"]
			D = filter_kmer(D, limit)
		#print("ho finito di leggere il kmer file di "+P)
		return D.to_dict()["abundance"]
	elif os.path.isfile(P+"/kmc_out/kmer_count_1.txt") and os.path.isfile(P+"/kmc_out/kmer_count_2.txt"):
		with open(P+"/kmc_out/kmer_count_1.txt", mode = "r") as uno, open(P+"/kmc_out/kmer_count_2.txt", mode = "r") as due:
			Uno = pd.read_csv(uno, sep='\t', index_col=0)
			Due = pd.read_csv(due, sep='\t', index_col=0)
			Tot = pd.concat([Uno, Due], axis = 1).fillna(0)
			Sum = Tot.sum(axis=1).to_frame()
			Sum.columns = ["abundance"]
			Sum = filter_kmer(Sum, limit)
		#print("ho finito di leggere il kmer file di "+P)
		return Sum.to_dict()["abundance"]
	else:
		print("Non ho trovato kmers file")
		return {"ho scazzato":0,"tutto":1}
	
def read_kmer_parallel(popolazioni):
	with Pool(10) as piscina:
		list_of_dictionary_counts = piscina.map(read_single_kmer_file, popolazioni)
	r = pd.DataFrame(list_of_dictionary_counts).fillna(0)
	r["index column"] = popolazioni
	return r.set_index("index column")

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
	if set(r_abundances) == {0,0.25}:
		l=[0 if x == 0 else 1 for x in r_abundances]
		r_abundances = l
	return {"Metadata":metadata, "Abundances":r_abundances}

def cv5(classifier, resistome, phenotype):
	predictions = list(cross_val_predict(classifier, resistome, phenotype, cv=5, n_jobs=-1))
	CROSS = cross_validate(classifier, resistome, phenotype, cv=5, n_jobs=-1, return_estimator=True)
	#return {"Scores":CROSS["test_score"], "Predictions and phenotype": [[predictions[i], phenotype[i]] for i in range(len(predictions))], "Models": CROSS["estimator"], "Feature importance" : "DA IMPLEMENTARE"}
	return {"Scores":CROSS["test_score"], "Predictions and phenotype": [[predictions[i], phenotype[i]] for i in range(len(predictions))], "Models":CROSS["estimator"]}
	
def plot_parameter_space(Data, abundances, classifier, ratio_param, intensity_param , param_space,filename):
	grid = GridSearchCV(classifier, param_space, n_jobs=-1)
	grid.fit(Data, abundances)
	scores = grid.cv_results_["mean_test_score"].reshape(len(param_space[intensity_param]), len(param_space[ratio_param]))
	print(scores)
	z_min, z_max = scores.min()-0.01, scores.max()	
	plt.pcolormesh(param_space[ratio_param],param_space[intensity_param],scores,vmin=z_min, vmax=z_max)
	plt.colorbar()
	plt.xlabel("L1 ratio")
	plt.ylabel("Penalization intensity")
	plt.savefig("M_learning_images/"+filename+".png")

def classification(Data, phenotype, filename=""):
	print("Starting logistic classification")
	El = LogisticRegression(penalty='elasticnet', solver="saga", n_jobs=-1, max_iter=100000, tol=0.005)
	print("Scanning paramether space")
	params = {"C":[0.01]+list(range(1, 100, 10)), "l1_ratio":[0, .1, .5, .7, .9, .95, .99, 1]}
	plot_parameter_space(Data, phenotype, El, "l1_ratio", "C", params, "Logistic_penalty_params")
	grid = GridSearchCV(El, params, n_jobs=-1)
	El_results= cv5(grid, Data, phenotype)
	print("Random forest classification")
	Forest = RandomForestClassifier(random_state = 0)
	Forest_results = cv5(Forest, Data, phenotype)
	print("Ada boost classification")
	Ada = AdaBoostClassifier(n_estimators=100, random_state=0)
	Ada_results = cv5(Ada, Data, phenotype)
	return {"EL_classification_results":El_results, "Forest_classification_results":Forest_results, "Ada_classification_results":Ada_results}

def regression(Data, abundances, filename=""):
	print("Starting elastic net linear regression")
	El = ElasticNet(max_iter=100000, tol=0.005)
	print("Scanning paramether space")
	params = {"alpha":[0.01]+list(range(1, 100, 10)), "l1_ratio":[0, .1, .5, .7, .9, .95, .99, 1]}
	plot_parameter_space(Data, abundances, El, "l1_ratio", "alpha", params, "Linear_penalty_params")
	grid = GridSearchCV(El, params, n_jobs=-1)
	El_results= cv5(grid, Data, abundances)
	print("Random forest regression")
	Forest = RandomForestRegressor(random_state = 0)
	Forest_results = cv5(Forest, Data, abundances)
	print("Ada boost regression")
	Ada = AdaBoostRegressor(n_estimators=100, random_state=0)
	Ada_results = cv5(Ada, Data, abundances)
	return {"EL_regression_results":El_results, "Forest_regression_results":Forest_results, "Ada_regression_results":Ada_results}

une_random = ["population_0","population_1","population_2","population_3","population_4","population_5","population_6","population_7","population_8","population_9","population_10","population_11","population_12","population_13","population_14","population_15","population_16","population_17","population_18","population_19"]
four_random = ["population_20","population_21","population_22","population_23","population_24","population_25","population_26","population_27","population_28","population_29","population_30","population_31","population_32","population_33","population_34","population_35","population_36","population_37","population_38","population_39"]
one_si_no = ["population_40","population_41","population_42","population_43","population_44","population_45","population_46","population_47","population_48","population_49","population_50","population_51","population_52","population_53","population_54","population_55","population_56","population_57","population_58","population_59"]
four_si_no = ["population_60","population_61","population_62","population_63","population_64","population_65","population_66","population_67","population_68","population_69","population_70","population_71","population_72","population_73","population_74","population_75","population_76","population_77","population_78","population_79"]

regression(read_bwt_2(one_random), read_metadata(one_random)["Abundances"], "1_random_bwt")
regression(read_bwt_2(four_random), read_metadata(four_random)["Abundances"], "4_random_bwt")
regression(read_kmer_2(one_random), read_metadata(one_random)["Abundances"], "1_random_kmer")
regression(read_kmer_2(four_random), read_metadata(four_random)["Abundances"], "4_random_kmer")
classification(read_bwt_2(one_si_no), read_metadata(one_si_no)["Abundances"], "1_si_no_bwt")
classification(read_bwt_2(four_si_no), read_metadata(four_si_no)["Abundances"], "4_si_no_bwt")
classification(read_kmer_2(one_si_no), read_metadata(one_si_no)["Abundances"], "1_si_no_kmer")
classification(read_kmer_2(four_si_no), read_metadata(four_si_no)["Abundances"], "4_si_no_kmer")