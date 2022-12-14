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
from sklearn.linear_model import ElasticNetCV
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

from sklearn.model_selection import train_test_split

average_read_length = 170

def read_kmer(popolazioni):
	"""This function reads the output of kmc_count
	The rest of the analysis is not already supported for this kind of data
	"""
	kmers = {}
	for P in popolazioni:
		if os.path.isfile(P+"/kmc_out/kmer_count.txt"):
			with open(P+"/kmc_out/kmer_count.txt") as d:
				D = pd.read_csv(d)
			kmers[P] = D.to_dict()
		elif os.path.isfile(P+"/kmc_out/kmer_count_1.txt") and os.path.isfile(P+"/kmc_out/kmer_count_2.txt"):
			with open(P+"/kmc_out/kmer_count_1.txt", mode = "r") as uno, open(P+"/kmc_out/kmer_count_2.txt", mode = "r") as due:
				Uno = pd.read_csv(uno, sep='\t', index_col=0)
				Due = pd.read_csv(due, sep='\t', index_col=0)
				Tot = pd.concat([Uno, Due], axis = 1).fillna(0)
				Sum = Tot.sum(axis=1)
				Tot = pd.concat([Tot, Sum], axis=1)
				Tot.columns = ["Forward", "Backward", "Sum"]
			kmers[P] = Tot.to_dict()["Sum"]
		else:
			print("Non ho trovato kmers file")
			return 0
	Return = pd.DataFrame.from_dict(kmers, orient="index").fillna(0)
	return Return
	
def read_bwt_2(popolazioni):
	"""Read the output of rgi
	Paramethers:
		bwt_dir: the path to the directory where the files to be read are stored
	Returns:
		a two dimensional dictionary: at first level we have a sequence of sample_name:resistome
		each resistome have the structure ARO_identity:FPKB_estimated_abundance
	"""
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
	"""Read the metadata file
	Paramethers:
		meta_dir:
			The path to the directory where the metadata are stored.
			Notice that each metadata file have to be named exaccly equal to the data file (from rgi_bwt or kmc_count) 
		key:
			This key is refered to the structure of the json file used as metadata.
			Specify which value to consider as a label in the machine learning procedure that will follow, "phenotype" in our case
	Returns:
		A dictionary with the structure ARO_identity: phenotype
	"""
	metadata = {}
	for P in popolazioni:
		M_pop = pd.read_csv(P+"/input.tsv", sep="\t")
		metadata[P] = M_pop
	return metadata


def view_database(dataset, filename, view_a_spectrum = None):
	"""Produces various plots of the dataset provided as argument and also of a specific sample
	if "view_a_spectrum" have a reasonable integer value
	"""
	#PLOT CON LA SOMMA DEI VALORI NELLE COLONNE
	V = dataset.sum(axis = 1)
	dic = {a : V[a] for a in V.index.tolist()}
	plt.figure(0)
	plt.bar(dic.keys(),dic.values())
	plt.xlabel("Samples")
	plt.ylabel("Sample integral")
	plt.xticks([])
	plt.savefig("M_learning_images/"+filename+"_Integral.png")

	#SCATTERPLOT SPARSITà E MEDIANA
	plot_zoom=1
	y = dataset.median(axis=0)
	y = [i+0.1 if i < 0.1 else i for i in list(y)]
	x = [(c_data == 0).astype(int).sum()*100/(dataset.shape[0]-1) for (c_name, c_data) in dataset.iteritems()]
	colors = np.random.rand(3)
	somma = dataset.sum(axis = 0)/50
	somma = [l*plot_zoom for l in list(somma)]
	plt.figure(1)
	#plt.scatter(x, y, s=somma, c=colors, alpha=0.5)
	plt.scatter(x, y, c=colors, alpha=0.5)
	plt.yscale("log")
	plt.ylabel("Median")
	plt.xlabel("Percentage of sparsity")
	plt.title("Filtering results")
	#	plt.plot([0,100],[median_filter,median_filter], color = "red") #median
	#	plt.plot([sparsity_filter,sparsity_filter],[0,1500], color = "red") #sparsity
	plt.savefig("M_learning_images/"+filename+"_Scatter.png")
	
	#BARPLOT DI UN SINGOLO GRAFICO
	if isinstance(view_a_spectrum, int):
		V = dataset.drop("phenotype", axis = 1)
		dic_2 = V.iloc[view_a_spectrum,:].to_dict()
		dic_2 = dict(sorted(dic_2.items(), key=lambda item: item[1]))
		plt.bar(dic_2.keys(),dic_2.values())
		plt.xlabel("ARO indentities")
		plt.ylabel("Abundance (FPKB)")
		plt.xticks([])
		plt.show()
		
def CLR_normalization(dataset):
	"""Performs the CLR normalization on the dataset in argument
	The original dataset is not modified	
	"""
	temporary = dataset.replace(to_replace = 0, value = 0.5)
	temporary = temporary.drop("phenotype", axis = 1)
	ph = dataset["phenotype"].tolist()
	G = pow(temporary.product(axis = 1, skipna = True), 1/dataset.shape[1])
	temporary = np.log(temporary.div(G, axis = 0))
	temporary["phenotype"] = ph
	return temporary

def PCA_maker(dataset, numero_componenti = None, plots = False):
	"""Performs the PCA on the dataset
	It can be set a maximal number of components to be kept, otherwise we will end up with a number of them equal to the number of samples
	If plot is set thrue the function produces also some usefull plots
	The original dataset is not modified
	"""
	pca = PCA()	
	pca.fit(dataset.drop("phenotype", axis = 1).T)
	D = pd.DataFrame(pca.components_, columns=dataset.index.values.tolist()).T
	D.columns = ["component "+str(i) for i in range(D.shape[1])]
	if plots:
		V = pd.DataFrame(pca.explained_variance_ratio_)
		a = [str(v) for v in list(V.index)]
		b = list(V[0])
		plt.bar(a,b)
		plt.xlabel("PCA components")
		plt.ylabel("Explained variance rateo")
		plt.show()	
		plt.scatter(D.iloc[:,0], D.iloc[:,1], color = ["red" if pheno == "resistant" else "blue" for pheno in list(dataset["phenotype"])])
		plt.xlabel("First component")
		plt.ylabel("Second component")
		plt.show()
	D["phenotype"]=dataset["phenotype"].tolist()
	return D
	
def std_scale(dataset):
	"""Performs the standard scaling of the valueds in the dataframe
	The original dataset is not modified
	"""
	scale = StandardScaler()
	nuovo = scale.fit_transform(dataset.drop("phenotype", axis = 1))
	nuovo = pd.DataFrame(nuovo, index = list(dataset.index), columns = list(dataset.drop("phenotype", axis = 1).columns))
	nuovo.insert(nuovo.shape[1], "phenotype", dataset.loc[:,"phenotype"])
	return nuovo

def loo(classifier, resistome, phenotype, plot = False, linear = False):
	predictions = list(cross_val_predict(classifier, resistome, phenotype, cv=LeaveOneOut(), n_jobs=-1))
	CROSS = cross_validate(classifier, resistome, phenotype, cv=LeaveOneOut(), n_jobs=-1, return_estimator=True)
	return {"Scores":CROSS["test_score"], "Predictions": predictions, "Models": CROSS["estimator"], "Feature importance" : "DA IMPLEMENTARE"}

def cv5(classifier, resistome, phenotype, plot = False, linear = False):
	predictions = list(cross_val_predict(classifier, resistome, phenotype, cv=5, n_jobs=-1))
	CROSS = cross_validate(classifier, resistome, phenotype, cv=5, n_jobs=-1, return_estimator=True)
	return {"Scores":CROSS["test_score"], "Predictions and phenotype": [[predictions[i], phenotype[i]] for i in range(len(predictions))], "Models": CROSS["estimator"], "Feature importance" : "DA IMPLEMENTARE"}

pop_1_specie_random_resistance = ["population_0","population_1","population_2","population_3","population_4","population_5","population_6","population_7","population_8","population_9","population_10","population_11","population_12","population_13","population_14","population_15","population_16","population_17","population_18","population_19"]
pop_4_specie_random_resistance = ["population_20","population_21","population_22","population_23","population_24","population_25","population_26","population_27","population_28","population_29","population_30","population_31","population_32","population_33","population_34","population_35","population_36","population_37","population_38","population_39"]
pop_1_specie_si_no_resistance = ["population_40","population_41","population_42","population_43","population_44","population_45","population_46","population_47","population_48","population_49","population_50","population_51","population_52","population_53","population_54","population_55","population_56","population_57","population_58","population_59"]
pop_4_specie_si_no_resistance = ["population_60","population_61","population_62","population_63","population_64","population_65","population_66","population_67","population_68","population_69","population_70","population_71","population_72","population_73","population_74","population_75","population_76","population_77","population_78","population_79"]

Data_bwt = read_bwt_2(pop_1_specie_random_resistance)
print("Data_bwt created")
Data_kmer = read_kmer(pop_1_specie_random_resistance)
print("Data_kmer created")
Metadata = read_metadata(pop_1_specie_random_resistance)

view_database(Data_bwt, "FPKB_sum_per_sample")
view_database(Data_kmer, "Kmer_abundance_per_sample")

for m in Metadata:
	Metadata[m]["Phenotype"].replace(["susceptible", "resistant"], [0,1], inplace=True)
R_abundances = pd.DataFrame.from_dict({pop_name : float(Metadata[pop_name][["Abundance", "Phenotype"]].prod(axis=1).sum(axis=0)) for pop_name in Metadata}, orient="index")
r_abundances = [float(Metadata[pop_name][["Abundance", "Phenotype"]].prod(axis=1).sum(axis=0)) for pop_name in Metadata]


print("\nElastic net logistic regression and classification, results in 'El_cv5_results _bwt _kmer' varible")
El_1 = LogisticRegression(penalty='elastic', l1_ratio= [.1, .5, .7, .9, .95, .99, 1], eps=0.001, n_alphas=100, n_jobs=-1, max_iter=100000, tol=0.005)
El_regressor_results_bwt= cv5(El_1, Data_bwt, r_abundances)
El_2 = LogisticRegression(l1_ratio= [.1, .5, .7, .9, .95, .99, 1], eps=0.001, n_alphas=100, n_jobs=-1, max_iter=100000, tol=0.005)
El_regressor_results_kmer= cv5(El_2, Data_kmer, r_abundances)

print("\nRandom forest classification, results in 'Forest_classifier_results _bwt _kmer' variable")
Forest_1 = RandomForestClassifier(random_state = 0)
Forest_classifier_results_bwt = cv5(Forest_1, Data_bwt, r_abundances)
Forest_2 = RandomForestClassifier(random_state = 0)
Forest_classifier_results_kmer= cv5(Forest_2, Data_kmer, r_abundances)

print("\nAda boost classification, results in 'Ada_classifier_results _bwt _kmer' varible")
Ada_1 = AdaBoostClassifier(n_estimators=100, random_state=0)
Ada_regressor_results_bwt = cv5(Ada_1, Data_bwt, r_abundances)
Ada_2 = AdaBoostClassifier(n_estimators=100, random_state=0)
Ada_regressor_results_kmer= cv5(Ada_2, Data_kmer, r_abundances)

print("\nElastic net linear regression, results in 'El_regressor_result _bwt _kmer' varible")
El_1 = ElasticNetCV(l1_ratio= [.1, .5, .7, .9, .95, .99, 1], eps=0.001, n_alphas=100, n_jobs=-1, max_iter=100000, tol=0.005)
El_regressor_results_bwt= cv5(El_1, Data_bwt, r_abundances)
El_2 = ElasticNetCV(l1_ratio= [.1, .5, .7, .9, .95, .99, 1], eps=0.001, n_alphas=100, n_jobs=-1, max_iter=100000, tol=0.005)
El_regressor_results_kmer= cv5(El_2, Data_kmer, r_abundances)

print("\nRandom forest regression, results in 'Forest_regressor_results _bwt _kmer' variable")
Forest_1 = RandomForestRegressor(random_state = 0)
Forest_regressor_results_bwt = cv5(Forest_1, Data_bwt, r_abundances)
Forest_2 = RandomForestRegressor(random_state = 0)
Forest_regressor_results_kmer= cv5(Forest_2, Data_kmer, r_abundances)

print("\nAda boost regression, results in 'Ada_regressor_results _bwt _kmer' varible")
Ada_1 = AdaBoostRegressor(n_estimators=100, random_state=0)
Ada_regressor_results_bwt = cv5(Ada_1, Data_bwt, r_abundances)
Ada_2 = AdaBoostRegressor(n_estimators=100, random_state=0)
Ada_regressor_results_kmer= cv5(Ada_2, Data_kmer, r_abundances)