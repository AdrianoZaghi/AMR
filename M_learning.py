#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 14:57:20 2022

@author: utente
"""

import os
import json
import pandas as pd
import modin.pandas as PD
import numpy as np
import matplotlib.pyplot as plt
import statistics
from os.path import isfile, join
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

from sklearn.model_selection import train_test_split

average_read_length = 170

def read_kmer(popolazioni):
	"""This function reads the output of kmc_count
	The rest of the analysis is not already supported for this kind of data
	"""
	kmers = {}
	for P in popolazioni:
		with open(P+"/kmc_out/kmer_count_1.txt", mode = "r") as uno, open(P+"/kmc_out/kmer_count_2.txt", mode = "r") as due:
			Uno = pd.read_csv(uno, sep='\t', index_col=0)
			Due = pd.read_csv(due, sep='\t', index_col=0)
			Tot = pd.concat([Uno, Due], axis = 1).fillna(0)
			Sum = Tot.sum(axis=1)
			Tot = pd.concat([Tot, Sum], axis=1)
			Tot.columns = ["Forward", "Backward", "Sum"]
		kmers[P] = Tot.to_dict()["Sum"]
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
#			print(P + " ha saltato " + str(len(dati_)-len(dati)) + " entrate su " + str(len(dati_)))
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

def filter_sparcity_median(dataset, sparsity_filter, median_filter, plot = False, plot_zoom = 1):
	"""Perform the filter operations on the database
	The values for the median and the sparsity filters are the two positional arguments, while the standard arguments regulates the plots of the function
	See my thesis for a bether explenation of the effects of those filters
	The original dataset is not modified.
	"""
	new = dataset.copy()
	if plot:
		y = dataset.median(axis=0)
		y = [i+0.1 if i < 0.1 else i for i in list(y)]
		x = [(c_data == 0).astype(int).sum()*100/(new.shape[0]-1) for (c_name, c_data) in new.iteritems()]
		colors = np.random.rand(3)
		somma = dataset.drop("phenotype", axis=1).sum(axis = 0)/50
		somma = [l*plot_zoom for l in list(somma)]
		plt.scatter(x, y, s=somma, c=colors, alpha=0.5)
		plt.yscale("log")
		plt.ylabel("Median")
		plt.xlabel("Percentage of sparsity")
		plt.title("Filtering results")
		plt.plot([0,100],[median_filter,median_filter], color = "red") #median
		plt.plot([sparsity_filter,sparsity_filter],[0,1500], color = "red") #sparsity
		plt.show()		
	for (c_name, c_data) in new.iteritems():
		zeros_percent = (c_data == 0).astype(int).sum()*100/new.shape[0]
		if zeros_percent > sparsity_filter:
			new = new.drop(c_name, axis=1)
	median = new.median(axis=0)
	median = median < median_filter
	for (c_name, c_data) in new.iteritems():
		if median[c_name]:
			new = new.drop(c_name, axis=1)
	print("filter_prevalence got read of " + str(dataset.shape[1]-new.shape[1]) + " columns over " +str(dataset.shape[1]))
	aro_left_out = [k for k in dataset.columns if k not in new.columns]
	left_out = dataset[aro_left_out]
	print("È stato rimosso " + str(left_out.to_numpy().sum() / dataset.to_numpy().sum()) + " del contenuto dei resistomi")
	X = list(dataset.index)
	Yprima = list(dataset.sum(axis = 1))
	Ydopo = list(new.sum(axis = 1))
	X_axis = np.arange(len(X))
	plt.bar(X_axis - 0.2, Yprima, 0.4, label = 'Before')
	plt.bar(X_axis + 0.2, Ydopo, 0.4, label = 'After')
	plt.xticks(X_axis, X)
	plt.xlabel("Populations")
	plt.ylabel("Spectrum integral for the population")
	plt.legend()
	plt.show()	
	return new

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
	

def preprocessing(dataframe, sparsity=70, median=5):
	"""Creates the database in a format usefull for machine learning
	Paramethers:
		dict_dict: a dictionary containing other dictionary, this argument is supposed to be the returning value from read_bwt
		phenotypes: a dictionary continining the ground throoth about sample's labels, labelled with the sample names
		The standard agruments are passed to the functino filter_median_sparsity, see the comment of it
	Returns:
		A pandas dataframe containing all the components of the data (samples) after various operations are performed: we have
		-Filtering the columns for their sparsity and median vales
		-Normalization with center log rateo
		-Standard scaling
		-PCA
	Each of the steps of this function are independent and can be commented it is required not to have them
	"""
#	view_database(dataset)
#FILTERING FOR PREVALENCE #toglgie quelle con occorrenza di zeri maggiore del tot%
#	dataframe = filter_sparcity_median(dataframe, sparsity, median, plot = False, plot_zoom=15)
#	view_database(dataset, 2)
#NORMALIZE #normalizza con center og rateo
#	dataset = CLR_normalization(dataset)
#MAKING STANDARDIZATION
#	dataset = std_scale(dataset)
#MAKING PCA #il secondo argomento standard determina quante componenti tiene, in questo caso tutte
#	dataset = PCA_maker(dataset, plots=False)
	return dataframe

def loo(classifier, resistome, phenotype, plot = False, linear = False):
	predictions = list(cross_val_predict(classifier, resistome, phenotype, cv=LeaveOneOut(), n_jobs=-1))
	CROSS = cross_validate(classifier, resistome, phenotype, cv=LeaveOneOut(), n_jobs=-1, return_estimator=True)
	#if linear:
	#	feat_imps = [i.coef_[0] for i in models["estimator"]]
	#else:
	#	feat_imps = [i.feature_importances_ for i in models["estimator"]]
	#feat_imps = pd.DataFrame(feat_imps, columns=resistome.columns)
	#feat_imps = {"Feature importance" : feat_imps.mean(axis=0),"Standard deviation": feat_imps.std(axis=0)}
	#if plot:
	#	plt.hist(scores)
	#	plt.xlabel("Scores")
	#	plt.ylabel("Abundance")
	#	plt.show()
	#	X_, Y_ = [i.split(" ")[1] for i in feat_imps["Feature importance"].index], feat_imps["Feature importance"].values
	#	if linear:
	#		Y_ = [abs(iii) for iii in Y_]
	#	plt.bar(X_, Y_)
	#	plt.xlabel('PCA component')
	#	plt.ylabel('Feature importance')
	#	plt.errorbar(X_, Y_, feat_imps["Standard deviation"].values, fmt=' ',color='Black')
	#	plt.show()
	#false_s , false_r = 0,0
	#for a in range(len(phenotype)):
	#	if phenotype[a] == "resistant" and results[a] == "susceptible":
	#		false_s +=1
	#	elif phenotype[a] == "susceptible" and results[a] == "resistant":
	#		false_r +=1
	#print("Leave one out report:")
	#print("False resistant = " + str(false_r))
	#print("False susceptible = " + str(false_s))
	#print("Scores = " + str(scores))
	#print("Score avrage " + str(sum(scores)/len(scores)) + ", variance " + str(np.var(scores)) + '\n')
	return {"Scores":CROSS["test_score"], "Predictions": predictions, "Models": CROSS["estimator"], "Feature importance" : "DA IMPLEMENTARE"}

def cv5(classifier, resistome, phenotype, plot = False, linear = False):
	predictions = list(cross_val_predict(classifier, resistome, phenotype, cv=5, n_jobs=-1))
	CROSS = cross_validate(classifier, resistome, phenotype, cv=5, n_jobs=-1, return_estimator=True)
	#if linear:
	#	feat_imps = [i.coef_[0] for i in models["estimator"]]
	#else:
	#	feat_imps = [i.feature_importances_ for i in models["estimator"]]
	#feat_imps = pd.DataFrame(feat_imps, columns=resistome.columns)
	#feat_imps = {"Feature importance" : feat_imps.mean(axis=0),"Standard deviation": feat_imps.std(axis=0)}
	#if plot:
	#	plt.hist(scores)
	#	plt.xlabel("Scores")
	#	plt.ylabel("Abundance")
	#	plt.show()
	#	X_, Y_ = [i.split(" ")[1] for i in feat_imps["Feature importance"].index], feat_imps["Feature importance"].values
	#	if linear:
	#		Y_ = [abs(iii) for iii in Y_]
	#	plt.bar(X_, Y_)
	#	plt.xlabel('PCA component')
	#	plt.ylabel('Feature importance')
	#	plt.errorbar(X_, Y_, feat_imps["Standard deviation"].values, fmt=' ',color='Black')
	#	plt.show()
	#false_s , false_r = 0,0
	#for a in range(len(phenotype)):
	#	if phenotype[a] == "resistant" and results[a] == "susceptible":
	#		false_s +=1
	#	elif phenotype[a] == "susceptible" and results[a] == "resistant":
	#		false_r +=1
	#print("5 fold report:")
	#print("False resistant = " + str(false_r))
	#print("False susceptible = " + str(false_s))
	#print("Scores = " + str(scores))
	#print("Score avrage " + str(sum(scores)/len(scores)) + ", variance " + str(np.var(scores)) + '\n')
	return {"Scores":CROSS["test_score"], "Predictions and phenotype": [[predictions[i], phenotype[i]] for i in range(len(predictions))], "Models": CROSS["estimator"], "Feature importance" : "DA IMPLEMENTARE"}



pop_1_specie_random_resistance = ["population_0","population_1","population_2","population_3","population_4","population_5","population_6","population_7","population_8","population_9","population_10","population_11","population_12","population_13","population_14","population_15","population_16","population_17","population_18","population_19"]
pop_4_specie_random_resistance = ["population_20","population_21","population_22","population_23","population_24","population_25","population_26","population_27","population_28","population_29","population_30","population_31","population_32","population_33","population_34","population_35","population_36","population_37","population_38","population_39"]
pop_1_specie_si_no_resistance = ["population_40","population_41","population_42","population_43","population_44","population_45","population_46","population_47","population_48","population_49","population_50","population_51","population_52","population_53","population_54","population_55","population_56","population_57","population_58","population_59"]
pop_4_specie_si_no_resistance = ["population_60","population_61","population_62","population_63","population_64","population_65","population_66","population_67","population_68","population_69","population_70","population_71","population_72","population_73","population_74","population_75","population_76","population_77","population_78","population_79"]

Data_bwt = read_bwt_2(pop_1_specie_random_resistance)
print("Data_bwt created")
Data_bwt = preprocessing(Data_bwt)
Data_kmer = read_kmer(pop_1_specie_random_resistance)
print("Data_kmer created")
Metadata = read_metadata(pop_1_specie_random_resistance)

view_database(Data_bwt, "FPKB_sum_per_sample")
view_database(Data_kmer, "Kmer_abundance_per_sample")

for m in Metadata:
	Metadata[m]["Phenotype"].replace(["susceptible", "resistant"], [0,1], inplace=True)
R_abundances = pd.DataFrame.from_dict({pop_name : float(Metadata[pop_name][["Abundance", "Phenotype"]].prod(axis=1).sum(axis=0)) for pop_name in Metadata}, orient="index")
r_abundances = [float(Metadata[pop_name][["Abundance", "Phenotype"]].prod(axis=1).sum(axis=0)) for pop_name in Metadata]

print("\nElastic net linear regression, results in 'El_cv5_result' varible")
El_1 = ElasticNet(alpha=1.0, l1_ratio=0.5)
El_cv5_results_bwt= cv5(El_1, Data_bwt, r_abundances)
El_2 = ElasticNet(alpha=1.0, l1_ratio=0.5, max_iter=1000000, tol=0.005)
El_cv5_results_kmer= cv5(El_2, Data_kmer, r_abundances)

print("\nRandom forest regression, results in 'Forest_cv5_results' variable")
Forest_1 = RandomForestRegressor(random_state = 0)
Forest_cv5_results_bwt = cv5(Forest_1, Data_bwt, r_abundances)
Forest_2 = RandomForestRegressor(random_state = 0)
Forest_cv5_results_kmer= cv5(Forest_2, Data_kmer, r_abundances)

print("\nAda boost regression, results in 'Ada_cv5_results' varible")
Ada = AdaBoostRegressor(n_estimators=100, random_state=0)
Ada_cv5_results_bwt = cv5(Ada, Data_bwt, r_abundances)
Ada_cv5_results_kmer= cv5(Ada, Data_kmer, r_abundances)