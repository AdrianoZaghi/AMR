from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
from sklearn.linear_model import ElasticNet
from operator import itemgetter
import random
import json
import seaborn as sns
import pandas as pd
import numpy as np
import os
average_read_length = 270

def plot_results(target, prediction, filename = None, R2 = None, x_up_lim = None, x_low_lim = None, y_up_lim = None, y_low_lim = None):
	"""Plots the target vs the prediction of a model
	The is a scatter plot: in an optimal case, the points should be on the line y=x
	Args:
		target (list): the target values
		prediction (list): the predicted values
		filename (str): the name of the file where the plot will be saved
		R2 (float): the R2 score of the model
	
	"""
	plt.plot(target, prediction, 'ro')
	x=[max([min(target)-0.01, min(prediction)-0.01]), min([max(target)+0.01,max(prediction)+0.01])]
	y=x
	plt.plot(x, y, 'k:', lw=2)
	if x_up_lim is not None or x_low_lim is not None:
		plt.xlim(x_low_lim, x_up_lim)
	if y_up_lim is not None or y_low_lim is not None:
		plt.ylim(x_low_lim, x_up_lim)
	plt.xlabel("Target")
	plt.ylabel("Predictions")
	title = "Target vs Prediction"
	if R2 is not None:
		title += ": R2 = " + str(round(R2,4))
	plt.title(title)
	if not filename == None:
		plt.savefig(filename)

def plot_controls(dataf, filename, n_columns=5):
	"""Non mi ricordo cosa fa bene, ma va usata con i campioni con una sola pecie per vedere le features dei geni e se coincidono con il batterio nel campione
	"""
	trans=dataf.T
	pdf = matplotlib.backends.backend_pdf.PdfPages(filename)
	for c in trans.columns:
		control = dict(trans[c])
		control_biggest = dict(sorted(control.items(), key = itemgetter(1), reverse = True)[:n_columns])
		fig,ax = plt.subplots()  #create a new figure
		ax.bar(list(control_biggest.keys()), list(control_biggest.values()))
		ax.set_title(c)
		ax.tick_params(labelrotation=12, labelsize=10)
		fig.savefig(pdf, format='pdf')
	pdf.close()
	return(control_biggest)

def read_bwt_2_aro(popolazioni, reduce_depth = 1):
	"""Reads the abundances of the ARO entries detected by the rgi bwt alignment
	The ARO entries are identifiers of the resistomes present in the rgi
	This function returns a database with the estimated abundances of each resistome for each sample
	Args:
		popolazioni (list): the list of the populations paths
	"""
	if reduce_depth <= 0 or reduce_depth > 1:
		raise ValueError("reduce_depth must be a number between 0 and 1, bigger than 0")
	resistomes = {}
	for P in popolazioni:
		with open(P+"/bwt_out/bwt.allele_mapping_data.json", mode = "r") as f:
			dati_ = json.load(f)
			dati = [di for di in dati_ if isinstance(di, dict)]
			R={}
			for aro in dati:
				aro_accession = str(aro["aro_accession"])
				aro_count = int(aro["reads"]["all"])
				if reduce_depth < 1:
					aro_count = np.random.binomial(n=aro_count, p=reduce_depth)
				if aro_accession in R.keys():
					R[aro_accession] += aro_count * average_read_length / int(aro["reference"]["sequence_length"])
				else:
					R[aro_accession] = aro_count * average_read_length / int(aro["reference"]["sequence_length"])
				
			resistomes[P] = R
	Return = pd.DataFrame.from_dict(resistomes, orient="index").fillna(0)
	Return = Return.reindex(popolazioni)
	return Return

def read_bwt_2_various(popolazioni):
	"""This function returns dataframeS containing features associated to the resistomes
	The abundances of these features are the sum of the abundances of the resistomes having that feature in the sample
	The features are stored in a dictionaty with the following keys:
		"Pathogens" Pathogen associated to the resistome
		"Antibiotics" Antibiotic for which resistance is developed
		"Drug_class" Drug for which resistance is developed
	Args:
		popolazioni (list): the list of the populations paths"""
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
	"""Reads the metadata of the populations
	The metadata is the file input.tsv in the populations folder
	Returns a dictionary with:
		"Metadata" A dataframe for each population containing all the metadata
		"Abundances_r" The abundance of the resistant phenotype for each population
		"Aureus_r" The abundance of the Aureus_R genome for each population
	Args:
		popolazioni (list): the list of the populations paths"""
	metadata = {}
	for P in popolazioni:
		M_pop = pd.read_csv(P+"/input.tsv", sep="\t")
		metadata[P] = M_pop
	[metadata[m]["Phenotype"].replace(["susceptible", "resistant"], [0,1], inplace=True) for m in metadata]
	r_abundances = {}
	Aureus_r_abundance = {}
	for pop_name in metadata:
		M = metadata[pop_name]
		tot = sum(M["Abundance"])
		r_abundances[pop_name] = M[["Abundance", "Phenotype"]].prod(axis=1).sum(axis=0)/tot
		if("Aureus_R" in M.values):
			Aureus_r_abundance[pop_name] = M[M["Genome_ID"]=="Aureus_R"].iloc[0]["Abundance"]/tot
		else: Aureus_r_abundance[pop_name] = 0
	return {"Metadata":metadata, "Abundances_r":r_abundances, "Aureus_r": Aureus_r_abundance}

def read_abundances_from_tax_ettore(popolazioni,  correspondance_SubSamp:[None, str]= None, correspondance_genomes = "Taxonomy_samples_de_cesare/correspondance_genomes.json", verbose=False):
	"""Reads the abundances of the genomes detected by the tax_ettore
	This function returns a database with the estimated abundances of each genome for each sample
	Args:
		popolazioni (list): the list of the populations paths
		correspondance_SubSamp (str) = Json file containing correspondance between original De Cesare samples and it's subsamples
			if None, the samples which contain the name of the original De Cesare sample will be considered to be subsamle of them
		correspondance_genomes (str) = "Taxonomy_samples_de_cesare/correspondance_genomes.json": the path of the file containing the correspondance between the genomes names and the Machine learining labels
	"""
	df_return = {}
	with open(correspondance_genomes, mode = "r") as f:
		info = json.load(f)
		correspondance = info["Contig_to_genome"]
		G_len = info["Genome_length"]
	if correspondance_SubSamp is not None:
		print("Questa versione è ancora da implementare, per ora non mi serve")
		with open(correspondance_SubSamp, mode = "r") as f:
			correspondance_SubSamp = json.load(f)
	else:
		ettore_tax = os.listdir("Taxonomy_samples_de_cesare")
		ettore_tax = [et for et in ettore_tax if "txt" in et]
		correspondance_SubSamp = {}
		for P in popolazioni:
			abundances = None
			for et in ettore_tax:
				if et.split("_table")[0] in P:
					#print("Correspondance between "+P.split("/")[1]+" and "+et)
					correspondance_SubSamp[P.split("/")[1]] = et
					abundances = pd.read_csv("Taxonomy_samples_de_cesare/"+et, sep="\t", index_col=0)
					break
			if abundances is None:
				if verbose: print("No correspondance for "+P)
			else:
				bacteria = {k:0 for k in set(correspondance.values())}
				for index, row in abundances.iterrows():
					bacteria[correspondance[index]] += row["mapped_reads"]
					if verbose: print(bacteria)
				df_return[P] = bacteria
	df_return = pd.DataFrame.from_dict(df_return, orient="index")
	df_return = df_return.div(G_len, axis=1)
	return df_return.div(df_return.sum(axis=1), axis=0)


def regression(Train, Test, Train_target, classifier, parameter_space, Keep_columns_that_appear_just_in_one_set = False):
	"""Performs a regression on the Train data and predicts the Test data
	Returns a dictionary with:
		"Best estimator": the best model found by the grid search
		"Predictions": the predictions of the model on the Test data
	Args:
		Train (DataFrame): the data used for training
		Test (DataFrame): the data used for testing
		Train_target (list): the target values for the training
		classifier (sklearn model): the model used for the regression, IT NEEDS TO HAVE THE .fit() METHOD
		parameter_space (dict): the hyperparameters of the model
		Keep_columns_that_appear_just_in_one_set (bool) = False: if True, the columns that appear only in Train or Test will be kept and filled with 0, if False, they will be dropped
	""" 
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

def log_odds(dict, upper_saturation = 0.99999, lower_saturation = 0.00001):
	#	Ci metto un -50 per evitare che il log di 0 dia errore
	adjust_boundaries  = {k : 0.99999 if v == 1 else v for k, v in dict.items()}
	adjust_boundaries  = {k : 0.00001 if v == 0 else v for k, v in adjust_boundaries.items()}
	ods = {k : v/(1 - v) for k, v in adjust_boundaries.items()}
	return {k : np.log(v) for k, v in ods.items()}
	

def log_shift(dict, shift = 1):
	return {k: np.log(v + shift) for k, v in dict.items()}

def normalize_counts(data):
	r = data.apply(lambda x: x / x.max())

def logistic(dict, min = 0.001, max = 0.999, scale = 1):
	d = {k: min if v < min else v for k, v in dict.items()}
	d = {k: max if v > max else v for k, v in d.items()}
	return {k : -np.log((1/v)-1)*scale for k, v in d.items()}

def feature_importance(Train_data, fitted_model, filename = None):
	"""Plots the feature importance of a model and returns a dictionary with the feature importance
	At the moment, it works only with the ElasticNet model and the ensambles models of sklearn: RandomForest, AdaBoost....
	Args:
		Train_data (DataFrame): the data used for training
		fitted_model (sklearn model): the model used for the regression, IT NEEDS TO HAVE THE .fit() METHOD
		filename (str) = None: the name of the file where the plot will be saved"""
	if(isinstance(fitted_model,ElasticNet)):
		F_imp = [abs(i) for i in fitted_model.coef_]
	else:
		F_imp = fitted_model.feature_importances_
	F = Train_data.columns
	F_imp = {F[i]:F_imp[i] for i in range(len(F))}
	F_imp = dict(sorted(F_imp.items(), key=lambda item: item[1]))
	F_imp_head = {i : F_imp[i] for i in list(F_imp.keys())[-15:]}
	plt.figure(figsize=(15, 6))
	plt.bar(F_imp_head.keys(), F_imp_head.values())
	if not filename == None:
		plt.savefig(filename)
	return F_imp

def heatmap_abundance(Data, Target = None, F_imp = None, Normalize = False, saturation = -1,rempove_empthy_col = True, filename = None, figsize = (15, 6)):
	"""Plots a heatmap of the data
	If Target is not None, the rows will be sorted by the target values
	If F_imp is not None, the columns will be sorted by the feature importance
	If Normalize is True, the values in each column will be normalized by the maximum value in the column
	If saturation is not >0, the values higher than saturation level will be set to saturation so that lower values are more visible
	Args:
		Data (DataFrame): the data used for the heatmap
		Target (list) = None: the target values
		F_imp (list) = None: the feature importance
		Normalize (bool) = False: if True, the data will be normalized
		saturation (int) = -1: the maximum value of the data
		filename (str) = None: the name of the file where the plot will be saved
	"""
	d = Data
	if Target is not None:
		sorted_rows = dict(sorted(Target.items(), key=lambda item: item[1]))
		d = d.reindex(sorted_rows)
	if F_imp is not None:
		sorted_columns = dict(sorted(F_imp.items(), key=lambda item: item[1]))
		d = d[sorted_columns]
	if saturation > 0:
		d = d.where(d < saturation, saturation)
	if rempove_empthy_col:
		print("Removing"+str(len(d.columns[(d == 0).all()]))+"empty columns")
		d = d.loc[:, (d != 0).any(axis=0)]
		if Normalize:
			d = d.apply(lambda x: x / x.max())
	elif Normalize:
		print("You can't normalize the data without removing the empty columns, please set rempove_empthy_col = True or Normalize = False")
	plt.figure(figsize=figsize)
	sns.heatmap(d, cmap='viridis')
	if not filename == None:
		plt.savefig(filename)
	return d