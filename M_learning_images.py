from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
from sklearn.linear_model import ElasticNet
from operator import itemgetter
import json
import seaborn as sns
import pandas as pd
average_read_length = 270

def plot_results(target, prediction, filename = None, R2 = 666):
	plt.plot(target, prediction, 'ro')
	x=[max([min(target)-0.01, min(prediction)-0.01]), min([max(target)+0.01,max(prediction)+0.01])]
	y=x
	plt.plot(x, y, 'k:', lw=2)
	plt.xlabel("Target")
	plt.ylabel("Predictions")
	plt.title("Target vs Prediction: R2 = " + str(round(R2,4)))
	if not filename == None:
		plt.savefig(filename)

def plot_controls(dataf, filename, n_columns=5):
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
	Return = Return.reindex(popolazioni)
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

def log_ods(vector):
	#	Ci metto un -50 per evitare che il log di 0 dia errore
	ods = [vv/(1 - vv) for vv in vector]
	log_ods = [np.log(vv) if vv < 20 else np.log(20,) for vv in ods]
	return([vv if vv > -4 else -4 for vv in log_ods])

def log_shift(vector):
	return(np.log([v+1 for v in vector]))

def normalize_counts(data):
	r = data.apply(lambda x: x / x.max())

def feature_importance(Train_data, fitted_model, filename = None):
	if(isinstance(fitted_model,ElasticNet)):
		F_imp = fitted_model.coef_
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

def heatmap_abundance(Data, Target = None, F_imp = None, Normalize = False, saturation = -1,filename = None):
	if Target is not None:
		sorted_rows = dict(sorted(Target.items(), key=lambda item: item[1]))
		d = Data.reindex(sorted_rows)
	if F_imp is not None:
		sorted_columns = dict(sorted(F_imp.items(), key=lambda item: item[1]))
		d = d[sorted_columns]
	if saturation > 0:
		Data = Data.where(Data < saturation, saturation)
	if Normalize:
		d_norm = d.apply(lambda x: x / x.max())
		plt.figure(figsize=(15, 6))
		sns.heatmap(d_norm, cmap='viridis')
	else:
		plt.figure(figsize=(15, 6))
		sns.heatmap(d, cmap='viridis')