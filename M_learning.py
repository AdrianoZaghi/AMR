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

def regression(Train, Test, Train_target, classifier, parameter_space):
	#GESTISCE LA GRIDSEARCH PER ALGORITMI DI REGRESSIONE 
	missing_in_Train = [c for c in Test.columns if c not in Train.columns]
	missing_in_Test = [c for c in Train.columns if c not in Test.columns]
	Test[missing_in_Test] = 0
	Train[missing_in_Train] = 0
	Test=Test.reindex(sorted(Test.columns), axis=1)
	Train=Train.reindex(sorted(Train.columns), axis=1)
	grid = GridSearchCV(classifier, parameter_space, n_jobs=-1, refit=True).fit(Train, Train_target)
	best = grid.best_estimator_
	Test_predictions = best.predict(Test)
	#print("Target \t Predictions")
	#print(Test_target, Test_predictions)
	return({"Best estimator":best, "Predictions":Test_predictions})

campioni = pd.read_csv("Database_campioni.csv", index_col=0)
WHICH = "Tutti gli aureus"
#Modelli allenati su tutti i campioni sintetici (0.01Gbp) (1 0 5 strain) che contengono esclusivamente Aureus S o R.
#Testati su tutti i campioni reali (0.01Gbp) che contengono esclusivamnte Aureus S o R
#7/11/23		R2 = 0.9584
#9/11/23		R2 = 0.9404
#WHICH = "Aureus 1 strain train"
#Modelli allenati su i campioni sintetici (0.01Gbp) da 1 sola strain che contengono esclusivamente Aureus S o R.
#Testati su tutti i campioni reali (0.01Gbp) che contengono esclusivamnte Aureus S o R
#7/11/23		R2 = 0.4297
#WHICH = "Aureus 5 strain train"
#Modelli allenati su i campioni sintetici (0.01Gbp) da 5 strain che contengono esclusivamente Aureus S o R.
#Testati su tutti i campioni reali (0.01Gbp) che contengono esclusivamnte Aureus S o R
#7/11/23		R2 = 0.9622
#WHICH = "Tutti gli aureus e Pseudomonas"
#Modelli allenati su i campioni sintetici (0.01Gbp) (1 0 5 strain) che contengono esclusivamente Aureus S o R e Pseudomonas.
#Testati su tutti i campioni reali (0.01Gbp) che contengono esclusivamnte Aureus S o R e Pseudomonas.
#Fiasco
#Da provare sintetici su sintetici
#WHICH = "Tutti gli aureus e Coli"
#Fiasco
#Da provare sintetici su sintetici
#WHICH = "0.01 GB Aureus e Arcobacter"
#Fiasco
#Da provare sintetici su sintetici
#WHICH = "0.1 GB Aureus Arcobacter"
#Fiasco
#Sintetici su sintetici funziona
WHICH = "0.1 GB Aureus Pseudomonas"
WHICH = "Confornto base"

print(WHICH)
if(WHICH == "Tutti gli aureus"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] == 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
elif(WHICH == "Aureus 1 strain train"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [1,0] for x in t["N strains"]]]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] == 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
elif(WHICH == "Aureus 5 strain train"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[[x in [5,0] for x in t["N strains"]]]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] == 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]	
elif(WHICH == "Tutti gli aureus e Pseudomonas"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[t["Coli"] == 0]
	t = t[t["Arcobacter"] == 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
elif(WHICH == "Tutti gli aureus e Coli"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] == 0]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
elif(WHICH == "0.01 GB Aureus e Arcobacter"):
	t = campioni[campioni["Size"] == 0.01]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	t = t[t["N strains"] < 7]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in t[t["Real/Simul"] == "Simul"].index]
elif(WHICH == "0.1 GB Aureus Arcobacter"):
	t = campioni[campioni["Size"] == 0.1]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] == 0]
	t = t[t["Arcobacter"] != 0]
	t = t[t["N strains"] < 7]
	Tr = t
	#Tr = t[["tatuaggio" in x for x in t.index]]
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in Tr[Tr["Real/Simul"] == "Simul"].index]
elif(WHICH == "0.1 GB Aureus Pseudomonas"):
	t = campioni[campioni["Size"] == 0.1]
	t = t[t["Coli"] == 0]
	t = t[t["Pseudomonas"] != 0]
	t = t[t["Arcobacter"] == 0]
	t = t[t["N strains"] < 7]
	Tr = t
	Test = ["samples/" + x for x in t[t["Real/Simul"] == "Real"].index]
	Train = ["samples/" + x for x in Tr[Tr["Real/Simul"] == "Simul"].index]
elif(WHICH == "Confornto base"):
	Train = ["samples/Patric_" + x for x in ["Aureus_R", "Aureus_S", "Coli", "Pseudomonas", "Arcobacter"]]
	Test = ["samples/Decesare_" + x for x in ["Aureus_R", "Aureus_S", "Coli", "Pseudomonas", "Arcobacter"]]

print("Quanti trains: " + str(len(Train)))
print("Quanti tests: " + str(len(Test)))

El = ElasticNet(max_iter=100000, tol=0.05)
params_El = {"alpha":[0.01]+list(range(1, 100, 10)), "l1_ratio":[0, .1, .5, .7, .9, .95, .99, 1]}

Forest = RandomForestRegressor(random_state = 0)
params_Forest = {"criterion":["squared_error", "absolute_error", "friedman_mse", "poisson"]}

Ada =  AdaBoostRegressor(n_estimators=1000, random_state=0)
params_Ada = {"loss":["linear", "square", "exponential"], "learning_rate":[0.5, 1, 1.5, 2 ,3, 5,10]}

#Info_controlli = read_bwt_2_various(Train)
#im.plot_controls(Info_controlli["Antibiotics"], images_dir+"Control_antibiotics.pdf")
#im.plot_controls(Info_controlli["Drugs_class"], images_dir+"Control_drugs.pdf")
#im.plot_controls(Info_controlli["Pathogens"], images_dir+"Control_pathogens.pdf")
#[Info_controlli[c].to_csv("Controllo_"+c+".tsv", sep="\t") for c in Info_controlli]

def log_ods(vector):
	return(np.log([vv/(1 - vv) for vv in vector]))

def log_shift(vector):
	return(np.log([v+1 for v in vector]))

Train_data = read_bwt_2_aro(Train)
Test_data = read_bwt_2_aro(Test)
#Train_data = Train_data.drop(columns=["3000178", "3007014", "3002835"])
#Test_data = Test_data.drop(columns=["3000178", "3007014"  , "3002835"])
Train_target= log_shift(read_metadata(Train)["Aureus_r"])
Test_target= log_shift(read_metadata(Test)["Aureus_r"])

#Data = read_bwt_2_aro(Train)
#Target = log_shift(read_metadata(Train)["Aureus_r"])
#Train_data, Test_data, Train_target, Test_target = train_test_split(Data, Target, test_size=0.3, random_state=42)

#Results = regression(Train_data, Test_data, Train_target, Test_target, El, params_El)
#Results = regression(Train_data, Test_data, Train_target, Test_target, Forest, params_Forest)
Results = regression(Train_data, Test_data, Train_target, Ada, params_Ada)
R2 = r2_score(Test_target, Results["Predictions"])
im.plot_results(Test_target, Results["Predictions"], "M_learning_images/prova", R2 = R2)

#plt.hist(read_metadata(Train)["Aureus_r"], bins=50)  # Choose an appropriate number of bins
if (True):
	Train_d = Train_data
	#Train_t = Train_target
	Train_d = Train_d.sort_index(axis=1, ascending=False)
	#Questo non serve per i confronti base
	#new_index_order = pd.DataFrame(Train_t, index=Train_d.index).sort_values(0).index
	#Train_d = Train_d.reindex(new_index_order)
	#Cols = [col for col in Train_d.columns if Train_d[col].mean() > 0]
	plt.figure(figsize=(15, 6))
	sns.heatmap(Train_d, cmap='viridis')

if (True):
	Saturation_LV = 10000
	Test_d = Test_data
	Test_d = Test_d.where(Test_d < Saturation_LV, Saturation_LV)
	#Test_t = Test_target
	Test_d = Test_d.sort_index(axis=1, ascending=False)
	#Cols = [col for col in Test_d.columns if Test_d[col].mean() > 0]
	plt.figure(figsize=(15, 6))
	sns.heatmap(Test_d, cmap='viridis')

#prova a fare train e test con campioni simulati