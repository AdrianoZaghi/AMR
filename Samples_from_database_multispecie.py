import pandas as pd
import os
import json
import csv

###########################################
#DA RUNNARE NELL'AMBIENTE DI BASE, CHE CONTIENE SEQTK

DATABASE = "Database_campioni_multispecie.csv"
SAMPLES_DIR = "samples_3/"
PROCESSORS = 25
FRAGMENT_SIZE = 255
FRAGMENT_STD_DEV = 27
genomes = [
        "Aureus_R",
        "Aureus_S",
		"Coli",
		"Pseudomonas",
		"Arcobacter",
        "Bacteroides_1_1_6",
        "Bacteroides_4_3_47",
        "Bacteroides_capillosus",
        "Bacteroides_coprocola",
        "Bacteroides_dorei",
        "Bacteroides_fragilis",
        "Bacteroides_plebeius",
        "Bacteroides_thetaiotaomicron",
        "Bacteroides_uniformis",
        "Bacteroides_vulgatus",
        "Bifidobacterium_adolescentis",
        "Bifidobacterium_longum",
        "Bleuta_hansenii",
        "Clostridium_asparagiforme",
        "Clostridium_leptum",
        "Clostridium_phytofermentans",
        "Clostridium_saccharolyticum",
        "Ethanoligenens_harbinense",
        "Faecalibacterium_prausnitzii",
        "Holdemania_filiformis",
        "Marvinbryantia_formatexigens",
        "Ruminococcaceae_bacterium",
        "Ruminococcus_albus",
        "Ruminococcus_torques",
        "Subdoligranulum_variabile"]

Database = pd.read_csv(DATABASE, index_col=0)
Samples_present = os.listdir(SAMPLES_DIR)
Samples_to_create = [x for x in list(Database.index) if x not in Samples_present]
Untracked_samples = [x for x in Samples_present if x not in list(Database.index)]
if len(Untracked_samples) > 0:
	print("Ci sono dei campioni non tracciati!!!")
	print(Untracked_samples)
print("Ci sarebbero "+str(len(Samples_to_create))+" campioni non ancora fatti")
print(Samples_to_create)

METADATA = {}
for G in genomes:
	with open("genomes/" + G + "_metadata.json") as j:
		dic = json.load(j)
		METADATA[G] = dic
METADATA = pd.DataFrame(METADATA).T	

def create_Simul_samples_from_database(Database, Samples_to_create, SAMPLES_DIR, genomes):
	for S in Samples_to_create:
		SEED = Database.loc[S]["Seed"]
		os.mkdir(SAMPLES_DIR + S)
		genomes_tot = 0
		with open(SAMPLES_DIR + S + "/input.tsv", mode="w") as Itsv:
			table = csv.writer(Itsv, delimiter = "\t")
			table.writerow(["Genome_ID",	"Species_Name",	"Abundance", "N_strains"	,"PATRIC_ID",	"Antibiotic",	"Phenotype",	"Genome_length",	"Novelty_category",	"Genome_file"])
			for G in genomes:
				if Database.loc[S][G] > 0:
					genomes_tot += int(Database.loc[S]["N_strains_for_all"])
					table.writerow([G, METADATA["Species_Name"][G], Database.loc[S][G], Database.loc[S]["N_strains_for_all"],METADATA["PATRIC_ID"][G],METADATA["Antibiotic"][G],METADATA["Phenotype"][G],METADATA["Genome_length"][G],"known_strain",METADATA["Genome_file"][G]])
				
		dic = {"simulation_directory": os.getcwd() + "/" + SAMPLES_DIR + S,
		"seed":str(SEED),
		"max_processor":str(PROCESSORS),
		"output_directory":SAMPLES_DIR + S + "/camisim_out",
		"size":str(Database.loc[S]["Size"]),
		"fragm_size_mean":str(Database.loc[S]["reads_length"]),
		"fragm_size_std_dev":str(FRAGMENT_STD_DEV),
		"number_of_samples":"1",
		"genomes_total":str(genomes_tot)}
		if Database.loc[S]["Real/Simul"] in ["Simul_new_error", "Simul_perfected"]:
			dic["error_profile"] = "own"
			dic["base_profile_name"] = "De_Cesare_error_profilesR"
			dic["profile_read_length"] = 251
		with open(SAMPLES_DIR + S + "/input.json", mode="w") as j:
			json.dump(dic,j, indent=4)
		with open(SAMPLES_DIR + S + "/USE_METAGESIM", mode="w") as j:
			pass

if False:
	print("I'm just whatching")
else:
#Nuove funzioni parallelizzate per i campioni reali
	create_Simul_samples_from_database(Database = Database, Samples_to_create = Samples_to_create, SAMPLES_DIR = SAMPLES_DIR, genomes = genomes)