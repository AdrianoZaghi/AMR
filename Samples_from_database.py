import pandas as pd
import os
import json
import csv

DATABASE = "Database_campioni.csv"
SAMPLES_DIR = "samples/"
REAL_SAMPLES_LOCATION = "Real_samples_de_Caesar/samples_trim/"
PROCESSORS = 25
FRAGMENT_SIZE = 270
FRAGMENT_STD_DEV = 27

Database = pd.read_csv(DATABASE, index_col=0)
Samples_present = os.listdir(SAMPLES_DIR)
Samples_to_create = [x for x in list(Database.index) if x not in Samples_present]
Untracked_samples = [x for x in Samples_present if x not in list(Database.index)]
if len(Untracked_samples) > 0:
	print("Ci sono dei campioni non tracciati!!!")
	print(Untracked_samples)

print("Ci sarebbero questi campioni non ancora fatti")
print(Samples_to_create)

def create_samples_from_database(Database, Samples_to_create):
	for S in Samples_to_create:
		SEED = Database.loc[S]["Seed"]
		os.mkdir(SAMPLES_DIR + S)
		if Database.loc[S]["Real/Simul"] == "Simul":
			genomes_tot = 0
			with open(SAMPLES_DIR + S + "/input.tsv", mode="w") as Itsv:
				table = csv.writer(Itsv, delimiter = "\t")
				table.writerow(["Genome_ID",	"Species_Name",	"Abundance",	"PATRIC_ID",	"Antibiotic",	"Phenotype",	"Genome_length",	"Novelty_category",	"Genome_file"])
				if Database.loc[S]["Aureus_S"] > 0:
					genomes_tot += 1
					table.writerow(["Aureus_S",	"Staphylococcus_aureus_ATCC_25923",	Database.loc[S]["Aureus_S"],	"1280.99999",	"Oxacillina",	"susceptible",	"2806340",	"known_strain",	"Aureus_S.fasta"])
				if Database.loc[S]["Aureus_R"] > 0:
					genomes_tot += 1
					table.writerow(["Aureus_R",	"Staphylococcus_aureus_ATCC_BAA-2312",	Database.loc[S]["Aureus_R"],	"1280.9999999",	"Oxacillina",	"resistant",	"2825116",	"known_strain",	"Aureus_R.fasta"])
				if Database.loc[S]["Coli"] > 0:
					genomes_tot += 1
					table.writerow(["Coli",	"Escherichia_coli_strain_ATCC_25922",	Database.loc[S]["Coli"],	"562.108487",	"Oxacillina",	"resistant",	"5209985",	"known_strain",	"Coli.fasta"])
				if Database.loc[S]["Pseudomonas"] > 0:
					genomes_tot += 1
					table.writerow(["Pseudomonas",	"Pseudomonas_aeruginosa_ATCC_27853",	Database.loc[S]["Pseudomonas"],	"287.26013",	"Oxacillina",	"resistant",	"6839777",	"known_strain",	"Pseudomonas.fasta"])
				if Database.loc[S]["Arcobacter"] > 0:
					genomes_tot += 1
					table.writerow(["Arcobacter",	"Arcobacter_butzleri_strain_NCTC_12481",	Database.loc[S]["Arcobacter"],	"28197.59",	"Oxacillina",	"susceptible",	"2347607",	"known_strain",	"Arcobacter.fasta"])

			dic = {"simulation_directory": os.getcwd() + "/" + SAMPLES_DIR + S,
			"seed":str(SEED),
			"max_processor":str(PROCESSORS),
			"output_directory":SAMPLES_DIR + S + "/camisim_out",
			"size":str(Database.loc[S]["Size"]),
			"fragm_size_mean":str(FRAGMENT_SIZE),
			"fragm_size_std_dev":str(FRAGMENT_STD_DEV),
			"number_of_samples":"1",
			"equally_distributed_strains":"True",
			"genomes_total":str(genomes_tot),
			"simulated_genomes_per_otu":str(Database.loc[S]["N strains"])}
			with open(SAMPLES_DIR + S + "/input.json", mode="w") as j:
				json.dump(dic,j, indent=4)

		elif Database.loc[S]["Real/Simul"] == "Real":
			if(Database.loc[S]["Size"] == 999):
				print("Salto "+S)
				continue
			Number_of_reads_to_peek = float(Database.loc[S]["Size"]) * 1000000000 / FRAGMENT_SIZE
			if(os.path.isfile(Database.loc[S]["Original sample"]+"_R1_001.fastq.gz")):
				os.system("seqtk sample -s "+str(SEED)+" "+Database.loc[S]["Original sample"]+"_R1_001.fastq.gz"+" "+str(int(Number_of_reads_to_peek))+" > " + SAMPLES_DIR + S + "/reads_1.fastq")
				os.system("seqtk sample -s "+str(SEED)+" "+Database.loc[S]["Original sample"]+"_R2_001.fastq.gz"+" "+str(int(Number_of_reads_to_peek))+" > " + SAMPLES_DIR + S + "/reads_2.fastq")
			else:
				print(Database.loc[S]["Original sample"]+"_R1_001.fastq.gz" + " pèroblema")
				continue
			with open(SAMPLES_DIR + S + "/input.tsv", mode="w") as Itsv:
				table = csv.writer(Itsv, delimiter = "\t")
				table.writerow(["Genome_ID",	"Species_Name",	"Abundance",	"PATRIC_ID",	"Antibiotic",	"Phenotype"])
				if Database.loc[S]["Aureus_S"] > 0:
					table.writerow(["Aureus_S",	"Staphylococcus_aureus_ATCC_25923",	Database.loc[S]["Aureus_S"],	"1280.99999",	"Oxacillina",	"susceptible"])
				if Database.loc[S]["Aureus_R"] > 0:
					table.writerow(["Aureus_R",	"Staphylococcus_aureus_ATCC_BAA-2312",	Database.loc[S]["Aureus_R"],	"1280.9999999",	"Oxacillina",	"resistant"])
				if Database.loc[S]["Coli"] > 0:
					table.writerow(["Coli",	"Escherichia_coli_strain_ATCC_25922",	Database.loc[S]["Coli"],	"562.108487",	"Oxacillina",	"resistant"])
				if Database.loc[S]["Pseudomonas"] > 0:
					table.writerow(["Pseudomonas",	"Pseudomonas_aeruginosa_ATCC_27853",	Database.loc[S]["Pseudomonas"],	"287.26013",	"Oxacillina",	"resistant"])
				if Database.loc[S]["Arcobacter"] > 0:
					table.writerow(["Arcobacter",	"Arcobacter_butzleri_strain_NCTC_12481",	Database.loc[S]["Arcobacter"],	"28197.59",	"Oxacillina",	"susceptible"])

create_samples_from_database(Database, Samples_to_create)

len(Database[Database["Pseudomonas"] > 0].index)