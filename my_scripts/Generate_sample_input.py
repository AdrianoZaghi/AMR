#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 17:34:05 2022

@author: utente
"""
import json
from numpy import random

def generate_input_tsv(input_name, Strain_name, Metadata_path, Resistance_percentage, Other_genomes_path = []):
	species_space = 1 / (len(Other_genomes_path)+1)
	with open(input_name, 'w') as tab:
		with open(Metadata_path+"_S.json", 'r') as S, open(Metadata_path+"_R.json", 'r') as R: 
			print("Genome_ID\tSpecies_Name\tAbundance\tPATRIC_ID\tAntibiotic\tPhenotype\tGenome_length\tNovelty_category\tGenome_file", file = tab)
			Su, Re = json.load(S), json.load(R)
			print(Su["Genome_ID"]+"\t"+Su["Species_Name"]+"\t"+str(round(species_space*(1 - Resistance_percentage),5))+"\t"+Su["PATRIC_ID"]+"\t"+Su["Antibiotic"]+"\t"+Su["Phenotype"]+"\t"+Su["Genome_length"]+"\t"+Su["Novelty_category"]+"\t"+Strain_name+"_S.fasta", file=tab)
			print(Re["Genome_ID"]+"\t"+Re["Species_Name"]+"\t"+str(species_space*Resistance_percentage)+"\t"+Re["PATRIC_ID"]+"\t"+Re["Antibiotic"]+"\t"+Re["Phenotype"]+"\t"+Re["Genome_length"]+"\t"+Re["Novelty_category"]+"\t"+Strain_name+"_R.fasta", file=tab)
		for names in Other_genomes_path:
			with open(names, 'r') as O:
				Ot = json.load(O)
				print(Ot["Genome_ID"]+"\t"+Ot["Species_Name"]+"\t"+str(species_space)+"\t"+Ot["PATRIC_ID"]+"\t"+Ot["Antibiotic"]+"\t"+Ot["Phenotype"]+"\t"+Ot["Genome_length"]+"\t"+Ot["Novelty_category"]+"\t"+names.split("Metadata_")[1].split(".")[0]+".fasta", file=tab)

how_many_samples = 20
distribuzione_random = random.uniform(low=0, high=1, size=how_many_samples)
distribuzione_binaria = [1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]
distribuzione_de_cesare = [0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.5,0.5,0.5,0.75,0.75,0.75,0.75,0.75,1,1,1,1,1]
for i in range(how_many_samples):
#	generate_input_tsv("input_p"+str(i), "Aureus","population_20/genomes/Metadata_Aureus", round(distribuzione_binaria[i],5), ["population_20/genomes/Metadata_Escherichia.json","population_20/genomes/Metadata_Acineto.json", "population_20/genomes/Metadata_Pseudomosa.json"])
	generate_input_tsv("input_p"+str(i), "Aureus","population_20/genomes/Metadata_Aureus", round(distribuzione_de_cesare[i],5))
