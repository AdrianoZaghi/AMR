#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 17:34:05 2022

@author: utente
"""
import json


def generate_input_tsv(Strains_path, Metadata_path, Resistance_percentage, Other_genomes_path = []):
	species_space = 1 / (len(Other_genomes_path)+1)
	print(species_space)
	with open("input.tsv", 'w') as tab:
		with open(Metadata_path+"_S.json", 'r') as S, open(Metadata_path+"_R.json", 'r') as R: 
			print("Genome_ID\tSpecies_Name\tAbundance\tPATRIC_ID\tAntibiotic\tPhenotype\tGenome_length\tNovelty_category\tGenome_file", file = tab)
			Su, Re = json.load(S), json.load(R)
			print(Su["Genome_ID"]+"\t"+Su["Species_Name"]+"\t"+str(round(species_space*(1 - Resistance_percentage),5))+"\t"+Su["PATRIC_ID"]+"\t"+Su["Antibiotic"]+"\t"+Su["Phenotype"]+"\t"+Su["Genome_length"]+"\t"+Su["Novelty_category"]+"\t"+Metadata_path+"_S.json", file=tab)
			print(Re["Genome_ID"]+"\t"+Re["Species_Name"]+"\t"+str(species_space*Resistance_percentage)+"\t"+Re["PATRIC_ID"]+"\t"+Re["Antibiotic"]+"\t"+Re["Phenotype"]+"\t"+Re["Genome_length"]+"\t"+Re["Novelty_category"]+"\t"+Metadata_path+"_R.json", file=tab)
		for path in Other_genomes_path:
			with open(path, 'r') as O:
				Ot = json.load(O)
				print(Ot["Genome_ID"]+"\t"+Ot["Species_Name"]+"\t"+str(species_space*Resistance_percentage)+"\t"+Ot["PATRIC_ID"]+"\t"+Ot["AntibioticPhenotype"]+"\t"+Ot["Genome_length"]+"\t"+Ot["Novelty_category"]+"\t"+path, file=tab)

pop = "population_prova_0"
generate_input_tsv("Roberto", pop+"/genomes/Metadata_Aureus", 0.7)