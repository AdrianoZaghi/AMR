#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 17:34:05 2022

@author: utente
"""
import json
import os

def generate_input_tsv(Filename, abundance, metadatas):
	with open(Filename, 'w') as tab:
		print("Genome_ID\tSpecies_Name\tAbundance\tPATRIC_ID\tAntibiotic\tPhenotype\tGenome_length\tNovelty_category\tGenome_file", file = tab)
		for i in range(len(abundance)):
			Ab=abundance[i]
			if Ab < 1 and Ab > 0:			
				with open(metadatas[i], 'r') as O:
					Ot = json.load(O)
					print(Ot["Genome_ID"]+"\t"+Ot["Species_Name"]+"\t"+str(Ab)+"\t"+Ot["PATRIC_ID"]+"\t"+Ot["Antibiotic"]+"\t"+Ot["Phenotype"]+"\t"+Ot["Genome_length"]+"\t"+"known_strain"+"\t"+Ot["Genome_file"], file=tab)
			elif Ab > 1:
				with open(metadatas[i], 'r') as O:
					Ot = json.load(O)
					print(Ot["Genome_ID"]+"\t"+Ot["Species_Name"]+"\t"+str(Ab*0.01)+"\t"+Ot["PATRIC_ID"]+"\t"+Ot["Antibiotic"]+"\t"+Ot["Phenotype"]+"\t"+Ot["Genome_length"]+"\t"+"known_strain"+"\t"+Ot["Genome_file"], file=tab)

def generate_input_json(Filename, abundance, processors=8, size=0.01, fragment_size=270, fragment_dev=27):
	dic = {"simulation_directory":Filename.split("/input.json")[0],
	"seed":"1",
	"max_processor":str(processors),
	"output_directory":Filename.split("/input.json")[0]+"camisim_out",
	"size":str(size),
	"fragm_size_mean":str(fragment_size),
	"fragm_size_std_dev":str(fragment_dev),
	"number_of_samples":"1",
	"equally_distributed_strains":"True",
	"genomes_total":str(sum([a>0 for a in abundance]))}
	with open(Filename, mode="w") as j:
		json.dump(dic,j)

Metadata_paths = ["genomes/Aureus_S_metadata.json", "genomes/Aureus_R_metadata.json", "genomes/Coli_metadata.json", "genomes/Pseudomonas_metadata.json", "genomes/Archeobacter_metadata.json"]
directory_tag = "Train_0"
abundances=[[27,73,0,0,0],\
			[86,14,0,0,0],\
			[92,8,0,0,0],\
			[38,62,0,0,0],\
			[37,63,0,0,0],\
			[14,86,0,0,0],\
			[80,20,0,0,0],\
			[28,72,0,0,0],\
			[32,68,0,0,0],\
			[75,25,0,0,0]]

for i in range(len(abundances)):
	path=directory_tag+str(i).zfill(2)
	os.mkdir(path)
	generate_input_tsv(path+"/input.tsv", abundances[i], Metadata_paths)
	generate_input_json(path+"/input.json", abundances[i])