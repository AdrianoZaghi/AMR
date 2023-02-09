#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 17:34:05 2022

@author: utente
"""
import json
import os

def generate_input_tsv(Filename, abundance, metadatas, are_percentages=True):
	with open(Filename, mode='w') as tab:
		print("Genome_ID\tSpecies_Name\tAbundance\tPATRIC_ID\tAntibiotic\tPhenotype\tGenome_length\tNovelty_category\tGenome_file", file = tab)
		for i in range(len(abundance)):
			Ab=abundance[i]
			if (not are_percentages) and Ab>0:			
				with open(metadatas[i], 'r') as O:
					Ot = json.load(O)
					print(Ot["Genome_ID"]+"\t"+Ot["Species_Name"]+"\t"+str(Ab)+"\t"+Ot["PATRIC_ID"]+"\t"+Ot["Antibiotic"]+"\t"+Ot["Phenotype"]+"\t"+Ot["Genome_length"]+"\t"+"known_strain"+"\t"+Ot["Genome_file"], file=tab)
			elif Ab>0:
				with open(metadatas[i], 'r') as O:
					Ot = json.load(O)
					print(Ot["Genome_ID"]+"\t"+Ot["Species_Name"]+"\t"+str(round(Ab*0.01, 14))+"\t"+Ot["PATRIC_ID"]+"\t"+Ot["Antibiotic"]+"\t"+Ot["Phenotype"]+"\t"+Ot["Genome_length"]+"\t"+"known_strain"+"\t"+Ot["Genome_file"], file=tab)

def generate_input_json(Filename, abundance, processors=8, size=0.01, fragment_size=270, fragment_dev=27):
	dic = {"simulation_directory":os.getcwd()+"/"+Filename.split("/input.json")[0],
	"seed":"1",
	"max_processor":str(processors),
	"output_directory":Filename.split("/input.json")[0]+"/camisim_out",
	"size":str(size),
	"fragm_size_mean":str(fragment_size),
	"fragm_size_std_dev":str(fragment_dev),
	"number_of_samples":"1",
	"equally_distributed_strains":"True",
	"genomes_total":str(sum([a>0 for a in abundance]))}
	with open(Filename, mode="w") as j:
		json.dump(dic,j, indent=4)

Metadata_paths = ["genomes/Aureus_S_metadata.json", "genomes/Aureus_R_metadata.json", "genomes/Coli_metadata.json", "genomes/Pseudomonas_metadata.json", "genomes/Arcobacter_metadata.json"]


Train_0=[[27,73,0,0,0],[86,14,0,0,0],[92,8,0,0,0],[38,62,0,0,0],[37,63,0,0,0],\
		[14,86,0,0,0],[80,20,0,0,0],[28,72,0,0,0],[32,68,0,0,0],[75,25,0,0,0]]

Test_0=[[95,5,0,0,0],[95,5,0,0,0],[99,1,0,0,0],[99,1,0,0,0],[99.5,0.5,0,0,0],\
			[99.5,0.5,0,0,0],[95,0,0,5,0],[95,0,0,5,0],[99,0,0,1,0],[99,0,0,1,0]]

Train_1=[[28,44,0,0,28],[19,62,0,0,19],[3.5,93,0,0,3.5],[39,22,0,0,39],[11,78,0,0,11],\
		[39.5,21,0,0,39.5],[3,94,0,0,3],[9,82,0,0,9],[9.5,81,0,0,9.5],[36.5,27,0,0,36.5]]

Test_1=[[47.5,5,0,0,47.5],[47.5,5,0,0,47.5],[49.5,1,0,0,49.5],[49.5,1,0,0,49.5],[49.75,0.5,0,0,49.75],\
		[49.75,0.5,0,0,49.75],[47.5,0,5,0,47.5],[47.5,0,5,0,47.5],[49.5,0,1,0,49.5],[49.5,0,1,0,49.5]]

Train_2=[[3,47,47,0,3],[26,24,24,0,26],[41,9,9,0,41],[32,18,18,0,32],[31,19,19,0,31],\
		[40,10,10,0,40],[48,2,2,0,48],[23,27,27,0,23],[15,35,35,0,15],[25,25,25,0,25]]

Test_2=[[47.5,5,0,0,47.5],[49.5,1,0,0,49.5],[49.75,0.5,0,0,49.75],[47.5,0,5,0,47.5],[49.5,0,1,0,49.5],[49.75,0,0.5,0,49.75],\
		[47.5,0,0,5,47.5],[49.5,0,0,1,49.5],[49.75,0,0,0.5,49.75],[47.5,2.5,2.5,0,47.5],[49.5,0.5,0.5,0,49.5],[47.5,2.5,0,2.5,47.5],\
		[49.5,0.5,0,0.5,49.5],[47.5,0,2.5,2.5,47.5],[49.5,0,0.5,0.5,49.5]]
directory_tag = "Train_0"
abundances=Train_0
for i in range(len(abundances)):
	path=directory_tag+str(i).zfill(2)
	if not os.path.exists(path):
		os.mkdir(path)
	generate_input_tsv(path+"/input.tsv", abundances[i], Metadata_paths, are_percentages=True)
	generate_input_json(path+"/input.json", abundances[i])