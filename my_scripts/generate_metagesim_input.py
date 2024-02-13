import sys
import os
import json
import pandas as pd

def make_abundance(input_tsv, output_path):
	with open(input_tsv, mode="r") as tsv:
		t = pd.read_csv(tsv, sep='\t')
	with open(output_path+"/abundance.tsv", mode="w") as abu:
		for index, row in t.iterrows():
			print(row["Genome_ID"]+"\t"+str(row["Abundance"]), file=abu)

def make_genome_to_id(input_tsv, output_path):
	with open(input_tsv, mode="r") as tsv:
		t = pd.read_csv(tsv, sep='\t')
	with open(output_path+"/genome_to_id.tsv", mode="w") as to_id:
		for index, row in t.iterrows():
			print(row["Genome_ID"]+"\t"+os.getcwd()+"/genomes/"+row["Genome_file"], file=to_id)

def make_metadata(input_tsv, output_path):
	with open(input_tsv, mode="r") as tsv:
		t = pd.read_csv(tsv, sep='\t')
	with open(output_path+"/metadata.tsv", mode="w") as meta:
		print("genome_ID\tOTU\tNCBI_ID\tnovelty_category\tN_strains", file=meta)
		for index, row in t.iterrows():
			print(row["Genome_ID"]+'\t'+str(index+1)+'\t'+str(row["PATRIC_ID"]).split(".")[0]+'\t'+row["Novelty_category"] + '\t' + str(row["N_strains"]),file=meta)

def make_config(input_tsv, input_json, output_path):
	with open(input_tsv, mode="r") as tsv, open(input_json, mode="r") as Json:
		t = pd.read_csv(tsv, sep='\t')
		j = json.load(Json)
	with open(output_path+"/config.ini", mode="w") as conf:
		print("[Main]", file=conf)
		print("seed = "+str(j["seed"]), file=conf)
		print("phase = 0", file=conf)
		print("max_processor = "+str(j["max_processor"]), file=conf)
		print("dataset_id = RL", file=conf)
		print("output_directory = "+j["output_directory"], file=conf)
		print("temp_directory = /tmp", file=conf)
		print("gsa = False", file=conf)
		print("pooled_gsa = False", file=conf)
		print("anonymous = True", file=conf)
		print("compress = 1", file=conf)

		print("\n[ReadSimulator]", file=conf)
		print("readsim = CAMISIM/tools/art_illumina-2.3.6/art_illumina", file=conf)
		print("error_profiles = CAMISIM/tools/art_illumina-2.3.6/profiles", file=conf)
		print("samtools = CAMISIM/tools/samtools-1.3/samtools", file=conf)
		print("profile = mbarc", file=conf)
		print("size = "+str(j["size"]), file=conf)
		print("type = art", file=conf)
		print("fragments_size_mean = "+str(j["fragm_size_mean"]), file=conf)
		print("fragment_size_standard_deviation = "+str(j["fragm_size_std_dev"]), file=conf)

		print("\n[CommunityDesign]", file=conf)
		print("ncbi_taxdump = CAMISIM/tools/ncbi-taxonomy_20170222.tar.gz", file=conf)
		print("strain_simulation_template = CAMISIM/scripts/StrainSimulationWrapper/sgEvolver/simulation_dir", file=conf)
		print("number_of_samples = 1", file=conf)

		print("\n[community0]", file=conf)
		print("metadata = "+j["simulation_directory"]+"/metadata.tsv", file=conf)
		print("id_to_genome_file = "+j["simulation_directory"]+"/genome_to_id.tsv", file=conf)
		print("id_to_gff_file = ", file=conf)
		print("path_to_abundance_file = "+j["simulation_directory"]+"/abundance.tsv", file=conf)
		print("genomes_total = " + j["genomes_total"], file=conf)
		print("num_real_genomes = " + str(t.shape[0]), file=conf)
#		print("max_strains_per_otu = " + j["simulated_genomes_per_otu"], file=conf)
		print("ratio = 1", file=conf)
		print("mode = known_distribution", file=conf)
#		print("equally_distributed_strains = True", file=conf)
		print("strain_style = strain_from_metadata", file=conf)
		print("input_genomes_to_zero = True", file=conf)
		print("log_mu = 1", file=conf)
		print("log_sigma = 2", file=conf)
		print("gauss_mu = 1", file=conf)
		print("gauss_sigma = 1", file=conf)
		print("view = False", file=conf)


input_tsv = sys.argv[1]
input_json = sys.argv[2]
output_path = input_tsv.split("/input.tsv")[0]

make_metadata(input_tsv, output_path)
make_abundance(input_tsv, output_path)
make_genome_to_id(input_tsv, output_path)
make_config(input_tsv, input_json, output_path)