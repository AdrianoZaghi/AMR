assembly_K_list = "27,47,67,87,107,127"
rgi_annotation_alignment_tool = "BLAST"
k_mer_size = 27
anno = "2023"
max_count = 120000000
wild_card_version = "4.0.0"
card_version = "3.2.6"
import os

Populations = os.listdir("samples/",)

rule terminated:
	input:
				alignments_test = expand("samples/{pop}/paired_alignment", pop=Populations),
# 				kmer_counts = expand("{pop}/single_kmer_count", pop=Populations),
#				single_end_reads = expand("samples/{pop}/reads.fastq", pop=Populations)
#				pop_configurations = expand("samples/{pop}/config.ini", pop=Populations)
	output:
		"TERMINATED"
	shell:
		"""
		date
		touch TERMINATED
		"""

rule population_set_up:
	input:
		i_1 = "{sample}/input.tsv",
		i_2 = "{sample}/input.json"
	output:
		"{sample}/config.ini"
	shell:
		"""
		python my_scripts/generate_camisim_input.py {input.i_1} {input.i_2}
		"""

rule camisim:
	input:
		C = "{sample}/config.ini"
	output:
		"{sample}/reads_1.fastq",
		"{sample}/reads_2.fastq",
		"{sample}/reads.fastq"
	conda: "env.yaml"
	threads: 6
	resources:
		mem_mb = 100000
	shell:
		"""
		rm -rf {wildcards.sample}/camisim_out
		mkdir -p tmp
		python CAMISIM/metagenomesimulation.py {input.C}
		gunzip {wildcards.sample}/camisim_out/{anno}*/reads/anonymous_reads.fq.gz
		python my_scripts/split_pair_end.py {wildcards.sample}/camisim_out/{anno}*/reads/anonymous_reads.fq {wildcards.sample}/
		mv {wildcards.sample}/camisim_out/{anno}*/reads/anonymous_reads.fq {wildcards.sample}/reads.fastq
		rm -r {wildcards.sample}/camisim_out/{anno}*/reads
		"""

rule kmer_count_paired:
	input:
		S1 = "{sample}/reads_1.fastq",
		S2 = "{sample}/reads_2.fastq"
	output:
		"{sample}/kmc_out/kmer_count_1.txt",
		"{sample}/kmc_out/kmer_count_2.txt",
		"{sample}/paired_kmer_count"
	conda: "env.yaml"
	threads: 20
	params:
		ram = 10
	shell:
		"""
		rm -r {wildcards.sample}/kmc_out
		mkdir {wildcards.sample}/kmc_out
		mkdir {wildcards.sample}/kmc_out/kmc_temp
		kmc -k{k_mer_size} -m{params.ram} -t{threads} -cs{max_count} -cx{max_count} {input.S1} {wildcards.sample}/kmc_out/k_1 {wildcards.sample}/kmc_out/kmc_temp
		kmc_tools transform {wildcards.sample}/kmc_out/k_1 dump {wildcards.sample}/kmc_out/kmer_count_1.txt
		
		kmc -k{k_mer_size} -m{params.ram} -t{threads} -cs{max_count} -cx{max_count} {input.S2} {wildcards.sample}/kmc_out/k_2 {wildcards.sample}/kmc_out/kmc_temp
		kmc_tools transform {wildcards.sample}/kmc_out/k_2 dump {wildcards.sample}/kmc_out/kmer_count_2.txt
		touch {wildcards.sample}/paired_kmer_count
		"""

rule kmer_count_single:
	input:
		S = "{sample}/reads.fastq"
	output:
		"{sample}/kmc_out/kmer_count.txt",
		"{sample}/single_kmer_count"
	conda: "env.yaml"
	threads: 20
	params:
		ram = 10
	shell:
		"""
		rm -r {wildcards.sample}/kmc_out
		mkdir {wildcards.sample}/kmc_out
		mkdir {wildcards.sample}/kmc_out/kmc_temp
		kmc -k{k_mer_size} -m{params.ram} -t{threads} -cs{max_count} -cx{max_count} {input.S} {wildcards.sample}/kmc_out/K {wildcards.sample}/kmc_out/kmc_temp
		kmc_tools transform {wildcards.sample}/kmc_out/K dump {wildcards.sample}/kmc_out/kmer_count.txt
		touch {wildcards.sample}/single_kmer_count
		"""

rule rgi_bowtie2_paired:
	input:
		S1 = "{sample}/reads_1.fastq",
		S2 = "{sample}/reads_2.fastq",
		r1 = "rgi/card_annotation.log",
		r2 = "rgi/wildcard_annotation.log"
	output:
		"{sample}/paired_alignment",
		"{sample}/bwt_out/bwt.allele_mapping_data.json"
	conda: "env.yaml"
	threads: 5
	shell:
		"""
		mkdir -p {wildcards.sample}/bwt_out
		cd rgi
		rgi bwt -1 ../{input.S1} -2 ../{input.S2} -a bowtie2 -n {threads} --clean --include_wildcard -o ../{wildcards.sample}/bwt_out/bwt --local 
		rm ../{wildcards.sample}/bwt_out/bwt.sorted.length_100.bam*
		cd ..
		touch {wildcards.sample}/paired_alignment
		"""

rule rgi_bowtie2_paired_compressed:
	input:
		S1 = "{sample}/reads_1.fastq.gz",
		S2 = "{sample}/reads_2.fastq.gz",
		r1 = "rgi/card_annotation.log",
		r2 = "rgi/wildcard_annotation.log"
	output:
		"{sample}/paired_alignment",
		"{sample}/bwt_out/bwt.allele_mapping_data.json"
	conda: "env.yaml"
	threads: 5
	shell:
		"""
		mkdir -p {wildcards.sample}/bwt_out
		cd rgi
		rgi bwt -1 ../{input.S1} -2 ../{input.S2} -a bowtie2 -n {threads} --clean --include_wildcard -o ../{wildcards.sample}/bwt_out/bwt --local 
		rm ../{wildcards.sample}/bwt_out/bwt.sorted.length_100.bam*
		cd ..
		touch {wildcards.sample}/paired_alignment
		"""
		
rule rgi_bowtie2_single:
	input:
		S = "{sample}/reads.fastq",
		r1 = "rgi/card_annotation.log",
		r2 = "rgi/wildcard_annotation.log"
	output:
		"{sample}/single_alignment",
		"{sample}/bwt_out/bwt.allele_mapping_data.json"
	conda: "env.yaml"
	threads: 20
	shell:
		"""
		mkdir -p {wildcards.sample}/bwt_out
		cd rgi
		rgi bwt -1 ../{input.S} -a bowtie2 -n {threads} --clean --include_wildcard -o ../{wildcards.sample}/bwt_out/bwt --local
		rm ../{wildcards.sample}/bwt_out/bwt.sorted.length_100.bam*
		cd ..
		touch {wildcards.sample}/single_alignment
		"""

rule rgi_env_setup:
	output:
		"rgi/card_annotation.log",
		"rgi/wildcard_annotation.log"
	conda: "env.yaml"
	shell:
		"""
		rm -r rgi
		mkdir rgi
		cd rgi
		#CARD DATABASE
		wget -O data https://card.mcmaster.ca/latest/data
		tar -xvf data ./card.json
		#rgi load --card_json card.json --local
		rgi card_annotation -i card.json > card_annotation.log 2>&1
		#rgi load -i card.json --card_annotation card_database_v{card_version}.fasta --local
		#WILDCARD DATABASE
		wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
		mkdir -p wildcard
		tar -xjf wildcard_data.tar.bz2 -C wildcard
		gunzip wildcard/*.gz
		rgi wildcard_annotation -i wildcard --card_json card.json -v{wild_card_version} > wildcard_annotation.log 2>&1
		rgi load --card_json card.json --debug --local --card_annotation card_database_v{card_version}.fasta --wildcard_annotation wildcard_database_v{wild_card_version}.fasta --wildcard_index wildcard/index-for-model-sequences.txt --wildcard_version {wild_card_version} --amr_kmers wildcard/all_amr_61mers.txt --kmer_database wildcard/61_kmer_db.json --kmer_size 61
		#rgi load --wildcard_annotation wildcard_database_v{wild_card_version}.fasta --wildcard_index wildcard/index-for-model-sequences.txt --card_annotation card_database_v{card_version}.fasta --local
		cd ..
		"""