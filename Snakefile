assembly_K_list = "27,47,67,87,107,127"
rgi_annotation_alignment_tool = "BLAST"
k_mer_size = 27
anno = "2025"
max_count = 120000000
wild_card_version = "4.0.0"
card_version = "3.2.6"
import os

Populations2 = os.listdir("samples_2/",)
Populations3 = os.listdir("samples_3/",)
Populations4 = os.listdir("samples_4/",)


rule terminated_2:
	input:
				alignments2 = expand("samples_2/{pop}/paired_alignment", pop=Populations2)
	output:
		"TERMINATED_2"
	shell:
		"""
		date
		touch TERMINATED_2
		"""

rule terminated_eliminate_reads_2:
	input:
				alignments2 = expand("samples_2/{pop}/paired_alignment_eliminate_reads", pop=Populations2)
	output:
		"TERMINATED_no_reads_2"
	shell:
		"""
		date
		touch TERMINATED_no_reads_2
		"""

rule terminated_3:
	input:
				alignments3 = expand("samples_3/{pop}/paired_alignment", pop=Populations3)
	output:
		"TERMINATED_3"
	shell:
		"""
		date
		touch TERMINATED_3
		"""

rule terminated_eliminate_reads_4:
	input:
				alignments3 = expand("samples_4/{pop}/paired_alignment_eliminate_reads", pop=Populations4)
	output:
		"TERMINATED_no_reads_4"
	shell:
		"""
		date
		touch TERMINATED_no_reads_4
		"""

rule eliminate_reads:
	input:
		i_1 = "{sample}/paired_alignment"
	output:
		"{sample}/paired_alignment_eliminate_reads"
	shell:
		"""
		rm -f {input.i_1}
		rm -f {wildcards.sample}/reads_1.fastq
		rm -f {wildcards.sample}/reads_2.fastq
		rm -f {wildcards.sample}/reads.fastq
		touch {wildcards.sample}/paired_alignment_eliminate_reads
		"""

genomes_dir = "Consensus_genomes_assembled"
rule population_set_up_and_Metagesim:
	input:
		per_distinguere_questa_regola_da_camisim = "{sample}/USE_METAGESIM",
		i_1 = "{sample}/input.tsv",
		i_2 = "{sample}/input.json",
	output:
		"{sample}/reads_1.fastq",
		"{sample}/reads_2.fastq",
		"{sample}/reads.fastq"
	conda: "env.yaml"
	threads: 4
	resources:
		mem_mb = 100000
	shell:
		"""
		python my_scripts/generate_metagesim_input.py {input.i_1} {input.i_2} {genomes_dir}
		rm -rf {wildcards.sample}/camisim_out
		mkdir -p tmp
		python MetaGeSim-AMR/metagenomesimulation2.py {wildcards.sample}/config.ini
		gunzip {wildcards.sample}/camisim_out/{anno}*/reads/anonymous_reads.fq.gz
		python my_scripts/split_pair_end.py {wildcards.sample}/camisim_out/{anno}*/reads/anonymous_reads.fq {wildcards.sample}/
		mv {wildcards.sample}/camisim_out/{anno}*/reads/anonymous_reads.fq {wildcards.sample}/reads.fastq
		rm -r {wildcards.sample}/camisim_out/{anno}*/reads
		rm -r {wildcards.sample}/camisim_out/{anno}*/bam/*.bam
		rm -r {wildcards.sample}/camisim_out/{anno}*/contigs
		rm -r {wildcards.sample}/camisim_out/source_genomes
		"""

rule population_set_up:
	input:
		per_distinguere_questa_regola_da_metagesim = "{sample}/USE_CAMISIM",
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
		per_distinguere_questa_regola_da_metagesim = "{sample}/USE_CAMISIM",
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
	threads: 4
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
	threads: 4
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

rule Kraken:
	input:
		S1 = "Real_samples_de_Caesar/samples_trim/{sample}_R1_001.fastq.gz",
		S2 = "Real_samples_de_Caesar/samples_trim/{sample}_R2_001.fastq.gz",
	output:
		"kraken2/Kraken_results/{sample}"
	conda: "env.yaml"
	threads: 20
	shell:
		"""
		./kraken2/kraken2 --db kraken2/My_reference --threads {threads} --paired --output kraken2/Kraken_results/{wildcards.sample} --report kraken2/Kraken_results/{wildcards.sample}_report.txt {input.S1} {input.S2}
		"""

rule Kraken_report_for_me:
	input:
		S="kraken2/Kraken_results/{sample}"
	output:
		"kraken2/Kraken_results/{sample}_report.csv"
	shell:
		"""
		python my_scripts/kraken_report.py {input.S} -o kraken2/Kraken_results/{wildcards.sample}_report.csv
		rm {input.S}
		"""