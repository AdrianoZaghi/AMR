assembly_K_list = "27,47,67,87,107,127"
rgi_annotation_alignment_tool = "BLAST"
k_mer_size = 27
anno = "2022"
max_count = 120000000

rule population_set_up:
	input:
		i_1 = "{sample}/input.json",
		i_2 = "{sample}/input.tsv",
		g = "{sample}/genomes"
	output:
		"{sample}/config.ini"
	shell:
		"""
		python my_scripts/generate_files.py {wildcards.sample}
		"""

rule camisim_paired:
	input:
		C = "{sample}/config.ini"
	output:
		"{sample}/reads_1.fastq",
		"{sample}/reads_2.fastq"
	conda: "env.yaml"
	threads: 20
	shell:
		"""
		rm -rf {wildcards.sample}/camisim_out
		mkdir -p {wildcards.sample}/tmp
		python CAMISIM/metagenomesimulation.py {input.C}
		gunzip {wildcards.sample}/camisim_out/{anno}*/reads/anonymous_reads.fq.gz
		python my_scripts/split_pair_end.py {wildcards.sample}/camisim_out/2022*/reads/anonymous_reads.fq {wildcards.sample}/
		"""

rule camisim_single:
	input:
		C = "{sample}/config.ini"
	output:
		"{sample}/reads.fastq"
	conda: "env.yaml"
	threads: 20
	shell:
		"""
		rm -rf {wildcards.sample}/camisim_out
		mkdir -p {wildcards.sample}/tmp
		python CAMISIM/metagenomesimulation.py {input.C}
		gunzip {wildcards.sample}/camisim_out/{anno}*/reads/anonymous_reads.fq.gz
		mv {wildcards.sample}/camisim_out/{anno}*/reads/anonymous_reads.fq {wildcards.sample}/reads.fastq
		"""

rule kmer_count_paired:
	input:
		S1 = "{sample}/reads_1.fastq",
		S2 = "{sample}/reads_2.fastq"
	output:
		"{sample}/kmc_out/kmer_count_1.txt.gz",
		"{sample}/kmc_out/kmer_count_2.txt.gz"
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
		gzip {wildcards.sample}/kmc_out/kmer_count_1.txt
		kmc -k{k_mer_size} -m{params.ram} -t{threads} -cs{max_count} -cx{max_count} {input.S2} {wildcards.sample}/kmc_out/k_2 {wildcards.sample}/kmc_out/kmc_temp
		kmc_tools transform {wildcards.sample}/kmc_out/k_2 dump {wildcards.sample}/kmc_out/kmer_count_2.txt
		gzip {wildcards.sample}/kmc_out/kmer_count_2.txt
		"""

rule kmer_count_single:
	input:
		S = "{sample}/reads.fastq"
	output:
		"{sample}/kmc_out/kmer_count.txt.gz"
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
		gzip {wildcards.sample}/kmc_out/kmer_count_1.txt
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
	threads: 20
	shell:
		"""
		mkdir -p {wildcards.sample}/bwt_out
		cd rgi
		rgi bwt -1 ../{input.S1} -2 ../{input.S2} -a bowtie2 -n {threads} --clean --include_wildcard -o ../{wildcards.sample}/bwt_out/bwt --local 
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
		mkdir -p rgi
		cd rgi
		#CARD DATABASE
		wget -O data https://card.mcmaster.ca/latest/data
		tar -xvf data ./card.json
		rgi load --card_json card.json --local
		rgi card_annotation -i card.json > card_annotation.log 2>&1
		rgi load -i card.json --card_annotation card_database* --local
		#WILDCARD DATABASE
		wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
		mkdir -p wildcard
		tar -xjf wildcard_data.tar.bz2 -C wildcard
		gunzip wildcard/*.gz
		rgi wildcard_annotation -i wildcard --card_json card.json -v_non_specificata > wildcard_annotation.log 2>&1
		rgi load --wildcard_annotation wildcard_database* --wildcard_index wildcard/index-for-model-sequences.txt --card_annotation card_database* --local
		cd ..
		"""
