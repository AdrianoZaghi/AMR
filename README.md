# PYPELINE TO GENERATE SYNTHETIC SAMPLES AND PERFORM ML ABOUT AMR

The Snakefile you can find in this directory was used to generate synthetic bulk DNA. Those are genrated in a classical way, without the emplyment of any ML generative model, but starging from reference genomes. Those synthetic samples have been used to train ML models that was tested on real data. The target of this inspection is quantify the AMR-bacteria abundance in the sample.
For the samples generation, the pipeline relies on the software MetaGeSim https://github.com/EttoreRocchi/MetaGeSim-AMR (and its older version CAMISIM https://github.com/CAMI-challenge/CAMISIM)
Note that the MetaGeSim version have been slightly modified for this project and the updates still needs to be pushed.
To explain the functioning of the piepeline is usefull to devide the process in three:

- Read library simulation
- Allignemnt
- Machine learning

## Read library generation

The read library generation is the procedure that we implemented in order to obtain simulated read libraries starting from genome sequences (.fasta files). This job is handled by MetaGeSim. The command used is reported in the Snakefile under the rule "camisim". In this rule are generated the reads in the anonimus mode of the software and then is used a script in order to separe the backward and forward reads in the two pair end files.
In order this command to work, you need to set up a population directory. In this folder will be stored all the data produced about that population. To automate the setting up, we implemented the rule "population set up". In this other rule is lunched a script that generates the input files required by MetaGeSim, starting from the information contained in input.json and input.tsv. If further automation is needed, the base informations can be stored in a row of "Database Campioni" and, if the folder is not alreadt present, a folder containing the input files is created lunching in python the script "Samples_form_database.py"
Following this procedure, many paramethers of MetaGeSim execution are set with default values that we found appropriate for our goal. Anyway, if someone intend to modify those parameters should modify the scripts "my_scripts/generate_files.py"  or "Samples_from_database.py" and follow camisim manual.

Once the directories are fully prepared, the following command can be used (in the AMR directory) to generate the read libraries:
Paired end
```
snakemake -c8 population_folder/reads_1.fastq --use-conda
```
Single end
```
snakemake -c8 population_folder/reads adriano
.fastq --use-conda
```

If you have real DNA bulk samples that should be involved in the ML process, you can subasample them so that they have the same size of the simulated ones. This is implemented in the script "Samples_from_database.py" and requires to adapt some of the columns values. Notice that, even if the samples are real, the script and the database format expect the sample composition to be declared, so that the metadata can be created. 

## Allignemnt

The allignment is the procedure that associates the read libraries to the references from CARD and WildCARD. This procedure relyes on the software rgi in its bwt mode. The command used to run the allignment is reported in the rule "rgi_bowtie2" and requires as input the read libraries resulting from the grneration part. Rgi requires to set up an envoirment that include the borrows weeler transform of the reference database. This procedure can be performed just once for all, usually if is the firt time running the pipeline. If it is required, the rule "rgi_bowtie2" calls the rule "rgi_env_setup" that contains the rgi commands used in order to prepare this envoirment.

To obtain the perform the allignmet, the following command can be used (in the AMR directory) to generate the read libraries:
```
snakemake -c8 population_*/bwt_out/bwt.allele_mapping_data.json --use-conda
```

## Rule all
If all the population folders are properly set up, you can call the output
```
snakemake -c8 TERMINATED --use-conda
```
To both symulate the samples and perform the alignement

## Machine learning

This part of the pipeline is not included in the Snakefile, but is handeld by the script "M_learning.py" and its dependencies. 
In this scripts, the population

The snakefiles contains also other rules, which are not already functioning, so do not realy on them.
