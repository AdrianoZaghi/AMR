suppressMessages(library("knitr"))
suppressMessages(library("BiocStyle"))
suppressMessages(library("ggplot2")) 
suppressMessages(library(gridExtra))
suppressMessages(library(phyloseq))
suppressMessages(library(DECIPHER))
suppressMessages(library(phangorn))
suppressMessages(library(RcppParallel))
suppressMessages(library(parallel))
suppressMessages(library(R.utils))

original_samples_directory <- "/home/PERSONALE/adriano.zaghi2/AMR/De_Caesar/samples_left"
trimmed_dir_path <- "/home/PERSONALE/adriano.zaghi2/AMR/De_Caesar/samples_trim"
seq_suffix <- "001" #999 usa i due campioni di prova, 001, 002 e 003 usa tutti i campioni delle run 3, 4 e 5
images_path="images/"
plots=TRUE

CUTADAPT = FALSE

cutadapt_adapter_F <- "CTGTCTCTTATA"
cutadapt_adapter_R <- "CTGTCTCTTATA"
trimmomatic_adapters_path <- "Trimmomatic/adapters/NexteraPE-PE.fa"
trash <- "trimmomatic_trash_files/"

#CTGTCTCTTATACACATCT


MinLen<-50		#MINIMA LUNGHEZZA PER NON SCARTARE UNa READ DOPO IL TRIMMING
############## ATTENZIONE: DADA2 REQUIRES NO Ns   ############################
MaxN <- 0		#MASSIMO NUMERO DI N TROVATE NELLA SEQUENZA PER NON SCARTARLA (proporzione)
##############
MaxErrors<-4 	#MASSIMO NUMERO DI ERRORI DA ASPETTARSI NELLA SEQUENZA PER NON RIGETTARLA (SONO STIMATI CON I PUNTEGGI DI QUALITà DELLE BASI)
Len_F <- 10		#Remove LEN bases from each read (or R1 if paired; use -U option for R2). If LEN is positive, remove bases from the beginning. If LEN is negative, remove bases from the end. Can be used twice if LENs have different signs. Applied *before* adapter trimming.
Len_R <- 10

search_files <- function(myseq_path, seq_suffix=""){
	#-------------------myseq_path
	#Path to dir. containing all sequence files
	#-------------------seq_suffix
	#The function will consider just those sequences that match that specific path
	#-------------------return a list of columns:
	#	"F"/"R": path to forward/reverse read files
	#	"Sample names": nomi dei campioni
	#-------------------Note
	#File names must match be in the following format: SampleName_R1_Whatever.fastq and SampleName_R2_Whatever.fastq 
	fnFs <- sort(list.files(myseq_path, pattern=paste0("R1_",seq_suffix,".fastq"), full.names=TRUE))
	fnRs <- sort(list.files(myseq_path, pattern=paste0("R2_",seq_suffix,".fastq"), full.names=TRUE))
	if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")
	Names <- sapply(strsplit(fnFs, "_R1_"), "[", 1)
	Names <- sapply(strsplit(Names, "/"), function(x) tail(x,n=1))
	r <- list(fnFs, fnRs, Names)
	names(r) <- c("F", "R", "Sample names")
	return(r)
}

generate_trimmed_path <- function(original_path, trimmed_dir_path){
	#Routine function for the function trim
	split <- unlist(strsplit(original_path, "/"))
	return(paste0(trimmed_dir_path, "/", split[length(split)]))
}

report_trimming <- function(Original_F_paths, Original_R_paths, Trimmed_F_paths, Trimmed_R_paths){
#NUMBER OF TOTAL INITIAL READS
	Original_F <- unlist(lapply(Original_F_paths, function(x) countLines(x)/4))
#NUMBER OF TOTAL TRIMMED SEQUENCE
	trim_F <- unlist(lapply(Trimmed_F_paths, function(x) countLines(x)/4))
	data <- data.frame(cbind(Original_F,trim_F))
	data$Rimaste <- data$trim_F / data$Original_F
	return(data)
}

print(paste(Sys.time(),"START"))
original_samples <- search_files(myseq_path=original_samples_directory, seq_suffix = "001")
print(paste(Sys.time(),"FILES READED"))

F_out <- sort(unlist(lapply(original_samples[["F"]], generate_trimmed_path, trimmed_dir_path)))
R_out <- sort(unlist(lapply(original_samples[["R"]], generate_trimmed_path, trimmed_dir_path)))
F_in <- sort(original_samples[["F"]])
R_in <- sort(original_samples[["R"]])
for(i in 1:length(F_out)){
		if(!file.exists(F_out[[i]])){
			if(CUTADAPT){
				command <- paste0("cutadapt --cores=0 -g ","\"",cutadapt_adapter_F,";max_errors=0.1\""," -a ","\"",cutadapt_adapter_R,";max_errors=0.1\""," -o ",F_out[[i]]," -p ",R_out[[i]]," ", F_in[[i]]," ", R_in[[i]]," --minimum-length ",MinLen," --max-n ",MaxN," --max-ee ",MaxErrors," --trim-n"," --poly-a")
			}else{
				command <- paste0("java -jar Trimmomatic/trimmomatic-0.39.jar PE -threads 25 ", F_in[[i]], " ", R_in[[i]], " ", F_out[[i]], " ", trash,"F_",i,".fastq ",R_out[[i]]," ",trash,"R_",i,".fastq ", "ILLUMINACLIP:",trimmomatic_adapters_path,":2:30:10:2:True SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:", MinLen)
				print(command)
			}
			suppressMessages(system(command))
		}
}
print(paste(Sys.time(),"SAMPLES TRIMMED"))
#print(report_trimming(F_in, R_in, F_out, R_out))
