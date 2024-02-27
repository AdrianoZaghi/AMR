#!/bin/bash

fastqc Real_samples_de_Caesar/samples_trim/DCWGS1_* -o fastqc_out/ --noextract --format fastq --threads 25 --dir fastqc_temp_dir
