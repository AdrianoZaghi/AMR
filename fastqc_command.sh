#!/bin/bash

fastqc De_Caesar/samples_trim/* -o fastqc_out/ --noextract --format fastq --threads 25 --dir fastqc_temp_dir
