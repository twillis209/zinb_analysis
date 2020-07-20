#!/bin/bash
#$ -wd "$HOME"/zinb_analysis/svenssonApproach/R
#$ -V
#$ -N fitZinbToKolodziejczyk
#$ -l h_rt=24:00:00,mem_free=40G,h_vmem=40G
#$ -pe smp 1
#$ -j y
#$ -o "$HOME"/logs/"$JOB_NAME"_"$JOB_ID"_o.txt
Rscript fitZinbToKolodziejczyk.R 
