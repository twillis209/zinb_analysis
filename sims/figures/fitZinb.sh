#!/bin/bash
#$ -wd "$HOME"/zinb_analysis/sims/figures
#$ -V
#$ -N fitZinb
#$ -l h_rt=6:00:00,mem_free=8G,h_vmem=12G
#$ -pe smp 8
#$ -e "$HOME"/logs/"$JOB_NAME"_"$JOB_ID"_e.txt
#$ -o "$HOME"/logs/"$JOB_NAME"_"$JOB_ID"_o.txt
Rscript fitZinb.R
