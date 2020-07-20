#!/bin/bash
#$ -wd "$HOME"/zinb_analysis/sims/figures/fig5-S10-S11-S15-S9
#$ -V
#$ -l h_rt=24:00:00,h_vmem=20G
#$ -N timeZinb
#$ -pe smp 1
#$ -j y
#$ -o "$HOME"/logs/"$JOB_NAME"_"$JOB_ID"_o.txt
Rscript timeZinb.R 
