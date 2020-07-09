#!/bin/bash
#$ -wd "$HOME"/zinb_analysis/sims/figures/fig5-S10-S11-S15-S9
#$ -V
#$ -l h_rt=12:00:00,h_vmem=12G
#$ -N timeZinb_multicoreParam
#$ -pe smp 6
#$ -e "$HOME"/logs/"$JOB_NAME"_e.txt
#$ -o "$HOME"/logs/"$JOB_NAME"_o.txt
Rscript timeZinb.R 
