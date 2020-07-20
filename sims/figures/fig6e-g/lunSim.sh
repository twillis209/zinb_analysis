#!/bin/bash
#$ -wd "$HOME"/zinb_analysis/sims/figures/fig6e-g
#$ -V
#$ -l h_rt=24:00:00,h_vmem=20G
#$ -N lunSim
#$ -pe smp 1
#$ -j y
#$ -o "$HOME"/logs/"$JOB_NAME"_"$JOB_ID"_o.txt
Rscript lunSim.R 
