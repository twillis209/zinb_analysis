#!/bin/bash
#$ -wd "$HOME"/zinb_analysis/sims/figures
#$ -V
#$ -N figuresPaper_S9
#$ -l h_rt=24:00:00,mem_free=40G,h_vmem=40G
#$ -pe smp 1
#$ -j y
#$ -o "$HOME"/logs/"$JOB_NAME"_"$JOB_ID"_o.txt
Rscript figuresPaper.R 
