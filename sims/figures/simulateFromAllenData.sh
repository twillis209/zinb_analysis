#!/bin/bash
#$ -wd "$HOME"/zinb_analysis/sims/figures
#$ -V
#$ -N simulateFromAllenData
#$ -l h_rt=24:00:00,mem_free=20G,h_vmem=20G
#$ -pe smp 1
#$ -j y
#$ -o "$HOME"/logs/"$JOB_NAME"_"$JOB_ID"_o.txt
Rscript -e "source('simFunction.R'); simulateFromAllenData(BPPARAM=SerialParam());"
