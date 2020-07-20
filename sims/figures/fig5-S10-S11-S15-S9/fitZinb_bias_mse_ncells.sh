#!/bin/bash
#$ -wd "$HOME"/zinb_analysis/sims/figures/fig5-S10-S11-S15-S9
#$ -V
#$ -l h_rt=24:00:00,h_vmem=50G
#$ -N fitZinb_bias_mse_ncells
#$ -pe smp 1
#$ -j y
#$ -t 1-6:1
#$ -tc 6
#$ -o "$HOME"/logs/"$JOB_NAME"_"$JOB_ID"_o.txt
INPUT=("50 Zeisel 1 2 1000" 
	"100 Zeisel 1 2 1000"
	"500 Zeisel 1 2 1000" 
	"1000 Zeisel 1 2 1000" 
	"5000 Zeisel 1 2 1000"	
	"10000 Zeisel 1 2 1000")

eval "Rscript -e $(printf "\"source('fitZinb_bias_mse_ncells.R'); fitSim(nc=%s, ds=\'%s\', b2=%s, offs=%s, eps=%s, BPPARAM=SerialParam())\"" ${INPUT[$SGE_TASK_ID - 1]})"
