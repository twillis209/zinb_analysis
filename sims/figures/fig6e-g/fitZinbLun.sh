#!/bin/bash
#$ -wd "$HOME"/zinb_analysis/sims/figures/fig6e-g
#$ -V
#$ -l h_rt=48:00:00,h_vmem=30G
#$ -N fitZinbLun
#$ -pe smp 1
#$ -t 1-3:1
#$ -tc 3
#$ -j y
#$ -o "$HOME"/logs/"$JOB_NAME"_"$JOB_ID"_o.txt

#INPUT=("100 ziadd0"
#	"100 ziadd0.33"
#	"100 ziadd0.67"
#	"1000 ziadd0"
#	"1000 ziadd0.33"
#	"1000 ziadd0.67"
#	"10000 ziadd0"
#	"10000 ziadd0.33"
#	"10000 ziadd0.67")

INPUT=("1000 ziadd0.67" 
	"10000 ziadd0"
	"10000 ziadd0.33"
	"10000 ziadd0.67")
eval "Rscript -e $(printf "\"source('fitZinbLun.R'); fitSimData(ncells=%s, add=\'%s\', BPPARAM=SerialParam())\"" ${INPUT[$SGE_TASK_ID - 1]})"
