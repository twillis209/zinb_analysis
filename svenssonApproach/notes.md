# Plan

The plan is to run Svenssons's analysis on the (real) data sets of Risso et al. 2018. Not sure about the simulated data yet.

I will start by running the code from `Gene-wise_dispersion.ipynb`. Not sure this actually needs to be in a notebook. 

Svensson looked at:

* negative controls
* homogeneous cells
* heterogeneous cells

Which data sets are suitable for this investigation? 

* All contain ERCC spike-ins
* Homogeneous cells:
	* Kolodziejczyk has cells separated by batch and plate; separating only by plate led to a very poor fit
	* Allen contains cells from the V1 visual cortex
	* Fletcher has cells separated by lineage
	* Zeisel has cells from S1 somatosensory cortex and the CA1 cell field of the hippocampus
	* A pity that we don't have access to the Patel counts as cells from each tumour would be suitably homogeneous (one would hope)

# Risso data formats

Svensson code is written for `h5ad` format, but the format of the Risso data sets varies.

* Allen: 
	* now in `csv` file
	* contains ERCC features, need filtering
	* header and rownames
* Fletcher:
	* in tab-separated `txt` file
	* contains ERCC features, need filtering
	* header and rownames
* Kolodziejczyk: 
	* in space-separated file with header and rownames
	* header and rownames
	* Over 38k features, so will need to apply some filter; what did Svensson do?
* Patel: 
	* in tab-separated `txt` file
	* just header?
	* See notes in `readPatelData.R`, will move those here at some point
	* We only have access to the log TPM values at the moment, so we can go no further
* Zeisel:
	* in tab-separated `txt` file
	* issues:
		* Risso used both RPKM and full-quantile normalised data; I do think the data we have here are the raw counts
	* seems easiest just to use Risso's R code to obtain the counts
* Zheng:
	* Svensson already handled this, just using his data
