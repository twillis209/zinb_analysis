# TODO (delete me)

* Write up details of technology used to generate all data
* Download a couple of droplet-based data sets and see if we can generate good fits without any filtering
* Find out how to fit ZINB; will try to stay in Python for now with `statsmodels`, but having difficulty atm
* Find a labelled data set

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
	* Kolodziejczyk has cells separated by batch and plate; separating only by plate led to a very poor fit, batch was no better
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
	* Over 38k features, so will need to apply some filter; what did Svensson do? Apparently no filter was applied
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

# Code

Can't get `statsmodels` ZINB model to fit, so trying `R`. First want to see if we can get a similar NB fit. 

# Analysis

## Spike-ins

Svensson got very good fits to the spike-in data using gene-wise dispersion coefficients, but we:

* got good fits to Zheng (droplet, used by Svensson) and Zeisel; Zeisel was done with STRT-Seq with UMIs
* poor fit to Allen, Fletcher, Kolodziejczyk

## Allen

## Kolodziejczyk

Fits were very poor when performed without filtering. Need to apply some kind of filter, should look at what was used in the original paper. Having now read it, it seems the data supplied on the Espresso website were already filtered.

ZI fits are good, with a bit of excessive zero prediction for means in the middle of the range. Notable that performance does not seem to differ between batches and media (i.e. aggregations of batches). What if we just aggregated all media, too? Doing that now.

So the fits of the three models (global dispersion, genewise dispersion, and ZI with genewise dispersion) seem not to vary regardless of whether the data set is:

* all batches aggregated
* batches aggregated by medium
* individual batches

Does this mean that these covariates have little influence on zero fraction? 

## Patel 

Contains no spike-ins, but does contain genes with symbols beginning with 'ERCC' :).

Let's try aggregating by patient and batch first.
