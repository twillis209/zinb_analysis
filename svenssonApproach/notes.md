# TODO (delete me)

* Write up details of technology used to generate all data
* Download a couple of droplet-based data sets and see if we can generate good fits without any filtering

# Plan

The plan is to run Svenssons's analysis on the (real) data sets of Risso et al. 2018. Not sure about the simulated data yet. 

Update 18/8/20: Also adding the MATQ-Seq dataset

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
	* No spike-ins
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

Fitting on the cluster without filters at the moment, but we should implement some. Risso used 

	filter <- rowSums(all.counts>50)>=50

## MATQ-Seq data set

MATQ-Seq uses random hexamer UMIs which the authors refer to as 'amplicon indexes'. The authors removed dependence on sequencing depth by normalising amplicon counts by the total number of unique amplicons in the cell, used amplicons per million amplicons. The authors state that absolute transcript numbers can be counted, but that this requires sufficient sequencing depth. 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78968

MATQ-Seq allows detection of both poly-A mature mRNA and non-poly-A pre-mRNA. Not obvious whether we should use intron, exon, or gene counts for our purposes.

Aspects of study:

* 6 averaged one-fifth MCF10A single-cell samples with ERCC spike-ins
* Sequenced 10 HEK293T cells individually
* Pooled 40 HEK293T cells, split their mixed lysate into 40 parts, and sequenced 10 of these single-cell averages (idea is that we can determine technical variation associated with each gene)
* Sequenced another HEK293T clone: 10 single cells and 10 single cell averages
* Sequenced breast cancer cell line MCF10A
* 6 MCF10A cells individually
* 6 MCF10A single cell averages (i.e. pool-and-split as with HEK293T)
* Additional 38 HEK293T single cells and 10 pool-and-split averaged sc samples
* F-testing for excess biological variation over technical variation

Pipeline: TopHat, HT-seq, 

Starting with ERCC spike-in data. Barcodes not mentioned in the files so I don't actually know if UMIs were used with these; only 'reads' are listed.

Do we have a sufficient number of replicates to get decent fits? What units are they in?

Poorly annotated data, not sure if read or UMI counts, assume UMIs on the basis of reads. Too few samples, really only 6. How many did Svensson have for his analyses?

* Macosko ERCC samples: 84 droplets (80 ERCC spike-ins)
* Svensson (1): 2000 droplets (24116 sequences)
* Svensson (2): 2000 droplets (24116 sequences)
* Zheng: 1015 (92 ERCC spike-ins)
* Klein: 953 droplets (25435 sequences)
* Padovan-Merhar: 96 cells, (92 ERCC spike-ins)
* Sheng: 6 cells (92 spike-ins)

Not sure if library size normalisation is going to make much of a difference. Did Svensson correct for sequencing depth? Apparently not, total counts per cell/sample varies widely within data sets. Assumption is that the 'raw' count distribution was of interest, wonder to what extent this causes problems with the fit. Attempt to correct sequencing depth with DESeq did not appear to change anything.

Although I have taken data thus far from the data made available directory, might be easier to use the data packaged up more neatly in a sub-directory of the `code` directory. It appears identical to the data published with the paper.

High number of sequences in averaged data, but high compared to the Svensson data? Could check this if necessary.

Experiencing the odd warning when trying to fit the negative binomial. Warning is raised by `scipy`'s `optimize` function. Was just one or two genes, not worth worrying about.

Analyses we have:

* ERCC spike-ins (6 samples)
* Set one pool-and-split 293T cells (10 samples)
* Set two pool-and-split 293T cells (10 samples)
* 293T cells (38 samples)

Svensson scaled the Poisson data in a way we don't seem to be doing. The code uses a single valued for the scaled_count_mean (the first in the matrix?), which is like fitting a single Poisson model to all the data. Still, the current code can be used to reproduce the Macosko results, so it does seem correct. But I'm not happy with the scaling in the Sheng case. Perhaps the low replicate count is giving us the odd result? Or perhaps it was highly sequenced? The counts do seem high, even when scaled.

How to reproduce the Townes-style histogram analysis? Not really enough replicates to do it for the ERCC and pool-and-split samples. Not sure how to proceed in a rigorous fashion. 

Whether or not we should use the 'noFilter' data should be decided by however Svensson did it. Apparently Svensson did not use a filter.

10x v3 PBMC
Genes:  33538
15147
10x v3 HEK293T
Genes:  57905
30684
10x v3 NIH3T3
Genes:  54232
31952
Padovan-Merhar et al 2015 (SMARTer)
Genes:  25590
8573
Klein et al 2015
Genes:  25435
169
Macosko et al 2015
Genes:  80
0
Zheng et al 2017
Genes:  92
1

Sheng et al. 2017 (P&S)
Genes:  57819
39131
Sheng et al. 2017 (293T single cell)
Genes:  57819
26682
