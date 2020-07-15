# Readme

Needed somewhere for notes organised by topic rather than by date (as in `diary.md`).

# Installation

Trying to install the version of zinbwave contemporary with the paper (0.1.1):

	> install_github("drisso/zinbwave@v0.1.1")
	Downloading GitHub repo drisso/zinbwave@v0.1.1
	Error in utils::download.file(url, path, method = method, quiet = quiet,  :
	  cannot open URL 'https://api.github.com/repos/drisso/zinbwave/tarball/v0.1.1'

There are neither tags nor releases at `drisso/zinbwave`, so I assume they have since removed access to v0.1.1.

Let's proceed with the current version. Starting by installing missing packages required for `allen_covariates_1000.Rmd`.
Ok, installed `zifa` under the `diss` environment. May be necessary to rename `diss`.

So far, calls to the zifa wrapper function take the longest, it seems.

# Parallel code

Not entirely sure we are making full use of ZINB-WaVE's parallel implementation. Will use the vignette example and the Allen data set to time execution of the `zinbwave` function.

Seems to be no difference between `MulticoreParam(1)` and `MulticoreParam(2)` in terms of execution time, at least for the data and parameters chosen. Of course, we have no idea what the parallel fraction is and whether it would differ on a different data set.

...so we are no further forward in knowing whether `zinbwave` is runnning in parallel!

The zinbwave vignette states that one needs to either:

1. invoke `register()`; or
2. pass a BPPARAM object to the zinbwave functions

Got the following output after running `time Rscript parallelTest.R`:

	real    3m6.920s
	user    4m3.574s
	sys     1m14.575s

As user+sys > real, we clearly had some parallel execution.

This was with `register(MulticoreParam(2))` at the top of the file. Let's see if it is necessary to specify this or if we can just rely on the default arguments to the zinbwave functions:

	real    1m22.070s
	user    1m39.792s
	sys     0m30.182s

Again user+sys > real, so we have parallel execution.

## Running scripts on ARC

Testing my use of the `BiocParallel::BatchtoolsParam` class by fitting ZINB to the Allen data set on the cluster, but it is taking far too long, which suggests it is probably running serially.

I was able to get a simple `bplapply` function running in parallel or at least it seemed so, as the wall time was less than the some of the reported user+sys time. The scripts that matter are apparently not running in parallel, though.

Tried using `mclapply` and it runs in parallel only on the head node, not when submitted using SGE.

See the ARC notes in the `math5871` repo for details of my investigation of parallelising `zinbwave` on the cluster.

### 'Dumb' parallel execution


# Analyses

## Data sets

### V1/Tasic/Allen

* Tasic et al. 
* Available in `scRNAseq` package
* The 'Allen' data set is a subset of the Tasic data. From the package manual: '`ReprocessedAllenData()` provides 379 cells from the mouse visual cortex, which is a subset of the data from Tasic et al. (2016).'
* The counts used in the analyses by Risso et al. are the default `tophat_counts` (confirmed this with use of `digest::digest`)
* Plate-based: sorted by FACS into a 96-well plate and amplified using the SMARTer kit
### S1/CA1 data set

* Zeisel et al.
* Downloaded.
* Fluidigm C1 microfluidics cell capture platform

### mESC/Kolodziejczyk data set

* Kolodziejczyk et al.
* Downloaded.
* May now be in `scRNAseq`
* Fluidigm C1 microfluidics cell capture platform

### OE/Fletcher data set

* Fletcher et al.
* Downloaded
* Identifying developmental trajectories of OE cells
* Figure 3
* FACS then Fluidigm C1 microfluidics cell capture platform

### PBMC/Zheng data set

* Zheng et al.
* Downloaded.
* Clustering of PBMCs
* Figure 3
* Droplet-based 10x Genomics platform 

### GBM/Patel data set

* Patel et al.
* SMART-Seq

### Applications of data sets

* Impact of normalization methods
	* All four data sets
	* Figure 4
* Accounting for batch effects
	* mESC and glioblastoma data sets
	* Figure 5
* Goodness-of-fit of the ZINB-WaVE model
	* All four data sets 
* ZINB-WaVE estimators are asymptotically unbiased and robust
	* Simulated data
* ZINB-WaVE is more accurate than state-of-the-art methods
	* Simulated data:
		* ZINB model
		* Lun and Marioni model (`cputime.R`)

R files in this repo:

	real_data/allen_covariates_1000.Rmd (problem with absent files in batch effect analysis)
	real_data/goodness_of_fit_allen.Rmd (running)
	real_data/goodness_of_fit_patel.Rmd 
	real_data/patel_covariates.Rmd
	real_data/espresso_covariates.Rmd
	real_data/goodness_of_fit_espresso.Rmd
	real_data/goodness_of_fit_zeisel.Rmd
	real_data/zeisel_covariates.Rmd
	real_data/zeisel_plots.R
	real_data/espresso_plots.R
	real_data/patel_plots.R
	real_data/allen_plots.R
	real_data/silhouette.R
	sims/figures/simFunction.R
	sims/figures/figuresPaper.Rmd
	sims/figures/fig6e-g/lunSim.R
	sims/figures/fig6e-g/fitZinbLun.R
	sims/figures/fig5-S10-S11-S15-S9/timeZinb.R
	sims/figures/fig5-S10-S11-S15-S9/fitZinb_bias_mse_ncells.R 
	sims/figures/fig5-S10-S11-S15-S9/fitZinb_bias_mse_allParam.R
	sims/figures/fig6ad-S13-S14/fitZifa_allen_10000.R
	sims/figures/fig6ad-S13-S14/fitZifa_zeisel_10000.R
	sims/figures/fig6ad-S13-S14/fitZifa.R
	sims/figures/fig6ad-S13-S14/fitZinb_corSilh.R
	sims/figures/fig6ad-S13-S14/fitZinb10000.R
	sims/cputime/cputime.R


### Leads to biologically meaningful clusters

* Figure 2: Three DR techniques applied to the V1 data set (Allen) and correlation of components with QC measures
	* Figure 2a-b: PCA (TC)
	* Figure 2c-d: ZIFA (TC)
	* Figure 2e-f: ZINB-WaVE	
* sFigure 1: Fig 2 analysis on S1/CA1 data set (Zeisel)
* sFigure 2: Fig 2 analysis on mESC data set (Kolodziejczyk/Espresso)
* sFigure 3: Fig 2 analysis on glioblastoma data set (Patel)
* sFigure 4: ZW-derived 2D representation of V1 data set (Allen) using 500, 2000, 5000, 10000 most variable genes
* sFigure 5: Average silhouette widths for clusters produced with PCA, ZIFA, and ZW in V1, S1/CA1, glioblastoma, mESC data sets

R files:

	real_data/zeisel_plots.R
	real_data/espresso_plots.R
	real_data/patel_plots.R
	real_data/allen_plots.R

### Leads to novel biological insights

* Figure 3a-c: lineage inference on the OE data set (Fletcher) with Slingshot 
	* Figure 3a: PCA without endpoint supervision 
	* Figure 3b: PCA with endpoint supervision 
	* Figure 3b: ZINB-WaVE without endpoint supervision 
* Figure 3d-e: PBMC data set (Zheng)
	* Figure 3d: 2D t-SNE depiction of clusters obtained from ZW-derived 10D representation of PBMC cells
	* Figure 3e: Heatmap of expression of marker genes for the 18 clusters identified
* sFigure 6: t-SNE representation of results of sequential k-means clustering on ZW-derived 10D representation of PBMC cells
* sFigure 7: Sequential k-means clustering on ZW-derived 2D representation of PBMC cells
	* sFigure 7a: 2D ZW plot
	* sFigure 7b: 2D PCA plot

### Impact of normalisation

* Figure 4: average silhouette width in the real data sets by DR and normalisation technique and 
	* Figure 4a: V1
	* Figure 4b: S1/CA1
	* Figure 4c: glioblastoma
	* Figure 4d: mESC
* sFigure 8: PCA on V1 data set 
* sFigure 9: ZIFA on V1 data set 
* sFigure 10: PCA on S1/CA1 data set 
* sFigure 11: ZIFA on S1/CA1 data set 
* sFigure 12: PCA on mESC data set 
* sFigure 13: ZIFA on mESC data set 
* sFigure 14: PCA on glioblastoma data set 
* sFigure 15: ZIFA on glioblastoma data set 

### Accounting for batch effects

* Figure 5: Accounting for batch effects with ZW:
	* Figure 5a-d: mESC data set (Espresso)
		* Figure 5a: 2D ZW DR with sample-level intercept
		* Figure 5b: 2D ZW DR with batch as sample-level covariate
		* Figure 5c: silhouette widths for medium clusters with and without batch covariate
		* Figure 5d: silhouette widths for batch clusters with and without batch covariate
	* Figure 5e-f: glioblastoma data set (Patel)
		* Figure 5e: 2D ZW with sample-level intercept
		* Figure 5f: 2D ZW with total no. of expressed genes as sample-level intercept
* sFigure 16: DR on mESC data set (Espresso), 'ZINB-WaVE and ComBat', comparing ZW to remove batch effect with normalisation method, ComBat:
	* sFigure 16a: PCA on FQ-normalised counts
	* sFigure 16b: PCA on ComBat-normalised counts
	* sFigure 16c: silhouette widths by medium and choice of PCA with raw, TC, and FQ, and with/out ComBat, ZW with/out batch covariate
	* sFigure 16d: silhouette widths by batch and choice of PCA with raw, TC, and FQ, and with/out ComBat, ZW with/out batch covariate
* sFigure 17: DR on glioblastoma data set (Patel):
	* sFigure 17a: PCA with FQ and ComBat 
	* sFigure 17b: ZW with batch as sample-level covariate (not satisfactory separation)

### Goodness-of-fit of ZINB-WaVE model

### ZINB-WaVE estimators are asymptotically unbiased and robust

### ZINB-WaVE is more accurate than state-of-the-art methods

## Scripts for data analysis

`run_zifa.py` provides a command line interface to ZIFA.

`silhouette.R` does silhouette calculations for the data sets in this directory.

### Workflow for real data scripts


* allen_covariates_1000.Rmd (missing metadata)
* allen_plots.R
* espresso_covariates.Rmd
* espresso_plots.R
* goodness_of_fit_allen.Rmd
* goodness_of_fit_espresso.Rmd
* goodness_of_fit_patel.Rmd
* goodness_of_fit_zeisel.Rmd
* patel_covariates.Rmd
* patel_plots.R
* zeisel_covariates.Rmd
* zeisel_plots.R

### Patel

Can't actually repeat this analysis as it requires alignment of the raw reads with TopHat etc. Don't have the resources for this at the moment.

#### `goodness_of_fit_patel.Rmd`

#### `patel_covariates.Rmd`

#### `patel_plots.R`

### Allen

Missing metadata not accessible from the library `scRNAseq`, which is providing the rest. The `allen_covariates_1000.Rmd` references `collection_date` which is not part of the metadata of the `SingleCellExperiment` object provided by `scRNAseq`.

#### `goodness_of_fit_allen.Rmd`

#### `allen_covariates_1000.Rmd`

The chunk `zinb_check_batch` references relative paths to directories not created during execution of the preceding chunks.

#### `allen_plots.R`

### Zeisel

#### `goodness_of_fit_zeisel.Rmd`

#### `zeisel_covariates.Rmd`

#### `zeisel_plots.R`
	
### Espresso

#### `espresso_covariates.Rmd`

#### `espresso_plots.R`

#### `goodness_of_fit_espresso.Rmd`

## Simulation scripts

### Workflow

* simFunction.R: (apparently complete, but will only really know once we come to use the simulated data and look at the figure)
	1. simulate from Allen data (complete)
	1. simulate from Zeisel data (complete)
	2. zeiselBiasMSECpuTime	(complete)
	2. zeiselMeanDifferencesS26 (complete)
* fig5-S10-S11-S15-S9:
	* fitZinb_bias_mse_allParam.R (running after filtering for zero-count samples)
	* fitZinb_bias_mse_ncells.R (running after filtering for zero-count samples)
	* timeZinb.R (complete; check output)
* fig6ad-S13-S14 (need to install ZIFA):
	* fitZifa_allen_10000.R
	* fitZifa.R
	* fitZifa_zeisel_10000.R
	* fitZinb10000.R (running on cluster, need to identify which jobs failed due to samples with zero counts)
	* fitZinb_corSilh.R
* fig6e-g (not complete):
	* fitZinbLun.R (running)
	* lunSim.R (complete; check output)
* figS12 (no scripts):

A recurring problem is that we are filtering the simulated data to remove genes with zero counts across all samples and samples with zero counts across all genes. How does this affect the results? I suppose it doesn't matter: we can't run it any other way.

### `simFunction.R`

`simFunction.log` shows that `zinbSimWrapper` only ran three times, rather than four as expected. Perhaps when it did run, it did make use of a cached result? All the files expected to be there, are there. 

Used:

	for i in $(grep "^# fig" simFunction.R | sed 's/# //'); do if [ ! -f $i ]; then echo "$i"; fi; done

to check that all files listed in `simFunction.R` comment were indeed created.

`zinbSimWrapper` simulates B data sets from a real data set and for a particular combination of parameters. Note that all B data sets are written out to one `rda` file.

One of the b^2 values used to generate Allen data set-based simulated data is 50, when it is given as 10 in the paper (and the corresponding value of b^2 used for the Zeisel data set is 10, too). 

Fitting ZINB to the Allen data set: 

	real    82m6.631s
	user    124m44.392s
	sys     24m15.627s

Fitting ZINB to the Zeisel data set: 

	real    185m17.744s
	user    292m20.967s
	sys     35m31.413s

Simulating data for the file `benchmark_simZeisel_nc%s_ratio%s_offs%s.rda` (test of `zinbSimWrapper`): 

	real    10m24.311s
	user    15m1.607s
	sys     1m0.416s

Testing `simulateFromAllenData`: `simulateFromAllenData(ncells=100, ratioSSW_SSB=1, gammapiOffset=0)`

Trying to run again and getting `Execution halted`. Tried

	loadAndFilterAllenData()

	zinbSimWrapper(...)

and got the error:

	Error in fitZinbAndCache(core, K = 2, epsilon = ngenes, commonDispersion = F) :
	  object 'zinb' not found
	Calls: zinbSimWrapper -> fitZinbAndCache
	Execution halted

Fiddled with the calls to `assign` and to `load`, `zinb` appears in the correct environment now.

#### `zeiselMeanDifferencesS26`

This function's memory usage keeps blowing up. 72GB at the last count. 

Ran to completion with max. mem. usage of 86.252GB.

### `figuresPaper.Rmd`

Produces the following figures:

* Figures 5, S10, S11
* Figure S12
* Figure S9
* Figures 6a-d, S13, S14
* Figures 6e-g
	
Data files required:

	fig5-S10-S11-S15-S9/simZeisel_nc1000_ratio1_offs2_fittedAll.rda
	fig5-S10-S11-S15-S9/simZeisel_nc1000_ratio1_offs2.rda

	figS12/simZeisel_nc1000_ratio1_offs2_fitted.rda
	figS12/simZeisel_nc1000_ratio1_offs2.rda

	fig5-S10-S11-S15-S9/simZeisel_nc50_ratio1_offs2_fitted.rda
	fig5-S10-S11-S15-S9/simZeisel_nc100_ratio1_offs2_fitted.rda
	fig5-S10-S11-S15-S9/simZeisel_nc500_ratio1_offs2_fitted.rda
	fig5-S10-S11-S15-S9/simZeisel_nc1000_ratio1_offs2_fitted.rda
	fig5-S10-S11-S15-S9/simZeisel_nc5000_ratio1_offs2_fitted.rda
	fig5-S10-S11-S15-S9/simZeisel_nc10000_ratio1_offs2_fitted.rda
	fig5-S10-S11-S15-S9/simZeisel_nc50_ratio1_offs2.rda
	fig5-S10-S11-S15-S9/simZeisel_nc100_ratio1_offs2.rda
	fig5-S10-S11-S15-S9/simZeisel_nc500_ratio1_offs2.rda
	fig5-S10-S11-S15-S9/simZeisel_nc1000_ratio1_offs2.rda
	fig5-S10-S11-S15-S9/simZeisel_nc5000_ratio1_offs2.rda
	fig5-S10-S11-S15-S9/simZeisel_nc10000_ratio1_offs2.rda

	fig6ad-S13-S14/simZeisel_nc100_ratio1_offs-1.5.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio1_offs0.5.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio1_offs2.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio5_offs-1.5.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio5_offs0.5.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio5_offs2.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio10_offs-1.5.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio10_offs0.5.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio10_offs2.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio1_offs-1.5.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio1_offs0.5.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio1_offs2.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio5_offs-1.5.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio5_offs0.5.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio5_offs2.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio10_offs-1.5.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio10_offs0.5.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio10_offs2.rda
	fig6ad-S13-S14/simAllen_nc100_ratio1_offs0.rda
	fig6ad-S13-S14/simAllen_nc100_ratio1_offs2.rda
	fig6ad-S13-S14/simAllen_nc100_ratio1_offs5.rda
	fig6ad-S13-S14/simAllen_nc100_ratio5_offs0.rda
	fig6ad-S13-S14/simAllen_nc100_ratio5_offs2.rda
	fig6ad-S13-S14/simAllen_nc100_ratio5_offs5.rda
	fig6ad-S13-S14/simAllen_nc100_ratio50_offs0.rda
	fig6ad-S13-S14/simAllen_nc100_ratio50_offs2.rda
	fig6ad-S13-S14/simAllen_nc100_ratio50_offs5.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio1_offs0.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio1_offs2.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio1_offs5.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio5_offs0.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio5_offs2.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio5_offs5.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio50_offs0.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio50_offs2.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio50_offs5.rda

	fig6ad-S13-S14/simZeisel_nc100_ratio1_offs-1.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio1_offs0.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio1_offs2_fitted.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio5_offs-1.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio5_offs0.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio5_offs2_fitted.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio10_offs-1.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio10_offs0.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc100_ratio10_offs2_fitted.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio1_offs-1.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio1_offs0.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio1_offs2_fitted.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio5_offs-1.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio5_offs0.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio5_offs2_fitted.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio10_offs-1.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio10_offs0.5_fitted.rda
	fig6ad-S13-S14/simZeisel_nc1000_ratio10_offs2_fitted.rda
	fig6ad-S13-S14/simAllen_nc100_ratio1_offs0_fitted.rda
	fig6ad-S13-S14/simAllen_nc100_ratio1_offs2_fitted.rda
	fig6ad-S13-S14/simAllen_nc100_ratio1_offs5_fitted.rda
	fig6ad-S13-S14/simAllen_nc100_ratio5_offs0_fitted.rda
	fig6ad-S13-S14/simAllen_nc100_ratio5_offs2_fitted.rda
	fig6ad-S13-S14/simAllen_nc100_ratio5_offs5_fitted.rda
	fig6ad-S13-S14/simAllen_nc100_ratio50_offs0_fitted.rda
	fig6ad-S13-S14/simAllen_nc100_ratio50_offs2_fitted.rda
	fig6ad-S13-S14/simAllen_nc100_ratio50_offs5_fitted.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio1_offs0_fitted.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio1_offs2_fitted.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio1_offs5_fitted.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio5_offs0_fitted.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio5_offs2_fitted.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio5_offs5_fitted.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio50_offs0_fitted.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio50_offs2_fitted.rda
	fig6ad-S13-S14/simAllen_nc1000_ratio50_offs5_fitted.rda

	fig6ad-S13-S14/simZeisel_nc10000_ratio5_offs2.rda
	fig6ad-S13-S14/simAllen_nc10000_ratio5_offs5.rda
	fig6ad-S13-S14/simZeisel_nc10000_ratio5_offs2_fitted.rda
	fig6ad-S13-S14/simAllen_nc10000_ratio5_offs5_fitted.rda

`eval_data` takes all `zifa` data with all four normalisation techniques, i.e. 

    load(paste0(pp, '_zifa.rda'))
    load(paste0(pp, '_zifaTC.rda'))
    load(paste0(pp, '_zifaTMM.rda'))
    load(paste0(pp, '_zifaFQ.rda'))

`sil` calls `eval_data` as well as loading some files of its own.

	fig6e-g/simLun_100.rda
	fig6e-g/simLun_100_ziadd0.33.rda
	fig6e-g/simLun_100_ziadd0.67.rda
	fig6e-g/simLun_1000.rda
	fig6e-g/simLun_1000_ziadd0.33.rda
	fig6e-g/simLun_1000_ziadd0.67.rda
	fig6e-g/simLun_10000.rda
	fig6e-g/simLun_10000_ziadd0.33.rda
	fig6e-g/simLun_10000_ziadd0.67.rda

	fig6e-g/simLun_100_fitted.rda
	fig6e-g/simLun_100_ziadd0.33_fitted.rda
	fig6e-g/simLun_100_ziadd0.67_fitted.rda
	fig6e-g/simLun_1000_fitted.rda
	fig6e-g/simLun_1000_ziadd0.33_fitted.rda
	fig6e-g/simLun_1000_ziadd0.67_fitted.rda
	fig6e-g/simLun_10000_fitted.rda
	fig6e-g/simLun_10000_ziadd0.33_fitted.rda
	fig6e-g/simLun_10000_ziadd0.67_fitted.rda

	fig6e-g/simLun_100_zifa.rda
	fig6e-g/simLun_100_ziadd0.33_zifa.rda
	fig6e-g/simLun_100_ziadd0.67_zifa.rda
	fig6e-g/simLun_1000_zifa.rda
	fig6e-g/simLun_1000_ziadd0.33_zifa.rda
	fig6e-g/simLun_1000_ziadd0.67_zifa.rda
	fig6e-g/simLun_10000_zifa.rda
	fig6e-g/simLun_10000_ziadd0.33_zifa.rda
	fig6e-g/simLun_10000_ziadd0.67_zifa.rda

	fig6e-g/simLun_100_zifaTC.rda
	fig6e-g/simLun_100_ziadd0.33_zifaTC.rda
	fig6e-g/simLun_100_ziadd0.67_zifaTC.rda
	fig6e-g/simLun_1000_zifaTC.rda
	fig6e-g/simLun_1000_ziadd0.33_zifaTC.rda
	fig6e-g/simLun_1000_ziadd0.67_zifaTC.rda
	fig6e-g/simLun_10000_zifaTC.rda
	fig6e-g/simLun_10000_ziadd0.33_zifaTC.rda
	fig6e-g/simLun_10000_ziadd0.67_zifaTC.rda

	fig6e-g/simLun_100_zifaFQ.rda
	fig6e-g/simLun_100_ziadd0.33_zifaFQ.rda
	fig6e-g/simLun_100_ziadd0.67_zifaFQ.rda
	fig6e-g/simLun_1000_zifaFQ.rda
	fig6e-g/simLun_1000_ziadd0.33_zifaFQ.rda
	fig6e-g/simLun_1000_ziadd0.67_zifaFQ.rda
	fig6e-g/simLun_10000_zifaFQ.rda
	fig6e-g/simLun_10000_ziadd0.33_zifaFQ.rda
	fig6e-g/simLun_10000_ziadd0.67_zifaFQ.rda

	fig6e-g/simLun_100_zifaTMM.rda
	fig6e-g/simLun_100_ziadd0.33_zifaTMM.rda
	fig6e-g/simLun_100_ziadd0.67_zifaTMM.rda
	fig6e-g/simLun_1000_zifaTMM.rda
	fig6e-g/simLun_1000_ziadd0.33_zifaTMM.rda
	fig6e-g/simLun_1000_ziadd0.67_zifaTMM.rda
	fig6e-g/simLun_10000_zifaTMM.rda
	fig6e-g/simLun_10000_ziadd0.33_zifaTMM.rda
	fig6e-g/simLun_10000_ziadd0.67_zifaTMM.rda

Assuming that 

### fig5-S10-S11-S15-S9

Figures 5a-d relate to the mESC data. Fig 5a and b depict two-dimensional plots of mESC cells after ZINB-WaVE has been used for DR:

	Figure 5a: default ZINB-WaVE with sample-level intercept

	Figure 5b: ZINB-WaVE with batch as known sample-level covariate

	Figure 5c: average silhouette widths by biological condition for ZINB-WaVE with and without batch covariate

	Figure 5d: average silhouette widths by batch for ZINB-WaVE with and without batch covariate

Figures 5e and f relate to the glioblastoma dataset:

	Figure 5e: default ZINB-WaVE with sample-level intercept

	Figure 5f: ZINB-WaVE with total no. of expressed genes as sample-level covariate

Figure S9 depicts two factors from ZIFA applied to the V1 dataset with: 

	a) no normalisation
	b) TC
	c) FQ 
	d) TMM 

Figure S10 depicts the first two PCs for the S1/CA1 dataset with:

	a) no normalisation
	b) TC
	c) FQ 
	d) TMM 

Figure S11 depicts the first two factors (i.e. of FA) for the S1/CA1 dataset with:

	a) no normalisation
	b) TC
	c) FQ 
	d) TMM 

Figure S15 depicts two factors from ZIFA applied to the glioblastoma dataset with: 

	a) no normalisation
	b) TC
	c) FQ 
	d) TMM 

#### `timeZinb.R`

`timeZinb.R` takes hours and hours to run, and was not completed when I terminated it, so I must investigate that 

Killed on cluster for exceeding requested 32GB. Now trying 64GB. Now trying 96GB.

Killed on cluster for reaching 104.032GB. What is going on with this? 

Trying with just the first element of the vector.

#### `fitZinb_bias_mse_ncells.R`

Completed with error conditions relating to zero counts in samples.

Now running on ARC4, get

	Error: BiocParallel errors
	  element index: 1, 2, 3, 4, 5, 6, ...
	  first error: Sample 2 has only 0 counts!
	Execution halted

This is 'thrown' (actually just calls `stop`) in `zinbInitialize`. Presumably something has changed in the codebase since the code in the paper repo was written. Perhaps we should filter it manually? `fitZinb_bias_mse_allParam` filters the counts, manually.

Now failing due to profligate virtual memory usage on the cluster, and that when run serially! Ran a single file which nonetheless failed with `h_vmem=50G`.

#### `fitZinb_bias_mse_allParam.R`

	Error: BiocParallel errors
	  element index: 1, 2, 3, 4, 5, 6, ...
	  first error: Registry must be writeable. See ?loadRegistry.
	Execution halted

Trying it serially (on a single core), instead.

Encountered the following error:

	Loading required package: SingleCellExperiment
	Error in zinbInitialize(m, Y, nb.repeat = nb.repeat.initialize, BPPARAM = BPPARAM) : 
	  Sample 770 has only 0 counts!
	Calls: source ... myZinbFit -> zinbFit -> zinbFit -> .local -> zinbInitialize
	Execution halted

### fig6ad-S13-S14

Figure 6 depicts bias and MSE estimates for the ZINB-WaVE estimation procedure as a function of the number of unknown covariates, K, with and without sample-level intercepts, and with common and genewise dispersion parameters. 10 bootstrap datasets were generated from data simulated from the ZINB-WaVE model based on the S1/CA1 data set with 1000 cells and 1000 genes.

	Figure 6a: bias of log(mu) estimate

	Figure 6b: bias of pi estimate

	Figure 6c: MSE of log(mu) estimate

	Figure 6d: MSE of pi estimate

Figure S13 depicts two factors from ZIFA applied to the mESC data set with: 

	a) no normalisation
	b) TC
	c) FQ 
	d) TMM 

Figure S14 depicts the first two PCs for the glioblastoma data set with:

	a) no normalisation
	b) TC
	c) FQ 
	d) TMM 

#### `fitZifa_allen_10000.R`

Running ZIFA is computationally intensive, finished with many errors.

#### `fitZifa_zeisel_10000.R`

Outcome is likely to be the same as the allen script.

#### `fitZifa.R`

As above, but not a single input file, many, many files instead.

#### `fitZinb_corSilh.R`

#### `fitZinb10000.R`

Not clear which files we have to run this on to replicate the paper.

### fig6e-g

Figure 6e-g is actually 7e-g in the publication.

#### `lunSim.R`

Appears to have run to completion on the cluster, all files expected are present.

#### `fitZinbLun.R`

`time`'d `fitZinbLun.R`:

real    2790m39.814s
user    4076m18.857s
sys     1292m50.092s

And it still didn't complete! From `fitZinbLun.log`:

	Error in result[[njob]] <- value : 
	  attempt to select less than one element in OneIndex
	Calls: lapply ... unlist -> bplapply -> bplapply -> bploop -> bploop.lapply
	In addition: There were 12 warnings (use warnings() to see them)
	Execution halted

On the last execution, wrote out

	simLun_100_ziadd0_fitted.rda
	simLun_100_ziadd0.33_fitted.rda
	simLun_100_ziadd0.67_fitted.rda
	simLun_1000_ziadd0_fitted.rda
	simLun_1000_ziadd0.33_fitted.rda
	simLun_1000_ziadd0.67_fitted.rda
	simLun_10000_ziadd0_fitted.rda

But it is not clear that procedure was completed for the last file, so I consider only the first six files written out to completion.

What difference would it make if we did not use the 10000 results?

Running `fitSimData(10000, '_ziadd0')`:

	real    1217m59.326s
	user    1921m18.701s
	sys     472m29.077s

Over 20 hours for a single file!

Running `fitSimData(10000, '_ziadd0.33')`, ran out of memory. Terminated with: 

	Error in mcfork(detached) :
	  unable to fork, possible reason: Cannot allocate memory
	Calls: fitSimData ... bploop.lapply -> .send_to -> .send_to -> <Anonymous> -> mcfork
	Execution halted
	Error while shutting down parallel: unable to terminate some child processes

and ran for:

	real    1094m51.486s
	user    1657m22.648s
	sys     372m50.779s

12 GB of RAM is clearly not enough, although I was running the `ziadd0.67` run at the same time, so it is possible that running these serially would work.

Running `fitSimData(10000, '_ziadd0.67')`:

real    1295m21.258s
user    1882m26.426s
sys     555m38.486s

Seems to have run to completion.

Tried again with `_ziadd0.33`, but terminated it earlier as I now have an ARC account. Running this will be my first job.

On the cluster on one core, the script ran for 8 hours without completion. Really will need to parallelise it.

### figS12

Figure S12 depicts the first two PCs for the mESC data set with:

	a) no normalisation
	b) TC
	c) FQ 
	d) TMM 
