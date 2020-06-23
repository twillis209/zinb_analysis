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

# Analyses

Datasets (from Methods section):

	* V1 data set:
		* Tasic et al. 
		* Available in `scRNAseq` package
		* The 'Allen' data set is a subset of the Tasic data. From the package manual: '`ReprocessedAllenData()` provides 379 cells from the mouse visual cortex, which is a subset of the data from Tasic et al. (2016).'
	* S1/CA1 data set:
		* Zeisel et al.
		* Downloaded.
	* mESC data set:
		* Kolodziejczyk et al.
		* Downloaded.
		* May now be in `scRNAseq`
		* Identifying developmental trajectories of OE cells
		* Figure 3
		* Clustering of PBMCs
		* Figure 3
	* Impact of normalization methods
		* All four data sets
		* Figure 4
	* Accounting for batch effects
		* mESC and glioblastome data sets
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
	sims/figures/simFunction.R (complete)
	sims/figures/figuresPaper.Rmd
	sims/figures/fig6e-g/lunSim.R (complete)
	sims/figures/fig6e-g/fitZinbLun.R (terminated in error condition)
	sims/figures/fig5-S10-S11-S15-S9/timeZinb.R (not running; taking too long)
	sims/figures/fig5-S10-S11-S15-S9/fitZinb_bias_mse_ncells.R (complete)
	sims/figures/fig5-S10-S11-S15-S9/fitZinb_bias_mse_allParam.R (complete)
	sims/figures/fig6ad-S13-S14/fitZifa_allen_10000.R
	sims/figures/fig6ad-S13-S14/fitZifa_zeisel_10000.R
	sims/figures/fig6ad-S13-S14/fitZifa.R
	sims/figures/fig6ad-S13-S14/fitZinb_corSilh.R
	sims/figures/fig6ad-S13-S14/fitZinb10000.R
	sims/cputime/cputime.R

## Scripts for data analysis

`run_zifa.py` provides a command line interface to ZIFA.

`silhouette.R` does silhouette calculations for the data sets in this directory.

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

### `simFunction.R`

`simFunction.log` shows that `zinbSimWrapper` only ran three times, rather than four as expected. Perhaps when it did run, it did make use of a cached result? All the files expected to be there, are there. 

Used:

	for i in $(grep "^# fig" simFunction.R | sed 's/# //'); do if [ ! -f $i ]; then echo "$i"; fi; done

to check that all files listed in `simFunction.R` comment were indeed created.

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

#### `fitZinb_bias_mse_ncells.R`

Complete but must check output.

#### `fitZinb_bias_mse_allParam.R`

Complete but must check output.

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

### fig6e-g

Figure 6e-g is actually 7e-g in the publication.

#### `lunSim.R`

Appears to have run to completion.

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

### figS12

Figure S12 depicts the first two PCs for the mESC data set with:

	a) no normalisation
	b) TC
	c) FQ 
	d) TMM 
