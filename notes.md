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
	sims/figures/fig6e-g/fitZinbLun.R (running)
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

### Patel

Can't actually repeat this analysis as it requires alignment of the raw reads with TopHat etc. Don't have the resources for this at the moment.

### Allen

Missing metadata not accessible from the library `scRNAseq`, which is providing the rest. The `allen_covariates_1000.Rmd` references `collection_date` which is not part of the metadata of the `SingleCellExperiment` object provided by `scRNAseq`.

The chunk `zinb_check_batch` in `allen_covariates_1000.Rmd` references relative paths to directories not created during execution of the preceding chunks.

### Zeisel

### Espresso

## Simulation scripts

### `simFunction.R`

`simFunction.log` shows that `zinbSimWrapper` only ran three times, rather than four as expected. Perhaps when it did run, it did make use of a cached result? All the files expected to be there, are there. 

Used:

	for i in $(grep "^# fig" simFunction.R | sed 's/# //'); do if [ ! -f $i ]; then echo "$i"; fi; done

to check that all files listed in `simFunction.R` comment were indeed created.

### `figuresPaper.Rmd`

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

`timeZinb.R` takes hours and hours to run, and was not completed when I terminated it, so I must investigate that 

#### `fitZifa_allen_10000.R`

### fig6ad-S13-S14

`time`'d `fitZinbLun.R`:

real    2790m39.814s
user    4076m18.857s
sys     1292m50.092s

And it still didn't complete!

### fig6e-g

### figS12


