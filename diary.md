# 15/6/20

Trying to install the version of zinbwave contemporary with the paper (0.1.1):

	> install_github("drisso/zinbwave@v0.1.1")
	Downloading GitHub repo drisso/zinbwave@v0.1.1
	Error in utils::download.file(url, path, method = method, quiet = quiet,  :
	  cannot open URL 'https://api.github.com/repos/drisso/zinbwave/tarball/v0.1.1'

There are neither tags nor releases at `drisso/zinbwave`, so I assume they have since removed access to v0.1.1.

Analyses do these correspond in the paper?

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

Let's proceed with the current version. Starting by installing missing packages required for `allen_covariates_1000.Rmd`.

Missing ZIFA which needs to be installed through `conda`. Sidetracked by learning how `conda` works. 

Ok, installed `zifa` under the `diss` environment. May be necessary to rename `diss`.

So far, calls to the zifa wrapper function take the longest, it seems.

The chunk `zinb_check_batch` in `allen_covariates_1000.Rmd` references relative paths to directories not created during execution of the preceding chunks.

Issues (gathered):
	* Old parallelism used arguments like `ncores`, now using `BiocParallel`. I've used `htop` to confirm that jobs are using multiple cores.
	* `allen_covariates_1000.Rmd` is referencing non-existent relative paths as noted above
	
# 16/6/20

Trying to generate the simulated dataset. Get the following error with simFunction:

	Error in gzfile(file, "wb") : cannot open the connection
	Calls: zinbSimWrapper -> save -> gzfile
	In addition: Warning message:
	In gzfile(file, "wb") :
	  cannot open compressed file 'fig6ad-S13-S14/simAllen_nc100_ratio1_offs0.rda', probable reason 'No such file or directory'
	Execution halted

It's clear that the Rscripts and Rmd files are not sufficient to reproduce the analysis: the data is not to be found in the libraries they use (i.e. it doesn't all come from `scRNAseq`).	

Patel is out as we can't repeat the alignment. Let's look at Allen again. Need to get the metadata from somewhere. The library function being used to access the data is deprecated, so perhaps there are new ways to access the metadata.

Spent some time looking at the `SingleCellExperiment` class and its interface. The Allen metadata we need, the 'collection date', is absent from both the old version of the Allen dataset (accessed with `data("allen")`) and the new version (accessed with `ReprocessedAllenData()`) in `scRNAseq`. Risso loaded it from some file I don't have.

Note that `simFunction.R` takes a long time to run and depends on the data sets.

# 18/6/20

Ran `simFunction.R`. 
 
`sims/figures/fig5-S10-S11-S15-S9/timeZinb.R` takes hours and hours to run, and was not completed when I terminated it, so I must investigate that 

# 19/6/20

Prepended `sim` scripts with comments detailing the files they are intended to write out.

Working on running `fig6ad-S13-14` atm

	zinb_analysis/sims/figures/fig6ad-S13-S14/fitZifa_allen_10000.R : not running atm, problems with zifa.py paths
	zinb_analysis/sims/figures/fig6ad-S13-S14/fitZifa.R
	zinb_analysis/sims/figures/fig6ad-S13-S14/fitZifa_zeisel_10000.R
	zinb_analysis/sims/figures/fig6ad-S13-S14/fitZinb10000.R
	zinb_analysis/sims/figures/fig6ad-S13-S14/fitZinb_corSilh.R

Might be worth looking into caching of results and putting them somewhere in this project (rather than in `tmp` where they seem never to be found again).

`simFunction.log` shows that `zinbSimWrapper` only ran three times, rather than four as expected. Perhaps when it did run, it did make use of a cached result? All the files expected to be there, are there. 

Used:

	for i in $(grep "^# fig" simFunction.R | sed 's/# //'); do if [ ! -f $i ]; then echo "$i"; fi; done

to check that all files listed in `simFunction.R` comment were indeed created.


