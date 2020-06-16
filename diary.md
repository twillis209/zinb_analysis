# 15/6/20

Trying to install the version of zinbwave contemporary with the paper (0.1.1):

	> install_github("drisso/zinbwave@v0.1.1")
	Downloading GitHub repo drisso/zinbwave@v0.1.1
	Error in utils::download.file(url, path, method = method, quiet = quiet,  :
	  cannot open URL 'https://api.github.com/repos/drisso/zinbwave/tarball/v0.1.1'

There are neither tags nor releases at `drisso/zinbwave`, so I assume they have since removed access to v0.1.1.

We have the following to work with: to which analyses do these correspond in the paper?

Datasets (from Methods section):

	* V1 data set:
		* Tasic et al. 
	* S1/CA1 data set:
		* Zeisel et al.
	* mESC data set:
		* Kolodziejczyk et al.
	* Glioblastoma data set:
		* Patel et al.
	* OE data set:
		* Fletcher et al.
	* PBMC data set:
		* Zheng et al.

Results (by heading in the section of the same name): 
	* ZINB-WaVE leads to biologically meaningful clusters
		* Correlation of components/factors from PC, ZIFA, ZINB with technical features in V1 data set:
		* Figure 2
		* Analysis was repeated on the S1/CA1 data set
		* supp. figs 1-3
	* ZINB-WaVE leads to novel biological insights
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
			* Lun and Marioni model

	



Analyses in the paper:

	real_data/allen_covariates_1000.Rmd
	real_data/goodness_of_fit_allen.Rmd (running)
	real_data/goodness_of_fit_patel.Rmd
	real_data/patel_covariates.Rmd
	real_data/espresso_covariates.Rmd
	real_data/goodness_of_fit_espresso.Rmd
	real_data/goodness_of_fit_zeisel.Rmd
	real_data/zeisel_covariates.Rmd

Let's proceed with the current version. Starting by installing missing packages required for `allen_covariates_1000.Rmd`.

Missing ZIFA which needs to be installed through `conda`. Sidetracked by learning how `conda` works. 

Ok, installed `zifa` under the `diss` environment. May be necessary to rename `diss`.

So far, calls to the zifa wrapper function take the longest, it seems.

The chunk `zinb_check_batch` in `allen_covariates_1000.Rmd` references relative paths to directories not created during execution of the preceding chunks.

Issues (gathered):

	* Old parallelism used arguments like `ncores`, now using `BiocParallel`
	* `allen_covariates_1000.Rmd` is referencing non-existent relative paths as noted above
	



# 16/6/20

Trying to generate the simulated dataset. Get the following error with simFunction:

	Error in gzfile(file, "wb") : cannot open the connection
	Calls: zinbSimWrapper -> save -> gzfile
	In addition: Warning message:
	In gzfile(file, "wb") :
	  cannot open compressed file 'fig6ad-S13-S14/simAllen_nc100_ratio1_offs0.rda', probable reason 'No such file or directory'
	Execution halted


