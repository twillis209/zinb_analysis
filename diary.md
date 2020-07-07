# 15/6/20

Trying to install the version of zinbwave contemporary with the paper (0.1.1).

Let's proceed with the current version. Starting by installing missing packages required for `allen_covariates_1000.Rmd`.

Missing ZIFA which needs to be installed through `conda`. Sidetracked by learning how `conda` works. 
	
# 16/6/20

Became clear that the simulation scripts require the data (duh!). Spent some time acquiring those.

Spent some time looking at the `SingleCellExperiment` class and its interface. The Allen metadata we need, the 'collection date', is absent from both the old version of the Allen dataset (accessed with `data("allen")`) and the new version (accessed with `ReprocessedAllenData()`) in `scRNAseq`. Risso loaded it from some file I don't have.

Note that `simFunction.R` takes a long time to run and depends on the data sets.

# 18/6/20

Worked on running simulation code.
 
# 19/6/20

Prepended `sim` scripts with comments detailing the files they are intended to write out.

Working on running `fig6ad-S13-14` atm.

Might be worth looking into caching of results and putting them somewhere in this project (rather than in `tmp` where they seem never to be found again).

# 20/6/20

...and `fitZinbLun.R` is still running (almost 4 hours of CPU time). 

It's clear at this point that this repo is far from being a 'one push-button' means of reproducing the analysis in the Risso paper. There is not a script for every figure or analysis as some scripts indicate that some parameters need to be changed and the script rerun in order to produce the entirety of what was presented in the paper. I'm not sure how to go about this. 

Might be worth extending the lists of output files we expect from each script. Could do so by examining the script but we could also list all the figures and analyses in the paper and supplement...

# 21/6/20

`fitZinbLun.R` terminated today after a couple of days of wall time and that with termination in an error condition!

# 22/6/20

What would it mean to reproduce this paper? Surely not everything? Should think about what we need to do as part of the reproduction.

We also need to think about what the current scope of the dissertation is: reproduce (some of) Risso 18 then look at what use it makes of ZI in light of Svensson's conclusions. Should not forget that writing the review section of the dissertation is an opportunity to clarify my thinking and what should be done. 

# 27/6/20

Once we have things working (or some of the things...), we could repeat the GOF comparisons with the droplet-based data set from the paper (PBMCS, Zheng) as the V1 data (Allen) set was processed using a plate-based protocol. We have any number of other droplet-based data sets for use, too.

# 28/6/20

Currently running `testSimFunctionChanges.R` to fit ZINB to the Zeisel data set.

Need to fix `fitZinb_bias_mse_allParam.R` and `fitZinb_bias_mse_ncells.R`.

# 29/6/20

Useful tidbit for `vim`: `s/.*\/\([^\/]*\.R\)/\1` replaces full path with Rscript filename

Simulation-related scripts to fix:

* `fitZinb_bias_mse_allParam.R`
* `fitZinb_bias_mse_ncells.R`
* `fitZinbLun.R` (last few files take a long time to fit):
	* `fig6e-g/simLun_10000_ziadd0_fitted.rda` (running)
	* `fig6e-g/simLun_10000_ziadd0.33_fitted.rda`
	* `fig6e-g/simLun_10000_ziadd0.67_fitted.rda`
* `timeZinb.R` (not running; taking too long)
* `fitZifa_allen_10000.R`
* `fitZifa_zeisel_10000.R`
* `fitZifa.R`
* `fitZinb_corSilh.R`
* `fitZinb10000.R`
* `cputime.R`

Got a response from Risso about the Patel counts:

	I just added a file with the data at https://github.com/drisso/zinb_analysis/blob/master/real_data/GSE57872_GBM_data_matrix.txt.gz

	I am almost sure that these are the data that we used for the paper, but in case you see some discrepancies please re-open this and I will double check.

Unfortunately, these are the counts provided by Patel, not those obtained by Risso et al. through the use of TopHat and htseq-count.

## AWS

Not sure if it will be worth it, but will experiment with AWS today and see if it is worthwhile... No it's not! Free-tier only offers 1 vCPU, 1 GiB RAM. Nowhere near enough to improve on my current workflow. 

# 1/7/20

ZIFA is not working in my 'fork' of the `bioconductor/bioconductor_docker:devel` image. Must remember to run the unit tests before using it on ARC.

Got another response from Patel. He and the team are looking into the real counts.

# 2/7/20

Contents of the Singularity container are read-only, so perhaps it would be better to clone my `zinb_analysis` into home as I will be modifying the scripts therein.

# 4/7/20

Having difficulties running the code in `simFunction.R` in parallel, fitting ZINB to the Allen data set is taking far too long despite running the job on 8 cores with 32GB. 

# 6/7/20

## Cluster work

Being held back by the fact that I cannot run code in parallel.

Simulation scripts not working:

* fig5-S10-S11-S15-S9:
	* fitZinb_bias_mse_allParam.R: killed for sample with 0 counts
	* fitZinb_bias_mse_ncells.R: killed for sample with 0 counts
	* timeZinb.R: killed on cluster for reaching 104.032GB
* fig6ad-S13-S14: ZIFA not passing unit tests
	* fitZifa_allen_10000.R
	* fitZifa.R
	* fitZifa_zeisel_10000.R
	* fitZinb10000.R
	* fitZinb_corSilh.R
* fig6e-g:
	* fitZinbLun.R: ran 8 hours on one core on cluster without completion

I really must get `BatchtoolsParam` working on the cluster. Will read the `BiocParallel` and `BatchtoolsParam` documentation in-depth for some solution.

Might be worth trying the MPI example featured near the end of the `BiocParallel` to see if that will work.

Ok, I have gotten jobs to run in parallel (albeit by abandoning Singularity), but memory usage is blowing up with `MulticoreParam`. I will try `BatchtoolsParam` instead.

# 7/7/20

Spent whole day working on getting jobs parallelised on the cluster to some avail.

# 8/7/20

Check on most recent jobs on cluster with `exam`.
