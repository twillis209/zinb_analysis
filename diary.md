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


