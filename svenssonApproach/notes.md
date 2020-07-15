# Plan

The plan is to run Svenssons's analysis on the (real) data sets of Risso et al. 2018. Not sure about the simulated data yet.

I will start by running the code from `Gene-wise_dispersion.ipynb`. Not sure this actually needs to be in a notebook. 


Svensson looked at:

* negative controls
* homogeneous cells
* heterogeneous cells

Which data sets are suitable for this investigation? 

* Several contain ERCC spike-ins:
	* Allen
	* Fletcher 

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
* Patel: 
	* in tab-separated `txt` file
	* just header?
* Zeisel:
	* in tab-separated `txt` file
	* issues:
		* Risso used both RPKM and full-quantile normalised data; I do think the data we have here are the raw counts
	* seems easiest just to use Risso's R code to obtain the counts
* Zheng:
	* Svensson already handled this, just using his data
