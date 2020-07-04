
fitZinbAllen: 
	Rscript -e "setwd('sims/figures'); source('simFunction.R'); loadAndFilterAllenData(); fitZinbAndCache(core);"	

fitZinbZeisel:	
	Rscript -e "setwd('sims/figures'); source('simFunction.R'); loadAndFilterZeiselData(); fitZinbAndCache(counts);"	

allenSim: fitZinbAllen
	Rscript -e "setwd('sims/figures'); source('simFunction.R'); simulateFromAllenData();"	

zeiselSim: fitZinbZeisel
	Rscript -e "setwd('sims/figures'); source('simFunction.R'); simulateFromZeiselData();"	
