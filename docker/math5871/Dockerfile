FROM bioconductor/bioconductor_docker:devel

# Cargo cult behaviour here, I'm afraid (saw it in the Bioconductor tutorial)
RUN apt-get update \
	&& apt-get install -y --no-install-recommends vim \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& apt install screen htop
		
RUN git clone https://github.com/epierson9/ZIFA \
	&& cd ZIFA \
	&& python setup.py install

RUN R -e 'BiocManager::install(c("zinbwave", "mclust", "scRNAseq"))'
RUN R -e 'BiocManager::install("EDASeq")'
RUN R -e 'install.packages("RColorBrewer")'
RUN R -e 'install.packages("batchtools")'
