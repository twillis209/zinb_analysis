FROM bioconductor/bioconductor_docker:devel

# Cargo cult behaviour here, I'm afraid (saw it in the Bioconductor tutorial)
RUN apt-get update \
	&& apt-get install -y --no-install-recommends vim \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*
		
RUN git clone https://github.com/epierson9/ZIFA \
	&& cd ZIFA \
	&& python setup.py install

COPY /home/tom/diss/zinb_analysis 


RUN R -e 'BiocManager::install("zinbwave")'
