FROM bioconductor/bioconductor_docker:devel

RUN Rscript -e "install.packages(c('tidyverse', 'parallelDist', 'MatchIt', 'energy', 'multcomp', 'survey', 'dplyr'))"
RUN Rscript -e "BiocManager::install('sva')"
RUN Rscript -e "install.packages('igraph')"
RUN Rscript -e "install.packages('stringr')"
RUN Rscript -e "install.packages('mgcv')"
RUN Rscript -e "install.packages('entropy')"
RUN Rscript -e "install.packages(c('mltools'))"
RUN Rscript -e "install.packages(c('cdcsis'))"
RUN python3 -m venv /opt/neuroharm
RUN git clone https://github.com/rpomponio/neuroHarmonize.git /neuroharm

RUN . /opt/neuroharm/bin/activate && pip install neuroHarmonize numpy
RUN Rscript -e "install.packages('reticulate')"

RUN . /opt/neuroharm/bin/activate && pip install pandas
RUN . /opt/neuroharm/bin/activate && pip install statsmodels nibabel
RUN Rscript -e "install.packages(c('GeneralisedCovarianceMeasure', 'weightedGCM'))"
