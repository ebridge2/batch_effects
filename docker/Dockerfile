FROM bioconductor/bioconductor_docker:devel

RUN Rscript -e "install.packages(c('tidyverse', 'parallelDist', 'MatchIt', 'energy', 'multcomp', 'survey', 'dplyr'))"
RUN Rscript -e "BiocManager::install('sva')"
RUN Rscript -e "install.packages('igraph')"
RUN Rscript -e "install.packages('stringr')"
RUN Rscript -e "install.packages('mgcv')"
RUN Rscript -e "install.packages('entropy')"
RUN Rscript -e "install.packages(c('mltools'))"
