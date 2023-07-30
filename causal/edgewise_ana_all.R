# docker run -ti --entrypoint /bin/bash -v /cis/project/ndmg/batch_effects/:/data -v /cis/home/ebridge2/Documents/research/graphstats_repos/batch_effects/:/base neurodata/batch_effects:0.0.2
# docker run -ti --entrypoint /bin/bash -v /mnt/nfs2/batch_effects/:/data -v /home/eric/Documents/research/graphstats-repos/batch_effects/:/base neurodata/batch_effects:0.0.2
# cd /base/
require(ggplot2)
require(energy)
require(parallel)
require(tidyverse)
require(mltools)
require(data.table)
require(cdcsis)
source('./causalComBat.R')
ncores <- parallel::detectCores() - 1

cohort <- "CoRR"
parcellation <- "AAL"
modality <- "fMRI"
in.file <- sprintf('/base/data/dcorr/inputs_%s_%s_%s.rds', modality, parcellation, cohort)
preproc <- readRDS(in.file)

methods <- c("ComBat", "conditional NeuroH")#, #("Raw", "conditional ComBat")#, )

preproc.cb <- preproc$`ComBat`
cov.tab <- preproc.cb$covariates

R <- 1000
datasets <- unique(cov.tab$Dataset)

pos2coord<-function(pos=NULL, coord=NULL, dim.mat=NULL){
  if(is.null(pos) & is.null(coord) | is.null(dim.mat)){
    stop("must supply either 'pos' or 'coord', and 'dim.mat'")
  }
  if(is.null(pos) & !is.null(coord) & !is.null(dim.mat)){
    pos <- ((coord[,2]-1)*dim.mat[1])+coord[,1] 
    return(pos)
  }
  if(!is.null(pos) & is.null(coord) & !is.null(dim.mat)){
    coord <- matrix(NA, nrow=length(pos), ncol=2)
    coord[,1] <- ((pos-1) %% dim.mat[1]) +1
    coord[,2] <- ((pos-1) %/% dim.mat[1]) +1
    return(coord)
  }
}
# AAL parcellation has 116 vertices
nv <- 116

res <- do.call(rbind, lapply(methods, function(meth) {
  print(meth)
  X <- preproc[[meth]]$graphs
  cov.tab <- preproc[[meth]]$covariates
  Y <- cov.tab$Sex
  Z <- as.matrix(cov.tab$Age, ncol=1)
  # if (meth == "Raw") {
  #    # if Raw, account for dataset in conditional DCorr
  #    Z <- data.frame(Z)
  #    Z.sc <- scale(Z)
  # } else {
  #  Z.sc <- Z
  # }
  d <- dim(X)[2]
  # for each edge, cdcorr(edge, sex | age) or pdcorr(edge, sex | age)
  test = do.call(rbind, mclapply(1:d, function(i) {
    i.coord <- pos2coord(i, dim.mat=c(nv, nv))
    # check if coordinate in upper triangle
    if (i.coord[1] > i.coord[2]) {
      Z <- Z[, apply(Z, 2, var) > 0, drop=FALSE]
      DX <- dist(X[,i])
      #test <- cond.dcorr(DX, Y, Z, R=R, distance=TRUE)
      test <- gcm(X[,i,drop=FALSE], as.matrix(ohe(Y)), Z, R=R, regr.method="gam")
      if(i %% 100 == 0) {
        print(i)
      }
      return(data.frame(Row=i.coord[1], Column=i.coord[2], Edge=i, Statistic=test$statistic,
                        p.value=test$p.value))
    } else {
      return(NULL)
    }
  }, mc.cores = ncores - 1)) %>% mutate(Method=meth)
})) %>% mutate(Set="All", Cohort=cohort)

saveRDS(res, sprintf('../data/dcorr/outputs_edgewise_all_%s_CN.rds', cohort))