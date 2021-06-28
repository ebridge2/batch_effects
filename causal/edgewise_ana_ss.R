require(ggplot2)
require(energy)
require(parallel)
require(tidyverse)
require(mltools)
ncores <- parallel::detectCores() - 1

cohort <- "CoRR"
parcellation <- "AAL"
modality <- "fMRI"
in.file <- sprintf('/base/data/dcorr/inputs_%s_%s_%s.rds', modality, parcellation, cohort)
preproc <- readRDS(in.file)

methods <- c("Raw", "conditional ComBat", "ComBat", "causal ComBat")
R <- 10000

cov.tab <- preproc.caus$covariates
datasets <- unique(cov.tab$Dataset)
# get subsets of data corresponding to the subset causal ComBat runs on
Xs <- lapply(methods, function(meth) {
  cov.tab.meth <- preproc[[meth]]$covariates
  idx.match <- (cov.tab.meth %>%
                  right_join(cov.tab %>% select(Subid, Session), by=c("Subid"="Subid", "Session"="Session")))$id
  preproc[[meth]]$graphs[idx.match,]
})
names(Xs) <- methods
d <- dim(Xs$Raw)[2]

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

Y = preproc.caus$covariates$Sex; Z.caus = preproc.caus$covariates$Age
DY <- dist(Y); DZ.caus <- dist(Z.caus)
res <- do.call(rbind, lapply(methods, function(meth) {
  if (meth == "Raw") {
    # if Raw, account for dataset in partial DCorr
    Z <- cbind(one_hot(data.table(Dataset=factor(preproc.caus$covariates$Dataset))), Age=Z.caus)
    Z.sc <- scale(Z)
    DZ <- dist(Z.sc)
  } else {
    DZ <- DZ.caus
  }
  X <- Xs[[meth]]
  # for each edge, pdcorr(edge, sex | age) or pdcorr(edge, sex | age, dataset)
  do.call(rbind, mclapply(1:d, function(i) {
    i.coord <- pos2coord(i, dim.mat=c(nv, nv))
    # check if coordinate in upper triangle
    if (i.coord[1] > i.coord[2]) {
      print(i)
      DX <- dist(X[,i])
      test <- pdcor.test(DX, DY, DZ, R=R)
      return(data.frame(Row=i.coord[1], Column=i.coord[2], Edge=i, Statistic=test$statistic,
                        # p-value computation from energy package doesn't make sense
                        p.value=(1 + sum(test$replicates >= test$statistic))/(1 + R)))
    } else {
      return(NULL)
    }
  }, mc.cores = ncores)) %>% mutate(Method=meth)
})) %>% mutate(Set="Causal Subset", Cohort=cohort)

saveRDS(res, sprintf('../data/dcorr/outputs_edgewise_ss_%s.rds', cohort))