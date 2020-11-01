require(tidyverse)
require(MatchIt)
require(dplyr)
require(multcomp)
require(parallel)
require(survey)
require(energy)
require(igraph)
require(stringr)
require(parallelDist)
require(sva)
require(mgcv)
require(entropy)
source('./causal_investigation_helpers.R')

select <- dplyr::select
in.path <- '/data/'
n.vertices <- 116
pheno.path <- file.path(in.path, 'phenotypic/CoRR_AggregatedPhenotypicData.csv')
ncores <- detectCores() - 1
parcellation <- "AAL"
modality <- "dMRI"
mri.path <- file.path(in.path, modality, parcellation)


cov.dat.dset <- cov.dat %>%
  filter(Dataset == dataset)

Dmtx.dat.dset <- as.dist(Dmtx.dat[cov.dat.dset$id, cov.dat.dset$id])
compute_effect(Dmtx.dat.dset, cov.dat.dset, "Sex", "Age")

datasets=c("UWM", "NYU_1", "Utah1", "MRN1", "IBATRT", "UPSM_1", "NYU_2",
           "BMB1", "IPCAS_8", "IPCAS_4", "SWU3", "SWU2", "IACAS_1", "JHNU",
           "IPCAS_5", "IPCAS_2", "IPCAS_3", "BNU1", "IPCAS_1", "BNU2",
           "SWU1", "IPCAS_7", "BNU3", "XHCUMS", "HNU1", "SWU4")

gr.names <- list.files(path=fmri.path, pattern="*.csv", recursive=TRUE)
gr.names <- gr.names[sapply(gr.names, function(name) any(sapply(datasets, function(dataset) {str_detect(name, dataset)})))]
vertices <- 1:n.vertices
fmt <- 'ncol'

gr.dat <- t(simplify2array(mclapply(gr.names, function(name) {
  tgr <- get.adjacency(igraph::read_graph(paste0(fmri.path,name), format=fmt, predef=vertices), type="both", attr="weight", sparse=FALSE)
  diag(tgr) <- 0
  return(as.vector(tgr))
}, mc.cores=ncores)))

cov.full <- read.csv(pheno.path)
spl.names <- strsplit(basename(gr.names), '_|-')
dset.names <- strsplit(gr.names, '/')
cov.dat <- do.call(rbind, mclapply(1:length(spl.names), function(i) {
  spl.name <- spl.names[[i]]; dset.name <- dset.names[[i]]
  subid <- as.integer(spl.name[2]); sesid <- as.integer(spl.name[4])
  dset <- dset.name[length(dset.name)-1]
  cov.sc=cov.full %>%
    mutate(ID=row_number()) %>%
    filter(SUBID == subid) %>%
    mutate(Invalid.Entries=as.numeric(as.character(SEX) == "#") +
             as.numeric(as.character(AGE_AT_SCAN_1) == "#")) %>%
    ungroup() %>%
    filter(Invalid.Entries == min(Invalid.Entries)) %>%
    filter(ID==min(ID)) %>%
    dplyr::select(SEX, AGE_AT_SCAN_1) %>%
    rename(Sex=SEX, Age=AGE_AT_SCAN_1) %>%
    mutate(Subid=subid, Session=sesid, Dataset=dset,
           Sex=as.numeric(as.character(Sex)), Age=as.numeric(as.character(Age)))
  if (dim(cov.sc)[1] == 0) {
    return(data.frame(Subid=subid, Session=sesid, Dataset=dset, Sex=NA, Age=NA))
  } else {
    return(cov.sc)
  }
}, mc.cores=ncores))

continent <- c("IBATRT"="North America", "Utah1"="North America", "IPCAS_2"="Asia", "SWU1"="Asia", "UWM"="North America", "XHCUMS"="Asia", "SWU4"="Asia",
               "BNU2"="Asia", "IPCAS_3"="Asia", "SWU3"="Asia", "IPCAS_4"="Asia", "NYU2"="North America", "IPCAS_1"="Asia",
               "IPCAS_7"="Asia", "UPSM_1"="North America", "IACAS_1"="Asia", "IPCAS_5"="Asia", "NYU_1"="North America", "NYU_2"="North America", "BNU1"="Asia",
               "MRN1"="North America", "BNU3"="Asia", "HNU1"="Asia", "SWU2"="Asia", "IPCAS_8"="Asia", "JHNU"="Asia", "IPCAS_6"="Asia",
               "BMB1"="Europe")
cov.dat$Continent <- as.character(continent[as.character(cov.dat$Dataset)])

# strip entries with no phenotypic data
retain.ids <- complete.cases(cov.dat)
cov.dat <- cov.dat[retain.ids,] %>%
  ungroup() %>% mutate(id=row_number(), Continent=factor(Continent),
                       Sex=factor(Sex))
gr.dat <- gr.dat[retain.ids,]

retain.dims <- sapply(1:dim(gr.dat)[2], function(j) var(gr.dat[,j])) > 0
gr.dat <- gr.dat[,retain.dims]

results.cont <- lapply(c("Raw", "Ranked", "Z-Score", "ComBat"), function(norm) {
  if (norm == "ComBat") {
    dat.norm <- t(ComBat(t(gr.dat), cov.dat$Dataset))
  } else if (norm == "Z-Score") {
    dat.norm <- apply.along.dataset(gr.dat, cov.dat$Dataset, zsc.batch)
  } else if (norm == "Ranked") {
    dat.norm <- apply.along.dataset(gr.dat, cov.dat$Dataset, ptr.batch)
  } else {
    dat.norm <- gr.dat
  }
  # exhaustively compute full distance matrix once since $$$
  Dmtx.norm <- as.matrix(parDist(dat.norm, threads=ncores))

  result.site <- causal_ana_site(Dmtx.norm, cov.dat, mc.cores=ncores)
  result.cov <- causal_ana_cov(Dmtx.norm, cov.dat, mc.cores=ncores)
  result.cov.cont <- causal_ana_cov_cont(Dmtx.norm, cov.dat, mc.cores=ncores)
  result.signal <- signal_ana(dat.norm, cov.dat, parcellation=parcellation, mc.cores=ncores,
                              retain.dims=retain.dims)
  return(list(Site=result.site %>% mutate(Method=norm),
              Covariate=result.cov %>% mutate(Method=norm),
              Covariate.Cont=result.cov.cont %>% mutate(Method=norm),
              Signal=result.signal %>% mutate(Method=norm)))
})

result.site <- do.call(rbind, lapply(results, function(res) res$Site)) %>%
  mutate(Modality=modality)
result.cov <- do.call(rbind, lapply(results, function(res) res$Covariate)) %>%
  mutate(Modality=modality)
# result.cov.cont <- do.call(rbind, lapply(results, function(res) res$Covariate.Cont)) %>%
#   mutate(Modality=modality)
result.signal <- do.call(rbind, lapply(results, function(res) res$Signal)) %>%
  mutate(Modality=modality)

saveRDS(list(Site=result.site, Covariate=result.cov, #Covariate.Cont=result.cov.cont,
             Signal=result.signal),
        file=sprintf('/base/data/dcorr/pdcorr_outputs_%s_%s.rds', modality, parcellation))
