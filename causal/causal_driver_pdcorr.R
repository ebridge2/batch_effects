# docker run -ti --entrypoint /bin/bash -v /cis/project/ndmg/batch_effects/:/data -v /cis/home/ebridge2/Documents/research/graphstats_repos/batch_effects/:/base neurodata/batch_effects:0.0.1

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
modality <- "fMRI"
mri.path <- file.path(in.path, modality, parcellation)

datasets=c("UWM", "NYU1", "Utah1", "MRN1", "IBATRT", "UPSM1", "NYU2",
           "BMB1", "IPCAS8", "IPCAS4", "SWU3", "SWU2", "IACAS1", "JHNU",
           "IPCAS5", "IPCAS2", "IPCAS3", "BNU1", "IPCAS1", "BNU2",
           "SWU1", "IPCAS7", "BNU3", "XHCUMS", "HNU1", "SWU4", "NKI24tr645", "NKI24tr1400",
           "NKI24tr2500")

gr.names <- list.files(path=mri.path, pattern="*.csv", recursive=TRUE)
gr.names <- gr.names[sapply(gr.names, function(name) any(sapply(datasets, function(dataset) {str_detect(name, dataset)})))]
vertices <- 1:n.vertices
fmt <- 'ncol'

list2array <- function(x) {
  good.ids <- !sapply(x, function(xi) is.null(xi))
  x <- x[good.ids]
  return(list(good.ids=good.ids, result=t(simplify2array(x))))
}

read.gr <- function(name, mri.path='', format='ncol') {
  tryCatch({
    g <- igraph::read_graph(file.path(mri.path, name), format=format)
    V.incl <- as.numeric(V(g)$name)
    V.notincl <- (1:116)[!sapply(1:116, function(i) i %in% V.incl)]
    if (length(V.notincl) > 0) {
      g <- add_vertices(g, length(V.notincl), attr=list(name=as.character(V.notincl)))
    }
    g.perm <- permute.vertices(g, as.numeric(V(g)$name))
    g.adj <- get.adjacency(g.perm, type="both", attr="weight", sparse=FALSE)
    diag(g.adj) <- 0
    return(as.vector(g.adj))
  }, error=function(e) {
    return(NULL)
  })
}
gr.out <- list2array(mclapply(gr.names, function(name) {
  read.gr(name,  mri.path=mri.path, format=fmt)
}, mc.cores=ncores))
gr.dat.full <- gr.out$result; good.ids <- gr.out$good.ids

cov.full <- read.csv(pheno.path)
spl.names <- strsplit(basename(gr.names[good.ids]), '_|-')
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

continent <- c("IBATRT"="North America", "Utah1"="North America", "IPCAS2"="Asia", "SWU1"="Asia", "UWM"="North America", "XHCUMS"="Asia", "SWU4"="Asia",
               "BNU2"="Asia", "IPCAS3"="Asia", "SWU3"="Asia", "IPCAS4"="Asia", "NYU2"="North America", "IPCAS1"="Asia",
               "IPCAS7"="Asia", "UPSM1"="North America", "IACAS1"="Asia", "IPCAS5"="Asia", "NYU1"="North America", "NYU2"="North America", "BNU1"="Asia",
               "MRN1"="North America", "BNU3"="Asia", "HNU1"="Asia", "SWU2"="Asia", "IPCAS8"="Asia", "JHNU"="Asia", "IPCAS6"="Asia",
               "BMB1"="Europe", "NKI24tr645"="North America", "NKI24tr1400"="North America", "NKI24tr2500"="North America")
cov.dat$Continent <- as.character(continent[as.character(cov.dat$Dataset)])

cov.dat <- cov.dat %>%
  mutate(Dataset=sub("_", "", Dataset))

# strip entries with no phenotypic data
retain.ids <- complete.cases(cov.dat)
cov.dat <- cov.dat[retain.ids,] %>%
  ungroup() %>% mutate(id=row_number(), Continent=factor(Continent),
                       Sex=factor(Sex))
gr.dat.full <- gr.dat.full[retain.ids,]

retain.dims <- sapply(1:dim(gr.dat.full)[2], function(j) {
  all(sapply(unique(cov.dat$Dataset), function(dataset) {
    var(gr.dat.full[cov.dat$Dataset == dataset,j]) > 0
  }))
})
gr.dat <- gr.dat.full[,retain.dims]

R=10000
norm.options <- c("Raw", "Ranked", "Z-Score", "ComBat", "cond. ComBat")

results <- single.norm.driver(gr.dat, gr.dat.full, cov.dat, norm.options=norm.options,
                              parcellation="AAL", retain.dims=retain.dims, mc.cores=mc.cores,
                              R=R)

gr.stats <- lapply(results, function(res) res$Stats)
names(gr.stats) <- norm.options
result.site <- do.call(rbind, lapply(results, function(res) res$Site)) %>%
  mutate(Modality=modality)
result.cov <- do.call(rbind, lapply(results, function(res) res$Covariate)) %>%
  mutate(Modality=modality)
# result.cov.cont <- do.call(rbind, lapply(results, function(res) res$Covariate.Cont)) %>%
#   mutate(Modality=modality)
result.signal <- do.call(rbind, lapply(results, function(res) res$Signal)) %>%
  mutate(Modality=modality)

preproc.obj <- lapply(results, function(res) {return(list(D=res$D, graphs=res$graphs.full))})
names(preproc.obj) <- norm.options

saveRDS(list(Site=result.site, Covariate=result.cov, #Covariate.Cont=result.cov.cont,
             Signal=result.signal, Stats=gr.stats, Covar.Tbl=cov.dat),
        file=sprintf('/base/data/dcorr/pdcorr_outputs_%s_%s.rds', modality, parcellation))
saveRDS(preproc.obj, file=sprintf('/base/data/dcorr/inputs_%s_%s.rds', modality, parcellation))

results.pairwise <- pairwise.driver(gr.dat, cov.dat, parcellation=parcellation, retain.dims=retain.dims,
                                    mc.cores=ncores, R=R)
saveRDS(results.pairwise, file=sprintf('/base/data/dcorr/pdcorr_pairwise_%s_%s.rds', modality, parcellation))