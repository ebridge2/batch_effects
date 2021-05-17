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
source('./causalComBat.R')

aal.homo.gr <- matrix(0, nrow=116, ncol=116)
sdiag(aal.homo.gr, k=1) <- rep(c(1, 0), 116/2)[1:(116-1)]
sdiag(aal.homo.gr, k=-1) <- rep(c(1, 0), 116/2)[1:(116-1)]
diag(aal.homo.gr) <- NaN
aal.homo.gr <- as.vector(aal.homo.gr)

aal.hetero.gr <- matrix(0, nrow=116, ncol=116)
for (i in 1:116) {
  for (j in 1:116) {
    aal.hetero.gr[i,j] <- ifelse((i %% 2) == (j %% 2), 1, 0)
  }
}
diag(aal.hetero.gr) <- NaN
aal.hetero.gr <- as.vector(aal.hetero.gr)

des.homo.gr <- matrix(0, nrow=70, ncol=70)
sdiag(des.homo.gr, k = 35) <- 1
sdiag(des.homo.gr, k=-35) <- 1
diag(des.homo.gr) <- NaN
des.homo.gr <- as.vector(des.homo.gr)

des.hetero.gr <- matrix(0, nrow=70, ncol=70)
des.hetero.gr[(70/2 + 1):70, (70/2 + 1):70] <- 1
des.hetero.gr[1:(70/2), 1:(70/2)] <- 1
diag(des.hetero.gr) <- NaN
des.hetero.gr <- as.vector(des.hetero.gr)

tri.gr <- matrix(0, nrow=116, ncol=116)
tri.gr[upper.tri(tri.gr)] <- 1
aal.comm <- list(Homotopic=aal.homo.gr, Homophilic=aal.hetero.gr, Upper.Tri=tri.gr)

tri.gr <- matrix(0, nrow=70, ncol=70)
tri.gr[upper.tri(tri.gr)] <- 1
des.comm <- list(Homotopic=des.homo.gr, Homophilic=des.hetero.gr, Upper.Tri=tri.gr)

parcel.comm <- list(AAL=aal.comm, Desikan=des.comm)


pairwise.driver <- function(graphs, cov.dat, parcellation="AAL", retain.dims=(as.vector(diag(116)) != 1),
                            mc.cores=1, R=1000) {
  datasets <- sort(unique((cov.dat %>% dplyr::select(Dataset))$Dataset))
  dset.pairs <- combn(datasets, 2)
  result <- mclapply(1:dim(dset.pairs)[2], function(x) {
    dset.1 <- dset.pairs[1,x]
    dset.2 <- dset.pairs[2,x]
    n.1 <- sum(cov.dat$Dataset == dset.1)
    n.2 <- sum(cov.dat$Dataset == dset.2)
    # Control is always larger of the two for all combn
    if (n.1 < n.2) {
      dset.i <- dset.1; dset.j <- dset.2
      n.i <- n.1; n.j <- n.2
    } else if (n.1 >= n.2) {
      dset.i <- dset.2; dset.j <- dset.1
      n.i <- n.2; n.j <- n.1
    }
    graphs.ij <- graphs[as.character(cov.dat$Dataset) %in% c(dset.i, dset.j),]
    cov.ij <- cov.dat %>%
      filter(Dataset %in% c(dset.i, dset.j)) %>%
      mutate(Treatment = as.numeric(Dataset == dset.i))
    tryCatch({
      if (length(unique(cov.ij$Sex)) > 1) {
        graphs.combt <- t(ComBat(t(graphs.ij), cov.ij$Dataset,
                                 mod = model.matrix(~as.factor(Sex) + Age, data=cov.ij)))
      } else {
        graphs.combt <- t(ComBat(t(graphs.ij), cov.ij$Dataset,
                                 mod = model.matrix(~ Age, data=cov.ij)))
      }
      Dmtx.norm <- as.matrix(dist(graphs.combt))
      
      result.site <- site_pair(Dmtx.norm, cov.ij, dset.i=dset.i, dset.j=dset.j, R=R)
      pdcorr.cov.sex <- pdcor.test(Dmtx.norm, y=cov.ij$Sex, z=cov.ij$Age, R=R)
      pdcorr.cov.age <- pdcor.test(Dmtx.norm, y=cov.ij$Age, z=cov.ij$Sex, R=R)
      result.cov <- data.frame(Dataset.Trt=dset.i, Dataset.Ctrl=dset.j, Effect.Name=c("Sex", "Age"),
                               Effect=c(pdcorr.cov.sex$estimate, pdcorr.cov.age$estimate),
                               p.value=c(pdcorr.cov.sex$p.value, pdcorr.cov.age$p.value))
      
      result.signal <- signal_ana(graphs.combt, cov.ij, parcellation=parcellation,
                                  mc.cores=1, retain.dims=retain.dims) %>%
        mutate(Dataset.Tgt=dset.i, Dataset.Ctrl=dset.j)
      
      return(list(Site=result.site, Covariate=result.cov, Signal=result.signal))
    }, error=function(e) {
      return(NULL)
    })
  }, mc.cores=mc.cores)
  
  res.site <- do.call(rbind, lapply(result, function(res) res$Site))
  res.cov <- do.call(rbind, lapply(result, function(res) res$Covariate))
  res.signal <- do.call(rbind, lapply(result, function(res) res$Signal))
  return(list(Site=res.site, Covariate=res.cov, Signal=res.signal))
}

singlenorm.driver <- function(gr.dat, gr.dat.full, cov.dat,
                              norm.options = c("Raw", "Ranked", "Z-Score", "ComBat", "causal ComBat"),
                              parcellation="AAL", retain.dims=(as.vector(diag(116)) != 1),
                              mc.cores=1, R=1000) {
  lapply(norm.options, function(norm) {
    print(norm)
    cov.post <- cov.dat
    if (norm == "ComBat") {
      dat.norm <- t(ComBat(t(gr.dat), cov.dat$Dataset))
    } else if (norm == "Z-Score") {
      dat.norm <- apply.along.dataset(gr.dat, cov.dat$Dataset, zsc.batch)
    } else if (norm == "Ranked") {
      dat.norm <- apply.along.dataset(gr.dat, cov.dat$Dataset, ptr.batch)
    } else if (norm == "Raw") {
      dat.norm <- gr.dat
    } else if (norm == "conditional ComBat") {
      mod <- model.matrix(as.formula("~as.factor(Sex) + Age"), data=cov.dat)
      dat.norm <- t(ComBat(t(gr.dat), cov.dat$Dataset, mod=mod))
    } else if (norm == "causal ComBat") {
      # asia.cohort <- which(cov.dat$Dataset %in% c("SWU4", "HNU1", "BNU3", "SWU1", "BNU2", "IPCAS1",
      #                                              "BNU1", "IPCAS6", "IPCAS3", "SWU2", "SWU3", "IPCAS4"))
      # am.cohort <- which(cov.dat$Dataset %in% c("NKI24tr645", "NKI24tr1400", "NKI24tr2500", "NYU1", "UWM",
      #                                           "Utah1", "MRN1", "IBATRT", "NYU2"))
      # norm.asia <- t(ComBat(t(gr.dat[asia.cohort,]), cov.dat$Dataset[asia.cohort]))
      # norm.am <- t(ComBat(t(gr.dat[am.cohort,]), cov.dat$Dataset[am.cohort]))
      # dat.norm <- rbind(norm.asia, norm.am)
      # cov.dat <- rbind(cov.dat[asia.cohort,], cov.dat[am.cohort,])
      
      am.cohort <- which(cov.dat$Dataset %in% c("NYU2", "IBATRT", "MRN1", "UWM", "NYU1"))
      
      caus.cb <- causal.ComBat(gr.dat[am.cohort,], cov.dat$Dataset[am.cohort], cov.dat[am.cohort,],
                               'as.factor(Sex) + Age', exact="Sex")
      dat.norm <- caus.cb$Data
      cov.post <- caus.cb$Covariates
      gr.dat.full <- gr.dat.full[am.cohort[caus.cb$Retained.Ids],]
    }
    cov.post <- cov.post %>% ungroup() %>% mutate(id=row_number())
    # exhaustively compute full distance matrix once since $$$
    Dmtx.norm <- as.matrix(parDist(dat.norm, threads=ncores))
    result.site <- causal_ana_site(Dmtx.norm, cov.post, mc.cores=ncores, R=R)
    result.cov <- causal_ana_cov(Dmtx.norm, cov.post, mc.cores=ncores, R=R)
    #result.cov.cont <- causal_ana_cov_cont(Dmtx.norm, cov.dat, mc.cores=ncores)
    result.signal <- signal_ana(dat.norm, cov.post, parcellation=parcellation, mc.cores=ncores,
                                retain.dims=retain.dims)
    gr.dat.norm <- gr.dat.full
    gr.dat.norm[,retain.dims] <- dat.norm
    gr.stats <- sum.stats(gr.dat.norm, cov.dat, n.vertices=n.vertices)
    
    return(list(Site=result.site %>% mutate(Method=norm),
                Covariate=result.cov %>% mutate(Method=norm),
                #Covariate.Cont=result.cov.cont %>% mutate(Method=norm),
                Signal=result.signal %>% mutate(Method=norm),
                Stats=gr.stats, D=Dmtx.norm, graphs.full=gr.dat.norm,
                Covariates=cov.post, dat.norm=dat.norm))
  })
}

compute_propensities <- function(df, form="Treatment ~ Sex + Age + Continent", trim=.01) {
  df$prop_scores <- glm(form, family=binomial(link="logit"), data=df)$fitted.values
  df <- df %>%
    mutate(weights=ifelse(Treatment == 1, 1, prop_scores/(1 - prop_scores))) %>%
    group_by(Treatment) %>%
    mutate(weights=weights/sum(weights)) %>%
    ungroup() %>%
    mutate(Stay=prop_scores >= trim)
  return(df)
}

adjusted.site_effect <- function(Dmtx.ij, cov.ij, form="as.factor(Treatment) ~ Age + as.factor(Sex)", dset.i="", dset.j="", ov.ij=0, R=1000) {
  tryCatch({
    cov.dat.ij.prop <- compute_propensities(cov.ij, form=form, trim=.01)
    # subset dmatrix by untrimmed data
    Dmtx.ij.trim <- Dmtx.ij[cov.dat.ij.prop$Stay, cov.dat.ij.prop$Stay]
    cov.dat.ij.trim <- cov.dat.ij.prop %>%
      filter(Stay == TRUE) %>%
      dplyr::select(Age, Sex, Continent, Treatment) %>%
      mutate(Age=as.numeric(Age), Sex=as.factor(Sex), Continent=as.factor(Continent))
    
    if (length(unique(cov.dat.ij.trim$Sex)) == 1) {
      form <- "as.factor(Treatment) ~ Age"
    }
    if ("Sex" %in% form) {
      exact="Sex"
    } else {
      exact=NULL
    }
    n.k <- max(floor(sum(cov.dat.ij.trim$Treatment == 0)/sum(cov.dat.ij.trim$Treatment == 1)), 1)
    match <- matchit(formula(form), data=cov.dat.ij.trim, method="nearest", exact=exact, ratio=n.k, caliper=.1)
    retain.ids <- as.numeric(names(match.data(match)$weights))
    Dmtx.cmp <- Dmtx.ij.trim[retain.ids, retain.ids]
    cov.dat.cmp <- cov.dat.ij.trim[retain.ids,]
    
    test.adj <- pdcor.test(as.dist(Dmtx.cmp), cov.dat.cmp$Treatment,
                           cov.dat.cmp %>% select(Sex, Age) %>% mutate(Sex=as.numeric(Sex), Age=as.numeric(Age)), R=R)
    return(data.frame(Data="Adjusted", Method="PDcor", Dataset.Trt=dset.i,
                      Dataset.Ctrl=dset.j, Effect.Name="Site", Effect=test.adj$estimate,
                      p.value=test.adj$p.value, Overlap=ov.ij))
  }, error=function(e) {
    print(e)
    return(NULL)
  })
}

site_pair <- function(Dmtx.dat.ij, cov.dat.ij, dset.i, dset.j, R=1000) {
  result <- list()
  ov.ij=compute_overlap(cov.dat.ij %>% filter(Dataset == dset.i), cov.dat.ij %>% filter(Dataset == dset.j))
  tryCatch({
    # Uncorrected
    test.uncor <- dcor.test(as.dist(Dmtx.dat.ij), cov.dat.ij$Treatment, R=R)
    result$uncor <- data.frame(Data="Associational", Method="Dcor", Dataset.Trt=dset.i,
                               Dataset.Ctrl=dset.j, Effect.Name="Site", Effect=test.uncor$estimate["dCor"],
                               p.value=test.uncor$p.value, Overlap=ov.ij)
    
    # If both datasets are NKI, Uncorrected == Causal Crossover
    if (grepl("NKI24", dset.i) & grepl("NKI24", dset.j)) {
      result$causal <- result$uncor
      result$causal$Data <- "Causal Cross."
    }
    
    ## Conditional
    test.cond <- pdcor.test(as.dist(Dmtx.dat.ij), cov.dat.ij$Treatment,
                             cov.dat.ij %>% dplyr::select(Continent, Sex, Age) %>%
                               mutate(Continent=as.numeric(Continent), Sex=as.numeric(Sex),
                                      Age=as.numeric(Age)), R=R)
    result$cond <- data.frame(Data="Conditional", Method="PDcor", Dataset.Trt=dset.i,
                               Dataset.Ctrl=dset.j, Effect.Name="Site", Effect=test.cond$estimate,
                               p.value=test.cond$p.value, Overlap=ov.ij)
    
    # Adjusted Approach. If there is only 1 sex amongst both datasets, omit
    if (length(unique((cov.dat.ij %>% dplyr::select(Sex))$Sex)) == 1) {
      form <- "as.factor(Treatment) ~ Age"
    } else {
      form <- "as.factor(Treatment) ~ Age + as.factor(Sex)"
    }
    if (unique((cov.dat.ij %>% filter(Dataset == dset.i))$Continent) == unique((cov.dat.ij %>% filter(Dataset == dset.j))$Continent)) {
      result$adj <- adjusted.site_effect(Dmtx.dat.ij, cov.dat.ij, form=form, dset.i=dset.i, dset.j=dset.j, R=R, ov.ij=ov.ij)
    }
    return(do.call(rbind, result) %>%
             mutate(Effect.Name="Site"))
  }, error=function(e) {
    print(e)
    return(data.frame(Data=c("Associational", "Causal Obs."), Method="PDcor", Dataset.Trt=dset.i,
              Dataset.Ctrl=dset.j, Effect.Name="Site", Effect=NA, p.value=NA,
              Overlap=compute_overlap(cov.dat.ij %>% filter(Dataset == dset.i),
                                      cov.dat.ij %>% filter(Dataset == dset.j))))
  })
}

causal_ana_site <- function(Dmtx.dat, cov.dat, trim=.01, R=1000, mc.cores=1) {
  datasets <- sort(unique((cov.dat %>% dplyr::select(Dataset))$Dataset))
  dset.pairs <- combn(datasets, 2)
  result <- do.call(rbind, mclapply(1:dim(dset.pairs)[2], function(x) {
    # get first and second dset
    dset.1 <- dset.pairs[1,x]
    dset.2 <- dset.pairs[2,x]
    n.1 <- sum(cov.dat$Dataset == dset.1)
    n.2 <- sum(cov.dat$Dataset == dset.2)
    # Control is always larger of the two for all combn
    if (n.1 < n.2) {
      dset.i <- dset.1; dset.j <- dset.2
      n.i <- n.1; n.j <- n.2
    } else if (n.1 >= n.2) {
      dset.i <- dset.2; dset.j <- dset.1
      n.i <- n.2; n.j <- n.1
    }
    cov.dat.ij <- cov.dat %>%
      filter(Dataset %in% c(dset.i, dset.j)) %>%
      mutate(Treatment = as.numeric(Dataset == dset.i))
    ov.ij=compute_overlap(cov.dat.ij %>% filter(Dataset == dset.i), cov.dat.ij %>% filter(Dataset == dset.j))
    Dmtx.dat.ij <- Dmtx.dat[cov.dat.ij$id, cov.dat.ij$id]
    return(site_pair(Dmtx.dat.ij, cov.dat.ij, dset.i, dset.j, R=R))
  }, mc.cores = mc.cores))
  return(result)
}

compute_effect <- function(D.i, cov.dat.i, E1.name, E2.name, R=1000, nboots=100, mc.cores=1, ...) {
  if (length(unique(cov.dat.i$Dataset)) == 1) {
    col.name <- "Dataset"
    col.id <- unique(as.character(cov.dat.i$Dataset))
  } else if (length(unique(cov.dat.i$Continent)) == 1) {
    col.name <- "Continent"
    col.id <- unique(as.character(cov.dat.i$Continent))
  } else {
    return(NULL)
  }
  normal.fn <- function(D.i, cov.dat.i, E1.name, E2.name, R=1000, nboots=1000, mc.cores=1) {
    n.i <- length(cov.dat.i$id)

    test.uncor <- pdcor.test(D.i, as.numeric(cov.dat.i[[E1.name]]),
                             as.numeric(cov.dat.i[[E2.name]]), R=R)

    jk.stat <- simplify2array(mclapply(1:n.i, function(i) {
      ids <- (1:n.i)[-i]
      cov.dat.boot <- cov.dat.i[ids,]
      Dmtx.dat.boot <- as.dist(as.matrix(D.i)[ids, ids])
      test.boot <- pdcor(Dmtx.dat.boot, as.numeric(cov.dat.boot[[E1.name]]),
                         as.numeric(cov.dat.boot[[E2.name]]))
      return(test.boot)
    }))
    psi <- test.uncor$estimate*n.i - (n.i - 1)*jk.stat
    ps.mean <- mean(psi); ps.var <- var(psi)
    ci.jk <- c(ps.mean + qnorm(.025)*sqrt(1/n.i*ps.var),
               ps.mean + qnorm(.975)*sqrt(1/n.i*ps.var))

    result <- data.frame(Method="PDcor", Effect.Name=E1.name, Effect=test.uncor$estimate, Effect.lwr.jk=ci.jk[1], Effect.upr.jk=ci.jk[2],
                         # Effect.lwr.npboot=ci.npboot[1], Effect.upr.npboot=ci.npboot[2],
                         p.value=test.uncor$p.value, Entropy=entropy(as.numeric(cov.dat.i[[E1.name]]), method="MM"),
                         Variance=var(as.numeric(cov.dat.i[[E1.name]])),
                         n=n.i, N=length(unique(cov.dat.i$Subid)))
    return(result)
  }
  result <- tryCatch(normal.fn(D.i, cov.dat.i, E1.name, E2.name, R=R, nboots=nboots, mc.cores=mc.cores), error=function(e) {
    result <- data.frame(Method="PDcor", Effect.Name=E1.name, Effect=NA, Effect.lwr.jk=NA, Effect.upr.jk=NA,
                      Effect.lwr.npboot=NA, Effect.upr.npboot=NA, Effect.lwr.adjboot=NA, Effect.upr.adjboot=NA,
                      p.value=NA, Entropy=NA, Variance=NA,
                      n=n.i, N=length(unique(cov.dat.i$Subid)))
    return(result)
  })
  result[[col.name]] <- col.id
  return(result)
}

causal_ana_cov <- function(Dmtx.dat, cov.dat, nboots=100, R=1000, mc.cores=1) {
  cov.dat <- cov.dat %>%
    ungroup() %>%
    mutate(id=row_number())
  datasets <- sort(unique((cov.dat %>% dplyr::select(Dataset))$Dataset))
  result.sex <- do.call(rbind, lapply(datasets, function(dataset) {
    tryCatch({
      cov.dat.dset <- cov.dat %>%
        filter(Dataset == dataset)

      n.i <- length(cov.dat.dset$id)
      Dmtx.dat.dset <- as.dist(Dmtx.dat[cov.dat.dset$id, cov.dat.dset$id])
      res <- compute_effect(Dmtx.dat.dset, cov.dat.dset, "Sex", "Age", R=R, nboots=nboots, mc.cores=mc.cores)

      return(res)
    }, error=function(e) {
      print(e)
      return(NULL)
    })

  }))

  result.age <- do.call(rbind, lapply(datasets, function(dataset) {
    tryCatch({
      cov.dat.dset <- cov.dat %>%
        filter(Dataset == dataset)

      n.i <- length(cov.dat.dset$id)
      Dmtx.dat.dset <- as.dist(Dmtx.dat[cov.dat.dset$id, cov.dat.dset$id])
      res <- compute_effect(Dmtx.dat.dset, cov.dat.dset, "Age", "Sex", R=R, nboots=nboots, mc.cores=mc.cores)

      return(res)
    }, error=function(e) {
      return(NULL)
    })
  }))
  return(rbind(result.sex, result.age))
}

causal_ana_cov_cont <- function(Dmtx.dat, cov.dat, nboots=1000, R=1000, mc.cores=1) {
  cov.dat <- cov.dat %>%
    ungroup() %>%
    mutate(id=row_number())
  continents <- sort(unique((cov.dat %>% dplyr::select(Continent))$Continent))
  result.sex <- do.call(rbind, lapply(continents, function(continent) {
    tryCatch({
      cov.dat.dset <- cov.dat %>%
        filter(Continent == continent)

      n.i <- length(cov.dat.dset$id)
      Dmtx.dat.dset <- as.dist(Dmtx.dat[cov.dat.dset$id, cov.dat.dset$id])
      res <- compute_effect(Dmtx.dat.dset, cov.dat.dset, "Sex", "Age", R=R, nboots=nboots,
                            mc.cores=mc.cores)

      return(res)
    }, error=function(e) {
      return(NULL)
    })
  }))

  result.age <- do.call(rbind, lapply(continents, function(continent) {
    tryCatch({
      cov.dat.dset <- cov.dat %>%
        filter(Continent == continent)

      n.i <- length(cov.dat.dset$id)
      Dmtx.dat.dset <- as.dist(Dmtx.dat[cov.dat.dset$id, cov.dat.dset$id])
      res <- compute_effect(Dmtx.dat.dset, cov.dat.dset, "Age", "Sex", R=R, nboots=nboots,
                            mc.cores=mc.cores)

      return(res)
    }, error=function(e) {
      return(NULL)
    })
  }))
  return(rbind(result.sex, result.age))
}

overlap_dist <- function(X) {
  datasets = levels(X$Dataset)
  D=sapply(unique(levels(X$Dataset)), function(dataseti) {
    sapply(unique(levels(X$Dataset)), function(datasetj) {
      suppressMessages(
        compute_overlap(X %>% filter(Dataset == dataseti) %>% ungroup(),
                        X %>% filter(Dataset == datasetj) %>% ungroup()))
    })
  })
  colnames(D) <- rownames(D) <- levels(X$Dataset)
  data.frame(Dataset1=colnames(D)[col(D)], Dataset2=rownames(D)[row(D)],
             Overlap=c(D)) %>%
    mutate(Dataset1=factor(Dataset1, levels=levels(X$Dataset), ordered=TRUE),
           Dataset2=factor(Dataset2, levels=levels(X$Dataset), ordered=TRUE)) %>%
    arrange(Dataset1, Dataset2)
}

compute_overlap <- function(X1, X2) {
  # probability of drawing two individuals with the same sex
  X1.sex <- X1 %>%
    group_by(Sex) %>%
    summarize(Per=n(), .groups="keep") %>%
    ungroup() %>%
    mutate(Per=Per/sum(Per))
  X2.sex <- X2 %>%
    group_by(Sex) %>%
    summarize(Per=n(), .groups="keep") %>%
    ungroup() %>%
    mutate(Per=Per/sum(Per))
  range.age <- c(min(X1$Age, X2$Age), max(X1$Age, X2$Age))
  per.ov <- sum(sapply(1:2, function(sex) {
    tryCatch({
      ov.sex.X1 <- (X1.sex %>% filter(Sex == sex))$Per
      ov.sex.X2 <- (X2.sex %>% filter(Sex == sex))$Per
      X1.sex.age <- (X1 %>% filter(Sex == sex) %>% select(Age))$Age
      X2.sex.age <- (X2 %>% filter(Sex == sex) %>% select(Age))$Age
      # obtain pdf for age
      X1.dens <- density(as.numeric(X1.sex.age), from=min(range.age), to=max(range.age))
      X1.dens$y <- X1.dens$y/sum(X1.dens$y)
      X2.dens <- density(as.numeric(X2.sex.age), from=min(range.age), to=max(range.age))
      X2.dens$y <- X2.dens$y/sum(X2.dens$y)
      ov.sex.age <- sum(pmin(X1.dens$y*ov.sex.X1, X2.dens$y*ov.sex.X2))
      return(ov.sex.age)
    }, error=function(e) {return(0)})
  }))
  return(as.numeric(unique((X1$Continent)) == unique(X2$Continent))*per.ov)
}

signal_ana <- function(data, cov.dat, parcellation="AAL", mc.cores=1, retain.dims=NULL) {
  if (is.null(retain.dims)) {
    stop("Pass removal indices.")
  }
  parc.c <- parcel.comm[[parcellation]]
  upper.tri.edges <- parcel.comm[[parcellation]]$Upper.Tri[retain.dims]
  do.call(rbind, lapply(c("Homotopic", "Homophilic"), function(community) {
    cmp.gr <- parc.c[[community]][retain.dims]
    do.call(rbind, mclapply(1:dim(data)[1], function(i) {
      tryCatch({
        dmeas <- data[i,]
        signal <- median(dmeas[cmp.gr == 1 & (upper.tri.edges == 1)]) - 
          median(dmeas[cmp.gr == 0 & (upper.tri.edges == 1)])
        test.res <- wilcox.test(dmeas[cmp.gr == 1 & (upper.tri.edges == 1)],
                                dmeas[cmp.gr == 0 & (upper.tri.edges == 1)], alternative = "greater")
        return(data.frame(Dataset=cov.dat[i,"Dataset"], Subid=cov.dat[i,"Subid"], Session=cov.dat[i,"Session"],
                          Signal=signal, Effect.Size=test.res$statistic, p.value=test.res$p.value,
                          Community=community, Parcellation=parcellation))
      }, error=function(e) {print(community)})
    }, mc.cores = mc.cores))
  }))
}

apply.along.dataset <- function(data, datasets, fn) {
  for (dataset in unique(datasets)) {
    this.dset <- data[datasets==dataset,]
    data[datasets == dataset,] <- apply(this.dset, c(2), fn)
  }
  return(data)
}

zsc.batch <- function(x, datasets) {
  return((x - mean(x))/sd(x))
}

ptr.batch <- function(x, datasets) {
  nz <- which(x != 0)
  x_rank = rank(x[nz])/(length(nz))
  x[nz] <- x_rank
  return(x)
}


summary_full_driver <- function(graphs, cov.dat, n.vertices) {
  avg.con <- matrix(apply(graphs, c(2), mean), nrow=n.vertices, ncol=n.vertices)
  
  avg.male <- matrix(apply(graphs[cov.dat$Sex == 2,], c(2), mean), nrow=n.vertices, ncol=n.vertices)
  avg.female <- matrix(apply(graphs[cov.dat$Sex == 1,], c(2), mean), nrow=n.vertices, ncol=n.vertices)
  
  age.cuts <- quantile(cov.dat$Age, probs=c(.2, .8))
  avg.young <- matrix(apply(graphs[cov.dat$Age <= age.cuts[1],], c(2), mean), nrow=n.vertices, ncol=n.vertices)
  avg.old <- matrix(apply(graphs[cov.dat$Age >= age.cuts[2],], c(2), mean), nrow=n.vertices, ncol=n.vertices)
  
  male.old <- matrix(graphs[which(cov.dat$Sex == 2 & cov.dat$Age > age.cuts[2])[1],], nrow=n.vertices, ncol=n.vertices)
  male.young <- matrix(graphs[which(cov.dat$Sex == 2 & cov.dat$Age <= age.cuts[1])[1],], nrow=n.vertices, ncol=n.vertices)
  
  female.old <- matrix(graphs[which(cov.dat$Sex == 1 & cov.dat$Age > age.cuts[2])[1],], nrow=n.vertices, ncol=n.vertices)
  female.young <- matrix(graphs[which(cov.dat$Sex == 1 & cov.dat$Age <= age.cuts[1])[1],], nrow=n.vertices, ncol=n.vertices)
  return(list(All=avg.con, Male=avg.male, Female=avg.female, Young=avg.young, Old=avg.old, 
              Female.Young=female.young, Female.Old=female.old, Male.Old=male.old, Male.Young=male.young))
}

summarize_over <- function(graphs, cov.dat, dimname, n.vertices) {
  unique.dim <- unique(cov.dat[[dimname]])
  res <- lapply(unique.dim, function(x) {
    gr.set <- graphs[cov.dat[[dimname]] == x,]
    cov.set <- cov.dat[cov.dat[[dimname]] == x,]
    
    return(summary_full_driver(gr.set, cov.set, n.vertices))
  })
  names(res) <- unique.dim
  return(res)
}

sum.stats <- function(graphs, cov.dat, n.vertices) {
  all.stats <- summary_full_driver(graphs, cov.dat, n.vertices)
  
  dset.stats <- summarize_over(graphs, cov.dat, "Dataset", n.vertices)
  cont.stats <- summarize_over(graphs, cov.dat, "Continent", n.vertices)
  return(list(All=all.stats, Dataset=dset.stats, Continent=cont.stats))
}
