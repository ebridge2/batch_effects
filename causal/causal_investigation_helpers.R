
aal.homo.gr <- matrix(0, nrow=116, ncol=116)
sdiag(aal.homo.gr, k=1) <- rep(c(1, 0), 116/2)[1:(116-1)]
sdiag(aal.homo.gr, k=-1) <- rep(c(1, 0), 116/2)[1:(116-1)]
aal.homo.gr <- as.vector(aal.homo.gr)
diag(aal.homo.gr) <- NaN

aal.hetero.gr <- matrix(0, nrow=116, ncol=116)
for (i in 1:116) {
  for (j in 1:116) {
    aal.hetero.gr[i,j] <- ifelse((i %% 2) == (j %% 2), 1, 0)
  }
}
diag(aal.hetero.gr) <- NaN

des.homo.gr <- matrix(0, nrow=70, ncol=70)
sdiag(des.homo.gr, k = 35) <- 1
sdiag(des.homo.gr, k=-35) <- 1
des.homo.gr <- as.vector(des.homo.gr)
diag(des.homo.gr) <- NaN

des.hetero.gr <- matrix(0, nrow=70, ncol=70)
des.hetero.gr[(70/2 + 1):70, (70/2 + 1):70] <- 1
des.hetero.gr[1:(70/2), 1:(70/2)] <- 1
des.hetero.gr <- as.vector(des.hetero.gr)
diag(des.hetero.gr) <- NaN

aal.comm <- list(Homotopic=aal.homo.gr, Homophilic=aal.hetero.gr)
des.comm <- list(Homotopic=des.homo.gr, Homophilic=des.hetero.gr)
parcel.comm <- list(AAL=aal.comm, Desikan=des.comm)

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

causal_ana_site <- function(Dmtx.dat, cov.dat, trim=.01, R=500, mc.cores=1) {
  datasets <- sort(unique((cov.dat %>% dplyr::select(Dataset))$Dataset))
  dset.pairs <- combn(datasets, 2)
  result <- do.call(rbind, mclapply(1:dim(dset.pairs)[2], function(x) {
    result <- list()
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
    tryCatch({
      Dmtx.dat.ij <- Dmtx.dat[cov.dat.ij$id, cov.dat.ij$id]
      ## Uncorrected
      test.uncor <- pdcor.test(as.dist(Dmtx.dat.ij), cov.dat.ij$Treatment,
                               cov.dat.ij %>% dplyr::select(Continent, Sex, Age) %>%
                                 mutate(Continent=as.numeric(Continent), Sex=as.numeric(Sex),
                                        Age=as.numeric(Age)), R=R)
      result$uncor <- data.frame(Data="Untrimmed", Method="PDcor", Dataset.Trt=dset.i,
                                 Dataset.Ctrl=dset.j, Effect.Name="Site", Effect=test.uncor$estimate,
                                 p.value=test.uncor$p.value, Overlap=ov.ij)

      ## Trimmed
      if (length(unique((cov.dat.ij %>% filter(Dataset == dset.i) %>% dplyr::select(Sex))$Sex)) == 1) {
        form <- "Treatment ~ Age"
      } else {
        form <- "Treatment ~ Age + Sex"
      }
      if (unique((cov.dat.ij %>% filter(Dataset == dset.i))$Continent) == unique((cov.dat.ij %>% filter(Dataset == dset.j))$Continent)) {
        cov.dat.ij.prop <- compute_propensities(cov.dat.ij, form=form, trim=.01)
        # subset dmatrix by untrimmed data
        Dmtx.dat.ij.trim <- Dmtx.dat.ij[cov.dat.ij.prop$Stay, cov.dat.ij.prop$Stay]
        cov.dat.ij.trim <- cov.dat.ij.prop %>%
          filter(Stay == TRUE) %>%
          dplyr::select(Age, Sex, Continent, Treatment) %>%
          mutate(Age=as.numeric(Age), Sex=as.numeric(Sex), Continent=as.numeric(Continent))
        if (dim(cov.dat.ij.trim)[1] > n.i) {
          test.trim <- pdcor.test(as.dist(Dmtx.dat.ij.trim), cov.dat.ij.trim$Treatment,
                                  as.matrix(cov.dat.ij.trim %>% dplyr::select(Continent, Sex, Age)), R=R)
          result$trim <- data.frame(Data="Trimmed", Method="PDcor", Dataset.Trt=dset.i,
                                    Dataset.Ctrl=dset.j, Effect.Name="Site", Effect=test.trim$estimate,
                                    p.value=test.trim$p.value, Overlap=ov.ij)
        }
      }
      return(do.call(rbind, result) %>%
               mutate(Effect.Name="Site"))
    }, error=function(e) {
      return(data.frame(Data=c("Untrimmed", "Trimmed"), Method="PDcor", Dataset.Trt=dset.i,
                        Dataset.Ctrl=dset.j, Effect.Name="Site", Effect=NA, p.value=NA, Overlap=compute_overlap(cov.dat.ij %>% filter(Dataset == dset.i),
                                                                                                               cov.dat.ij %>% filter(Dataset == dset.j))))
    })
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
  tryCatch({
    n.i <- length(cov.dat.i$id)

    test.uncor <- pdcor.test(D.i, as.numeric(cov.dat.i[[E1.name]]),
                             as.numeric(cov.dat.i[[E2.name]]), R=R)
    ## Bootstrapped CIs
    bootstr.ci <- simplify2array(mclapply(1:nboots, function(i) {
      ids <- sample(1:n.i, size=n.i, replace = TRUE)
      cov.dat.boot <- cov.dat.i[ids,]
      Dmtx.dat.boot <- as.dist(as.matrix(D.i)[ids, ids])
      test.boot <- pdcor(Dmtx.dat.boot, as.numeric(cov.dat.boot[[E1.name]]),
                         as.numeric(cov.dat.boot[[E2.name]]))
      return(test.boot)
    }, mc.cores=mc.cores))
    ci.npboot <- quantile(bootstr.ci, c(.025, .975))

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
                      Effect.lwr.npboot=ci.npboot[1], Effect.upr.npboot=ci.npboot[2],
                      p.value=test.uncor$p.value, Entropy=entropy(as.numeric(cov.dat.i[[E1.name]]), method="MM"), Variance=var(as.numeric(cov.dat.i[[E1.name]])),
                      n=n.i, N=length(unique(cov.dat.i$Subid)))
  }, error=function(e) {
    result <- data.frame(Method="PDcor", Effect.Name=E1.name, Effect=NA, Effect.lwr.jk=NA, Effect.upr.jk=NA,
                      Effect.lwr.npboot=NA, Effect.upr.npboot=NA, Effect.lwr.adjboot=NA, Effect.upr.adjboot=NA,
                      p.value=NA, Entropy=NA, Variance=NA,
                      n=n.i, N=length(unique(cov.dat.i$Subid)))
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

causal_ana_cov_cont <- function(Dmtx.dat, cov.dat, nboots=100, R=1000, mc.cores=1) {
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
  if (is.null(retain.ids)) {
    stop("Pass removal indices.")
  }
  parc.c <- parcel.comm[[parcellation]]
  do.call(rbind, lapply(names(parc.c), function(community) {
    cmp.gr <- parc.c[[community]][retain.dims]
    do.call(rbind, mclapply(1:dim(data)[1], function(i) {
      dmeas <- data[i,]
      signal <- median(dmeas[cmp.gr == 1]) - median(dmeas[cmp.gr == 0])
      test.res <- wilcox.test(dmeas[cmp.gr == 1], dmeas[cmp.gr == 0], alternative = "greater")
      return(data.frame(Dataset=cov.dat[i,"Dataset"], Subid=cov.dat[i,"Subid"], Session=cov.dat[i,"Session"],
                        Signal=signal, Effect.Size=test.res$statistic, p.value=test.res$p.value,
                        Community=community, Parcellation=parcellation))
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
