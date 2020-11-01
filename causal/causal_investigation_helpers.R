
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

causal_ana_site <- function(Dmtx.dat, cov.dat, trim=.01, R=500) {
  datasets <- sort(unique((cov.dat %>% dplyr::select(Dataset))$Dataset))
  dset.pairs <- combn(datasets, 2)
  result <- do.call(rbind, mclapply(1:dim(dset.pairs)[2], function(x) {
    tryCatch({
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
      Dmtx.dat.ij <- Dmtx.dat[cov.dat.ij$id, cov.dat.ij$id]
      ## Uncorrected
      test.uncor <- pdcor.test(as.dist(Dmtx.dat.ij), cov.dat.ij$Treatment,
                               cov.dat.ij %>% dplyr::select(Continent, Sex, Age) %>%
                                 mutate(Continent=as.numeric(Continent), Sex=as.numeric(Sex),
                                        Age=as.numeric(Age)), R=R)
      result$uncor <- data.frame(Data="Untrimmed", Method="PDcor", Dataset.Trt=dset.i,
                                 Dataset.Ctrl=dset.j, Effect=test.uncor$estimate,
                                 PValue=test.uncor$p.value)

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
                                    Dataset.Ctrl=dset.j, Effect=test.trim$estimate,
                                    PValue=test.trim$p.value)
        }
      }
      return(do.call(rbind, result) %>%
               mutate(Effect.Name="Site"))
    }, error=function(e) {
      return(data.frame(Data=c("Untrimmed", "Trimmed"), Method="PDcor", Dataset.Trt=dset.i,
                        Dataset.Ctrl=dset.j, Effect.Name="Site", Effect=NA, PValue=NA))
    })
  }, mc.cores = ncores))
  return(result)
}

compute_effect <- function(D.i, cov.dat.i, E1.name, E2.name, R=1000, nboots=100, ...) {
  n.i <- length(cov.dat.i$id)

  test.uncor <- pdcor.test(D.i, as.numeric(cov.dat.i[[E1.name]]),
                           as.numeric(cov.dat.i[[E2.name]]), R=R)
  ## Bootstrapped CIs
  bootstr.ci <- sapply(1:nboots, function(i) {
    ids <- sample(1:n.i, size=n.i, replace = TRUE)
    cov.dat.boot <- cov.dat.i[ids,]
    Dmtx.dat.boot <- as.dist(as.matrix(D.i)[ids, ids])
    test.boot <- pdcor(Dmtx.dat.boot, as.numeric(cov.dat.boot[[E1.name]]),
                       as.numeric(cov.dat.boot[[E2.name]]))
    return(test.boot)
  })
  ci.npboot <- quantile(bootstr.ci, c(.025, .975))

  jk.stat <- sapply(1:n.i, function(i) {
    ids <- (1:n.i)[-i]
    cov.dat.boot <- cov.dat.i[ids,]
    Dmtx.dat.boot <- as.dist(as.matrix(D.i)[ids, ids])
    test.boot <- pdcor(Dmtx.dat.boot, as.numeric(cov.dat.boot[[E1.name]]),
                       as.numeric(cov.dat.boot[[E2.name]]))
    return(test.boot)
  })
  psi <- test.uncor$estimate*n.i - (n.i - 1)*jk.stat
  ps.mean <- mean(psi); ps.var <- var(psi)
  ci.jk <- c(ps.mean + qnorm(.025)*sqrt(1/n.i*ps.var),
             ps.mean + qnorm(.975)*sqrt(1/n.i*ps.var))

  return(data.frame(Method="PDcor", Effect.Name=E1.name, Continent=continent, Effect=res$estimate, Effect.lwr.jk=ci.jk[1], Effect.upr.jk=ci.jk[2],
                    Effect.lwr.npboot=ci.npboot[1], Effect.upr.npboot=ci.npboot[2], Effect.lwr.adjboot=ci.adjboot[1], Effect.upr.adjboot=ci.adjboot[2],
                    PValue=test.uncor$p.value, Entropy=entropy(as.numeric(cov.dat.i[[E1.name]]), method="MM"), Variance=var(as.numeric(cov.dat.i[[E1.name]])),
                    n=n.i, N=length(unique(cov.dat.i$Subid))))
}

causal_ana_cov <- function(Dmtx.dat, cov.dat, nboots=100, R=1000) {
  cov.dat <- cov.dat %>%
    ungroup() %>%
    mutate(id=row_number())
  datasets <- sort(unique((cov.dat %>% dplyr::select(Dataset))$Dataset))
  result.sex <- do.call(rbind, mclapply(datasets, function(dataset) {
    tryCatch({
      cov.dat.dset <- cov.dat %>%
        filter(Dataset == dataset)

      n.i <- length(cov.dat.dset$id)
      Dmtx.dat.dset <- as.dist(Dmtx.dat[cov.dat.dset$id, cov.dat.dset$id])
      res <- compute_effect(Dmtx.dat.dset, cov.dat.dset, "Sex", "Age", R=R, nboots=nboots)

      return(res)
    }, error=function(e) {
      print(e)
      return(NULL)
    })

  }, mc.cores=ncores))

  result.age <- do.call(rbind, mclapply(datasets, function(dataset) {
    tryCatch({
      cov.dat.dset <- cov.dat %>%
        filter(Dataset == dataset)

      n.i <- length(cov.dat.dset$id)
      Dmtx.dat.dset <- as.dist(Dmtx.dat[cov.dat.dset$id, cov.dat.dset$id])
      res <- compute_effect(Dmtx.dat.dset, cov.dat.dset, "Age", "Sex", R=R, nboots=nboots)

      return(res)
    }, error=function(e) {
      return(NULL)
    })
  }, mc.cores=ncores))
  return(rbind(result.sex, result.age))
}

causal_ana_cov_cont <- function(Dmtx.dat, cov.dat, nboots=100, R=1000) {
  cov.dat <- cov.dat %>%
    ungroup() %>%
    mutate(id=row_number())
  continents <- sort(unique((cov.dat %>% dplyr::select(Continent))$Continent))
  result.sex <- do.call(rbind, mclapply(continents, function(continent) {
    tryCatch({
      cov.dat.dset <- cov.dat %>%
        filter(Continent == continent)

      n.i <- length(cov.dat.dset$id)
      Dmtx.dat.dset <- as.dist(Dmtx.dat[cov.dat.dset$id, cov.dat.dset$id])
      res <- compute_effect(Dmtx.dat.dset, cov.dat.dset, "Sex", "Age", R=R, nboots=nboots)

      return(res)
    }, error=function(e) {
      return(NULL)
    })
  }, mc.cores=ncores))

  result.age <- do.call(rbind, mclapply(continents, function(continent) {
    tryCatch({
      cov.dat.dset <- cov.dat %>%
        filter(Continent == continent)

      n.i <- length(cov.dat.dset$id)
      Dmtx.dat.dset <- as.dist(Dmtx.dat[cov.dat.dset$id, cov.dat.dset$id])
      res <- compute_effect(Dmtx.dat.dset, cov.dat.dset, "Age", "Sex", R=R, nboots=nboots)

      return(res)
    }, error=function(e) {
      return(NULL)
    })
  }, mc.cores=ncores))
  return(rbind(result.sex, result.age))
}
