require(tidyverse)
require(energy)
require(parallel)

ncores <- parallel::detectCores() - 1

mgc.sims.2sphere <- function(n, d, r, cov.scale=0) {
  u <- matrix(mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=diag(d)), nrow=n)
  unorm <- diag(sqrt(apply(u^2, 1, sum)))
  pts <- r*(ginv(unorm) %*% u)
  pts <- pts + matrix(mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=cov.scale*diag(d)), nrow=n)
  return(pts)
}

associational.sim <- function(n, d, eps, s) {
  Y <- sample(c(0, 1), n, prob=c(0.5, 0.5), replace=TRUE) # the site assignment vector
  
  # whether the individual is from sphere or unit 2-ball
  which.distn <- sample(c(0, 1), n, prob=c(eps, (1 - eps)), replace=TRUE)
  
  gaus <- t(matrix(mvrnorm(n=n, mu=array(0, dim=c(d, 1)), Sigma=diag(d)), nrow=d))
  sphere <- mgc.sims.2sphere(n, d, r=4, cov.scale=0.1)
  
  X <- which.distn*gaus + (1 - which.distn)*sphere + Y*s/sqrt(d)
  return(list(X=X, Y=Y))
}

mvlm.associational <- function(X, Y, ...) {
  if((dim(X)[1] - 1) < dim(X)[2]) {
    X <- prcomp(X, center=TRUE, scale=TRUE)$x[, 1:(dim(X)[1] - 2)]
  }
  man.res <- summary(manova(X ~ Y))
  p.val <- man.res$stats["Y", "Pr(>F)"]
  return(p.val)
}

dcorr.associational <- function(X, Y, R=1000) {
  p.val <- dcor.test(X, Y, R=1000)$p.value
  return(p.val)
}

n <- 100
nbreaks <- 10
ds <- round(2^seq(log2(2), log2(1000), length.out=nbreaks))
d.fixed <- 101
epss <- seq(0, 1, length.out=nbreaks)
eps.fixed <- 0.5
ss <- seq(0, 1, length.out=nbreaks)*2
s.fixed <- 1
R <- 1000
nrep <- 100

alpha <- .05

methods <- list("MANOVA"=mvlm.associational, "DCorr"=dcorr.associational)

res.ss <- do.call(rbind, lapply(ss, function(s) {
  do.call(rbind, lapply(1:nrep, function(i) {
    sim <- associational.sim(n=n, d=d.fixed, eps=eps.fixed, s=s)
    X <- sim$X; Y <- sim$Y
    do.call(rbind, lapply(names(methods), function(method) {
      pval <- do.call(methods[[method]], list(X=X, Y=Y, R=R))
      return(data.frame(s=s, d=d.fixed, eps=eps.fixed, Method=method, i=i, p.value=pval,
                        Simulation="Associational", Setting="Effect Size"))
    }))
  }))#, mc.cores=ncores))
}))

res.ds <- do.call(rbind, lapply(ds, function(d) {
  do.call(rbind, mclapply(1:nrep, function(i) {
    sim <- associational.sim(n=n, d=d, eps=eps.fixed, s=s)
    X <- sim$X; Y <- sim$Y
    do.call(rbind, lapply(names(methods), function(method) {
      pval <- do.call(methods[[method]], list(X=X, Y=Y, R=R))
      return(data.frame(s=s.fixed, d=d, eps=eps.fixed, Method=method, i=i, p.value=pval,
                        Simulation="Associational", Setting="Dimensionality"))
    }))
  }))#, mc.cores=ncores))
}))

res.epss <- do.call(rbind, lapply(epss, function(eps) {
  do.call(rbind, mclapply(1:nrep, function(i) {
    sim <- associational.sim(n=n, d=d.fixed, eps=eps, s=s.fixed)
    X <- sim$X; Y <- sim$Y
    do.call(rbind, lapply(names(methods), function(method) {
      pval <- do.call(methods[[method]], list(X=X, Y=Y, R=R))
      return(data.frame(s=s.fixed, d=d.fixed, eps=eps, Method=method, i=i, p.value=pval,
                        Simulation="Associational", Setting="Non-Parametricity"))
    }))
  }))#, mc.cores=ncores))
}))

saveRDS(rbind(res.ss, res.ds, res.eps), file="./associational.rds")