require(tidyverse)
require(energy)
require(parallel)
require(MASS)
require(irlba)
require(cdcsis)

ncores <- parallel::detectCores() - 1


n <- 100
nbreaks <- 10
ds <- round(2^seq(log2(2), log2(1000), length.out=nbreaks))
d.fixed.high <- 101
d.fixed.low <- 20
epss <- seq(0, 1, length.out=nbreaks)
eps.fixed <- 0.5
ss <- seq(0, 1, length.out=nbreaks)*2
s.fixed <- 1
pis <- seq(0, 1, length.out=nbreaks)
pi.fixed <- 0.75
R <- 1000
nrep <- 100

alpha <- .05

## Associational Simulations ##
print("Associational...")
methods.ass <- list("MANOVA"=mvlm.associational, "DCorr"=dcorr.associational)

print("Effect Size...")
res.ss.ass <- do.call(rbind, lapply(ss, function(s) {
  print(s)
  do.call(rbind, lapply(c(d.fixed.low, d.fixed.high), function(d.fixed) {
    do.call(rbind, mclapply(1:nrep, function(i) {
      sim <- associational.sim(n=n, d=d.fixed, eps=eps.fixed, s=s)
      X <- sim$X; Y <- sim$Y
      do.call(rbind, lapply(names(methods.ass), function(method) {
        pval <- do.call(methods.ass[[method]], list(X=X, Y=Y, R=R))
        return(data.frame(s=s, d=d.fixed, eps=eps.fixed, Method=method, i=i, p.value=pval,
                          Simulation="Associational", Setting=sprintf("Effect Size, D=%d", d.fixed)))
      }))
    }, mc.cores=ncores))
  }))
}))

print("Dimensionality...")
res.ds.ass <- do.call(rbind, lapply(ds, function(d) {
  print(d)
  do.call(rbind, mclapply(1:nrep, function(i) {
    sim <- associational.sim(n=n, d=d, eps=eps.fixed, s=s.fixed)
    X <- sim$X; Y <- sim$Y
    do.call(rbind, lapply(names(methods.ass), function(method) {
      pval <- do.call(methods.ass[[method]], list(X=X, Y=Y, R=R))
      return(data.frame(s=s.fixed, d=d, eps=eps.fixed, Method=method, i=i, p.value=pval,
                        Simulation="Associational", Setting="Dimensionality"))
    }))
  }, mc.cores=ncores))
}))

print("Non-Parametricity...")
res.epss.ass <- do.call(rbind, lapply(epss, function(eps) {
  print(eps)
  do.call(rbind, lapply(c(d.fixed.low, d.fixed.high), function(d.fixed) {
    do.call(rbind, mclapply(1:nrep, function(i) {
      sim <- associational.sim(n=n, d=d.fixed, eps=eps, s=s.fixed)
      X <- sim$X; Y <- sim$Y
      do.call(rbind, lapply(names(methods.ass), function(method) {
        pval <- do.call(methods.ass[[method]], list(X=X, Y=Y, R=R))
        return(data.frame(s=s.fixed, d=d.fixed, eps=eps, Method=method, i=i, p.value=pval,
                          Simulation="Associational", Setting="Non-Parametricity, D=%d", d.fixed))
      }))
    }, mc.cores=ncores))
  }))
}))

## Conditional Simulations ##
print("Conditional Simulations...")
methods.cond <- list("MANOVA"=mvlm.cond, "CDCov"=cdcov.cond, "DCorr"=dcorr.associational)

print("Effect Size...")
res.ss.cond <- do.call(rbind, lapply(ss, function(s) {
  print(s)
  do.call(rbind, lapply(c(d.fixed.low, d.fixed.high), function(d.fixed) {
    do.call(rbind, mclapply(1:nrep, function(i) {
      sim <- cond.sim(n=n, d=d.fixed, eps=eps.fixed, s=s, pi=0.5)
      X <- sim$X; Y <- sim$Y; Z <- sim$Z
      do.call(rbind, lapply(names(methods.cond), function(method) {
        pval <- do.call(methods.cond[[method]], list(X=X, Y=Y, Z=Z, R=R))
        return(data.frame(s=s, d=d.fixed, eps=eps.fixed, Method=method, i=i, p.value=pval,
                          Simulation="Conditional", Setting=sprintf("Effect Size, D=%d", d.fixed)))
      }))
    }, mc.cores=ncores))
  }))
}))

print("Dimensionality...")
res.ds.cond <- do.call(rbind, lapply(ds, function(d) {
  print(d)
  do.call(rbind, mclapply(1:nrep, function(i) {
    sim <- cond.sim(n=n, d=d, eps=eps.fixed, s=s.fixed, pi=0.5)
    X <- sim$X; Y <- sim$Y; Z <- sim$Z
    do.call(rbind, lapply(names(methods.cond), function(method) {
      pval <- do.call(methods.cond[[method]], list(X=X, Y=Y, Z=Z, R=R))
      return(data.frame(s=s.fixed, d=d, eps=eps.fixed, Method=method, i=i, p.value=pval,
                        Simulation="Conditional", Setting="Dimensionality"))
    }))
  }, mc.cores=ncores))
}))

print("Non-Parametricity...")
res.epss.cond <- do.call(rbind, lapply(epss, function(eps) {
  print(eps)
  do.call(rbind, lapply(c(d.fixed.low, d.fixed.high), function(d.fixed) {
    do.call(rbind, mclapply(1:nrep, function(i) {
      sim <- cond.sim(n=n, d=d.fixed, eps=eps, s=s.fixed, pi=0.5)
      X <- sim$X; Y <- sim$Y; Z <- sim$Z
      do.call(rbind, lapply(names(methods.cond), function(method) {
        pval <- do.call(methods.cond[[method]], list(X=X, Y=Y, Z=Z, R=R))
        return(data.frame(s=s.fixed, d=d.fixed, eps=eps, Method=method, i=i, p.value=pval,
                          Simulation="Conditional", Setting="Non-Parametricity, D=%d", d.fixed))
      }))
    }, mc.cores=ncores))
  }))
}))


## Causal Simulations ##
print("Causal Simulations...")
methods.causal <- list("MANOVA"=mvlm.causal, "Causal CDCov"=cdcov.causal, "CDCov"=cdcov.cond, "DCorr"=dcorr.associational)

print("Effect Size...")
res.ss.caus <- do.call(rbind, lapply(ss, function(s) {
  print(s)
  do.call(rbind, lapply(c(d.fixed.low, d.fixed.high), function(d.fixed) {
    do.call(rbind, mclapply(1:nrep, function(i) {
      sim <- cond.sim(n=n, d=d.fixed, eps=eps.fixed, s=s, pi=pi.fixed, causal=TRUE)
      X <- sim$X; Y <- sim$Y
      do.call(rbind, lapply(names(methods.causal), function(method) {
        pval <- do.call(methods[[method]], list(X=X, Y=Y, Z=Z, R=R))
        return(data.frame(s=s, d=d.fixed, eps=eps.fixed, pi=pi.fixed, Method=method, i=i, p.value=pval,
                          Simulation="Causal", Setting=sprintf("Effect Size, D=%d", d.fixed)))
      }))
    }, mc.cores=ncores))
  }))
}))

print("Dimensionality...")
res.ds.caus <- do.call(rbind, lapply(ds, function(d) {
  print(d)
  do.call(rbind, mclapply(1:nrep, function(i) {
    sim <- associational.sim(n=n, d=d, eps=eps.fixed, s=s.fixed, pi=pi.fixed, causal=TRUE)
    X <- sim$X; Y <- sim$Y
    do.call(rbind, lapply(names(methods.causal), function(method) {
      pval <- do.call(methods[[method]], list(X=X, Y=Y, Z=Z, R=R))
      return(data.frame(s=s.fixed, d=d, eps=eps.fixed, pi=pi.fixed, Method=method, i=i, p.value=pval,
                        Simulation="Causal", Setting="Dimensionality"))
    }))
  }, mc.cores=ncores))
}))

print("Non-Parametricity...")
res.epss.caus <- do.call(rbind, lapply(epss, function(eps) {
  print(eps)
  do.call(rbind, lapply(c(d.fixed.low, d.fixed.high), function(d.fixed) {
    do.call(rbind, mclapply(1:nrep, function(i) {
      sim <- associational.sim(n=n, d=d.fixed, eps=eps, s=s.fixed, pi=pi.fixed, causal=TRUE)
      X <- sim$X; Y <- sim$Y
      do.call(rbind, lapply(names(methods.causal), function(method) {
        pval <- do.call(methods[[method]], list(X=X, Y=Y, Z=Z, R=R))
        return(data.frame(s=s.fixed, d=d.fixed, eps=eps, pi=pi.fixed, Method=method, i=i, p.value=pval,
                          Simulation="Causal", Setting="Non-Parametricity, D=%d", d.fixed))
      }))
    }, mc.cores=ncores))
  }))
}))

print("Irregularity...")
res.epss.caus <- do.call(rbind, lapply(pis, function(pi) {
  print(eps)
  do.call(rbind, lapply(c(d.fixed.low, d.fixed.high), function(d.fixed) {
    do.call(rbind, mclapply(1:nrep, function(i) {
      sim <- associational.sim(n=n, d=d.fixed, eps=eps.fixed, s=s.fixed, pi=pi, causal=TRUE)
      X <- sim$X; Y <- sim$Y
      do.call(rbind, lapply(names(methods.causal), function(method) {
        pval <- do.call(methods[[method]], list(X=X, Y=Y, Z=Z, R=R))
        return(data.frame(s=s.fixed, d=d.fixed, eps=eps.fixed, pi=pi, Method=method, i=i, p.value=pval,
                          Simulation="Causal", Setting="Irregularity, D=%d", d.fixed))
      }))
    }, mc.cores=ncores))
  }))
}))

saveRDS(list(Associational=rbind(res.ss.ass, res.ds.ass, res.eps.ass),
             Conditional=rbind(res.ss.cond, res.ds.cond, res.eps.cond),
             Causal=rbind(res.ss.caus, res.ds.caus, res.eps.caus)), file="./sim_results.rds")