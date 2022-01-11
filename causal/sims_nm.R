require('./causalComBat.R')

sim_nosupp <- function(n=100, pi=.5, null=FALSE) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  xs <- sapply(batches, function(batch) {
    ifelse(batch==0, runif(1, min=0, max=1), -runif(1, min=0, max=1))
  })
  
  eps <- rnorm(n=n, mean=0, sd=1)
  
  ys <- sapply(xs, function(x) {
    ifelse(x > 0, x, 2*x)
  }) + eps
  
  if (!null) {
    ys <- ys + batches
  }
  
  return(list(Y=matrix(ys, ncol=1), Batch=matrix(batches, ncol=1), X=matrix(xs, ncol=1),
              Eps=matrix(eps, ncol=1)))
}

sim_no_ov <- function(n=100, pi=.5, null=FALSE) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  xs <- sapply(batches, function(batch) {
    ifelse(batch==0, rbeta(1, 2, 8), rbeta(1, 8, 2))
  })
  
  eps <- 1/2*rnorm(n=n, mean=0, sd=1)
  
  ys <- sapply(xs, function(x) {
    ifelse(x > 0, (x - .5)^2, 2*(x - .5)^2)
  }) + eps
  
  if (!null) {
    ys <- ys + batches
  }
  
  
  return(list(Y=matrix(ys, ncol=1), Batch=matrix(batches, ncol=1), X=matrix(xs, ncol=1),
              Eps=matrix(eps, ncol=1)))
}

sim_sig <- function(n=100, pi=.5, null=FALSE) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  xs <- runif(n)
  
  eps <- rnorm(n=n, mean=0, sd=1)
  
  ys <- sapply(xs, function(x) {
    ifelse(x > 0, (x - .5)^2, 2*(x - .5)^2)
  }) + eps
  
  if (!null) {
    ys <- ys + batches
  }
  
  return(list(Y=matrix(ys, ncol=1), Batch=matrix(batches, ncol=1), X=matrix(xs, ncol=1),
              Eps=matrix(eps, ncol=1)))
}

nrep <- 200
ncores <- parallel::detectCores() - 1

sim = list("No Support"=sim_nosupp, "No Overlap" = sim_no_ov, "Signal"=sim_sig)

results <- do.call(rbind, lapply(names(sim), function(name) {
  do.call(rbind, lapply(c(TRUE, FALSE), function(null) {
    do.call(rbind, mclapply(1:nrep, function(i) {
      sim <- do.call(sim[[name]], list(null=null))
      
      testcd <- cdcov.test(sim$Y, sim$Batch, sim$X, num.bootstrap=1000)$p.value
      testcausal <- tryCatch(causal.cdcov(sim$Y, sim$Batch, data.frame(X=sim$X), match.form="X",
                                          match.args=list(method="nearest", exact=NULL, replace=FALSE, caliper=.1))$Test$p.value,
                             error=function(e) {
                               return(NaN)
                             })
      
      return(data.frame(Simulation=name, Null=ifelse(null, "Null", "Signal"), i=i, p.value=c(testcd, testcausal),
                        approach=c("Conditional", "Causal")))
    }, mc.cores=ncores))
  }))
}))
