require(tidyverse)
require(energy)
require(parallel)
require(MASS)
require(irlba)
require(cdcsis)

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

cond.sim <- function(n, d, eps, s, pi, causal=FALSE) {
  Y <- sample(c(0, 1), n, prob=c(0.5, 0.5), replace=TRUE) # the site assignment vector
  
  # whether the individual is from sphere or unit 2-ball
  which.distn <- sample(c(0, 1), n, prob=c(eps, (1 - eps)), replace=TRUE)
  
  gaus <- t(matrix(mvrnorm(n=n, mu=array(0, dim=c(d, 1)), Sigma=diag(d)), nrow=d))
  sphere <- mgc.sims.2sphere(n, d, r=4, cov.scale=0.1)
  
  pi_vec <- c(pi, (1 - pi))
  alpha_vec <- c(2, 1)
  beta_vec <- c(1, 2)
  Z1 <- sapply(Y, function(y) {sample(c(0, 1), 1, prob=c(pi_vec[as.numeric(!(as.logical(y))) + 1], pi_vec[y + 1]), replace=TRUE)})
  Z2 <- sapply(Y, function(y) {
    alpha <- alpha_vec[as.numeric(!(as.logical(y))) + 1]
    beta <- beta_vec[as.numeric(!(as.logical(y))) + 1]
    rbeta(1, alpha, beta)
  })
  
  X <- which.distn*gaus + (1 - which.distn)*sphere + Y*s/sqrt(d) + Z1*s/sqrt(d) + causal*Z2*s/sqrt(d)
  if(causal) {
    Z <- data.frame(Z1=Z1, Z2=Z2)
  } else {
    Z <- Z1
  }
  return(list(X=X, Y=Y, Z=Z))
}

getElbows <- function(dat, n = 2, threshold=FALSE, ...) {
  
  if (is.matrix(dat)) {
    d <- sort(apply(dat,2,sd), decreasing=TRUE)
  } else {
    d <- sort(dat,decreasing=TRUE)
  }
  
  if (!is.logical(threshold))
    d <- d[d > threshold]
  
  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")
  
  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  if (n > 1 && q < (p-1)) {
    q <- c(q, q + getElbows(d[(q+1):p], n-1))
  }
  
  return(q)
}

pca.dimselect <- function(X, max.dim=97, ...) {
  X.sv <- svd(scale(X))
  zg.2 <- getElbows(X.sv$d)[2]
  d <- ifelse(zg.2 < max.dim, zg.2, max.dim)
  
  X.red <- X.sv$u[,1:d] %*% diag(X.sv$d[1:d])
}

mvlm.associational <- function(X, Y, ...) {
  if((dim(X)[1] - 2) < dim(X)[2]) {
    X <- pca.dimselect(X, max.dim=dim(X)[1] - 2)
  }
  man.res <- summary(manova(X ~ Y))
  p.val <- man.res$stats["Y", "Pr(>F)"]
  return(p.val)
}

dcorr.associational <- function(X, Y, R=1000, ...) {
  p.val <- dcor.test(X, Y, R=1000)$p.value
  return(p.val)
}

mvlm.cond <- function(X, Y, Z, ...) {
  if((dim(X)[1] - 4) < dim(X)[2]) {
    X <- pca.dimselect(X, max.dim=dim(X)[1] - 3)
  }
  man.res <- summary(manova(X ~ Y + Z))
  p.val <- man.res$stats["Y", "Pr(>F)"]
}

cdcov.cond <- function(X, Y, Z, R=1000, ...) {
  p.val <- cdcov.test(X, Y, Z, num.bootstrap=1000)$p.value
  return(p.val)
}

mvlm.causal <- function(X, Y, Z, match.form="as.factor(Z1) + Z2",
                        match.args=list(method="nearest", caliper=0.1, replace=FALSE, exact="Z1"), ...) {
  retain.ids <- unique(balance.batches(batches, covariates, match.form, match.args=match.args))
  X.tilde <- X[retain.ids,]; Y.tilde <- covariates[retain.ids,]; t.tilde <- batches[retain.ids]
  
  if((dim(X.tilde)[1] - 4) < dim(X.tilde)[2]) {
    X.tilde <- pca.dimselect(X.tilde, max.dim=dim(X.tilde)[1] - 4)
  }
  man.res <- summary(manova(X.tilde ~ t.tilde + Y.tilde[,1] + Y.tilde[,2]))
  p.val <- man.res$stats["t.tilde", "Pr(>F)"]
}

cdcov.causal <- function(X, Y, Z, match.form="as.factor(Z1) + Z2",
                         match.args=list(method="nearest", caliper=0.1, replace=FALSE, exact="Z1"), R=1000, ...) {
  p.val <- causal.cdcov(X, Y, Z, match.form,  match.args, R)$Test$p.value
  return(p.val)
}