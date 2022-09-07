
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
    ifelse(batch==0, rbeta(1, 2, 7), rbeta(1, 7, 2))
  })
  
  eps <- 1/2*rnorm(n=n, mean=0, sd=1)
  
  ys <- -sapply(xs, function(x) {
    ifelse(x > 0, 2*(x - .5), 2*(x - .5))
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
  
  eps <- 1/2*rnorm(n=n, mean=0, sd=1)
  
  ys <- -sapply(xs, function(x) {
    ifelse(x > 0, 2*(x - .5), 2*(x - .5))
  }) + eps
  
  if (!null) {
    ys <- ys + batches
  }
  
  return(list(Y=matrix(ys, ncol=1), Batch=matrix(batches, ncol=1), X=matrix(xs, ncol=1),
              Eps=matrix(eps, ncol=1)))
}
