# Causal ComBat
#
# A function for performing ComBat across multiple datasets, achieved via
# matching.
# Inputs
#    X: An n row matrix (observations) with d columns (features)
#    batches: A n vector (observations) of batch annotations
#    covariates: A n row matrix (observations) with d columns (covariates)
#    match.form: the matching formula to use for matching individuals across
#       datasets by the covariates.
#    match.args: hyper parameters for performing matching of individuals
#       across batches (passed directly to matchit)
require(MatchIt)
require(sva)
require(cdcsis)
require(tidyverse)
require(nnet)

causal.ComBat <- function(X, batches, covariates, match.form, match.args=list(method="nearest", exact=NULL, replace=FALSE, caliper=.1)) {
  retain.ids <- unique(do.call(match_batches, list(batches, covariates, match.form, match.args=match.args)))
  X.tilde <- X[retain.ids,]; Y.tilde <- covariates[retain.ids,]; t.tilde <- batches[retain.ids]
  
  mod <- model.matrix(as.formula(sprintf("~%s", match.form)), data=Y.tilde)
  dat.norm <- t(ComBat(t(X.tilde), t.tilde, mod = mod))
  return(list(Data=dat.norm,
              Batches=t.tilde,
              Covariates=Y.tilde,
              Retained.Ids=retain.ids))
}

# Propensity Trim across multiple exposures via vector matching
# adapted from Lopez et al. 2017
#
# Ts denotes (1 of K) categorical treatments in an n (samples) vector
# Xs is a n x d matrix of covariates (columns) for each sample (rows)
#
# returns a boolean array for whether to include/exclude samples from
# subsequent analysis
vm_trim <- function(Ts, covariates) {
  covariates = as.data.frame(covariates)
  
  # Fitting the Multinomial Logistic Regression Model
  K <- length(unique(Ts))
  Ts <- as.numeric(factor(Ts, levels=unique(Ts)))
  Ts_unique = unique(Ts)

  m <- nnet::multinom(factor(Ts) ~ ., data = as.data.frame(covariates))
  
  # Making predictions using the fitted model
  pred <- predict(m, newdata = as.data.frame(covariates), type = "probs")
  
  # if only binary treatment levels, add a column for the reference
  if (K == 2) {
    pred <- cbind(1 - pred, pred)
    colnames(pred) <- Ts_unique
  }
  
  # Function to calculate the range of predicted probabilities for each treatment
  calculate_Rtable <- function(Tval) {
    Rtab = matrix(0, nrow=2, ncol=length(Ts_unique))
    for (Tp in Ts_unique) {
      Rtab[, Tp] = c(min(pred[Ts == Tp, Tval]), max(pred[Ts == Tp, Tval]))
    }
    c(max(Rtab[1,]), min(Rtab[2,]))
  }
  
  # Creating the Rtable to store the range of predicted probabilities for each treatment
  Rtable <- t(sapply(Ts_unique, calculate_Rtable))
  rownames(Rtable) = Ts_unique
  
  # Function to check if each observation satisfies balance condition for each treatment
  check_balance <- function(i) {
    sapply(as.character(Ts_unique), function(Tval) pred[i, Tval] >= Rtable[Tval, 1] & pred[i, Tval] <= Rtable[Tval, 2])
  }
  
  # Creating the balance_check matrix to check if each observation satisfies balance condition for each treatment
  balance_check <- t(sapply(1:nrow(covariates), check_balance))
  
  # Finding observations that satisfy balance condition for all treatments
  balanced_ids <- apply(balance_check, 1, all)
  return(balanced_ids)
}

zero_one_dist <- function(Ts) {
  # Convert input vector to a factor
  Ts_fact <- factor(Ts)
  
  # Perform one-hot encoding using model.matrix
  encoded_matrix <- model.matrix(~Ts_fact - 1)
  
  # Convert the matrix to a data frame (optional)
  Ts_ohe <- as.data.frame(encoded_matrix)
  
   return(dist(Ts_ohe))
}

# Conditional distance correlation from Wang et al., 2015
cond.dcorr <- function(X, Ts, covariates, R=1000, dist.method="euclidean", distance = FALSE, seed=1, num.threads=1) {
  covariates <- as.data.frame(covariates)
  
  # vector match for propensity trimming, and then reduce sub-sample to the
  # propensity matched subset
  if (length(retain.ids) == 0) {
    stop("No samples remain after balancing.")
  }
  if (isTRUE(distance)) {
    DX <- as.dist(X)
  } else {
    DX = dist(X, method=dist.method)
  }
  
  DT = zero_one_dist(Ts)
  
  # run statistical test
  test.out <- cdcov.test(DX, DT, covariates, num.bootstrap = R,
                         seed=seed, num.threads=num.threads, distance=TRUE)
  return(list(Test=test.out))
}

# R implementation of Bridgeford et al., 2023
causal.cdcorr <- function(X, Ts, covariates, R=1000, dist.method="euclidean", distance = FALSE, seed=1, num.threads=1) {
  covariates <- as.data.frame(covariates)
  
  # vector match for propensity trimming, and then reduce sub-sample to the
  # propensity matched subset
  retain.ids <- which(vm_trim(Ts, covariates))
  if (length(retain.ids) == 0) {
    stop("No samples remain after balancing.")
  }
  if (isTRUE(distance)) {
    DX.tilde <- as.dist(as.matrix(X)[retain.ids, retain.ids])
  } else {
    X.tilde <- X[retain.ids,]
    DX.tilde = dist(X.tilde, method=dist.method)
  }
  Y.tilde <- covariates[retain.ids,]
  
  DT.tilde = zero_one_dist(Ts[retain.ids])
  
  # run statistical test
  test.out <- cdcov.test(DX.tilde, DT.tilde, Y.tilde, num.bootstrap = R,
                        seed=seed, num.threads=num.threads, distance=TRUE)
  return(list(Test=test.out,
              Retained.Ids=retain.ids)) 
}

match_batches <- function(batches, covariates, match.form, match.args=NULL) {
  # obtain the smallest batch
  batches <- as.character(batches)
  covariates <- cbind(data.frame(Batch=batches), covariates)
  batch.sum <- batches %>% table()
  batch.names <- names(batch.sum)
  tx.batch <- batch.names[which.min(batch.sum)]
  rownames(covariates) <- 1:nrow(covariates)
  covar.tx <- covariates[batches == tx.batch,,drop=FALSE]
  
  
  paired.matches <- lapply(batch.names[batch.names != tx.batch], function(batch) {
    covar.cont <- covariates[batches == batch,]
    covariate.match(covar.tx, covar.cont, match.form=match.form, match.args=match.args)
  })
  
  I.mat <- which(apply(sapply(paired.matches, function(x) x$I.mat.k), c(1), sum) > 0)
  M.mat <- apply(
    do.call(cbind, lapply(paired.matches, function(x) x$M.mat.k))[I.mat,],
    c(2), sum)
  retain.ids <- as.numeric(c(names(I.mat), names(which(M.mat != 0))))
  return(retain.ids)
}

covariate.match <- function(covar.tx, covar.cont, match.form, match.args=NULL) {
  n.kprime <- dim(covar.tx)[1]; n.k <- dim(covar.cont)[1]
  n.matches <- floor(n.k/n.kprime)
  match <- do.call(matchit, c(list(formula(sprintf("as.factor(Treatment) ~ %s", match.form)),
                   data=rbind(covar.tx %>% mutate(Treatment = 1),
                              covar.cont %>% mutate(Treatment = 0)),
                   ratio=n.matches), match.args))
  mat.mtx <- match$match.matrix
  I.mat.k <- as.numeric(apply(mat.mtx, c(1), function(x) {sum(!is.na(x))}) > 0)
  names(I.mat.k) <- names(match$weights[1:n.kprime])
  M.mat.k <- matrix(FALSE, nrow=n.kprime, ncol=n.k)
  ctrl.names <- names(match$weights[(n.kprime + 1):(n.kprime + n.k)])
  for (i in 1:n.k) {
    for (j in 1:n.kprime) {
      M.mat.k[j,i] <- ifelse(ctrl.names[i] %in% mat.mtx[j,], 1, 0)
    }
  }
  colnames(M.mat.k) <- as.numeric(ctrl.names)
  
  return(list(I.mat.k=I.mat.k, M.mat.k=M.mat.k))
}

require(reticulate)
use_virtualenv("/opt/neuroharm/", required=TRUE)
py_config()
neuroharm <- import("neuroHarmonize")
causal.NeuroH <- function(X, batches, covariates, match.form, match.args=list(method="nearest", exact=NULL, replace=FALSE, caliper=.1)) {
  retain.ids <- unique(balance.batches(batches, covariates, match.form, exact=exact))
  X.tilde <- X[retain.ids,]; Y.tilde <- covariates[retain.ids,]; t.tilde <- batches[retain.ids]
  
  covars <- data.frame(SITE=t.tilde, AGE=Y.tilde$Age, SEX_M=as.numeric(Y.tilde$Sex == 2))
  dat.norm <- neuroharm$harmonizationLearn(X.tilde, covars)[[2]]
  return(list(Data=dat.norm,
              Batches=t.tilde,
              Covariates=Y.tilde,
              Retained.Ids=retain.ids))
}