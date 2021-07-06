# Causal ComBat
#
# A function for performing ComBat across multiple datasets.
# Inputs
#    X: An n row matrix (observations) with d columns (features)
#     

#
require(MatchIt)
require(sva)
use_virtualenv("/opt/neuroharm/", required=TRUE)
py_config()
neuroharm <- import("neuroHarmonize")

causal.ComBat <- function(X, batches, covariates, match.form, exact=NULL) {
  retain.ids <- unique(balance.batches(batches, covariates, match.form, exact=exact))
  X.tilde <- X[retain.ids,]; Y.tilde <- covariates[retain.ids,]; t.tilde <- batches[retain.ids]
  
  mod <- model.matrix(as.formula(sprintf("~%s", match.form)), data=Y.tilde)
  dat.norm <- t(ComBat(t(X.tilde), t.tilde, mod = mod))
  return(list(Data=dat.norm,
              Batches=t.tilde,
              Covariates=Y.tilde,
              Retained.Ids=retain.ids))
}

causal.NeuroH <- function(X, batches, covariates, match.form, exact=NULL) {
  retain.ids <- unique(balance.batches(batches, covariates, match.form, exact=exact))
  X.tilde <- X[retain.ids,]; Y.tilde <- covariates[retain.ids,]; t.tilde <- batches[retain.ids]
  
  covars <- data.frame(SITE=t.tilde, AGE=Y.tilde$Age, SEX_M=as.numeric(Y.tilde$Sex == 2))
  dat.norm <- neuroharm$harmonizationLearn(X.tilde, covars)[[2]]
  return(list(Data=dat.norm,
              Batches=t.tilde,
              Covariates=Y.tilde,
              Retained.Ids=retain.ids))
}

balance.batches <- function(batches, covariates, match.form, exact=NULL) {
  # obtain the smallest batch
  batch.sum <- batches %>% table()
  batch.names <- names(batch.sum)
  tx.batch <- batch.names[which.min(batch.sum)]
  rownames(covariates) <- 1:nrow(covariates)
  covar.tx <- covariates[batches == tx.batch,]
  
  
  paired.matches <- lapply(batch.names[batch.names != tx.batch], function(batch) {
    covar.cont <- covariates[batches == batch,]
    covariate.match(covar.tx, covar.cont, match.form=match.form, exact=exact)
  })
  
  I.mat <- which(apply(sapply(paired.matches, function(x) x$I.mat.k), c(1), sum) > 0)
  M.mat <- apply(
    do.call(cbind, lapply(paired.matches, function(x) x$M.mat.k))[I.mat,],
    c(2), sum)
  retain.ids <- as.numeric(c(names(I.mat), names(which(M.mat != 0))))
  return(retain.ids)
}

covariate.match <- function(covar.tx, covar.cont, match.form, exact=NULL) {
  n.kprime <- dim(covar.tx)[1]; n.k <- dim(covar.cont)[1]
  n.matches <- floor(n.k/n.kprime)
  match <- matchit(formula(sprintf("as.factor(Treatment) ~ %s", match.form)),
                   data=rbind(covar.tx %>% mutate(Treatment = 1),
                              covar.cont %>% mutate(Treatment = 0)),
                   method="nearest", exact=exact, ratio=n.matches, replace=FALSE,
                   caliper=.1)
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