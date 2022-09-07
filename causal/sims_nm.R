source('./causalComBat.R')
require(tidyverse)
require(parallel)
require(cdcsis)
source('./sim_helpers.R')

nrep <- 200
ncores <- parallel::detectCores() - 1

sims = list(#"No Support"=sim_nosupp, 
  "No Overlap" = sim_no_ov, "Overlap"=sim_sig)

results <- do.call(rbind, lapply(names(sims), function(name) {
  do.call(rbind, lapply(c(TRUE, FALSE), function(null) {
    do.call(rbind, mclapply(1:nrep, function(i) {
      print(i)
      sim <- do.call(sims[[name]], list(null=null))
      
      testcd <- cdcov.test(sim$Y, sim$Batch, sim$X, num.bootstrap=1000)$p.value
      testcausal <- tryCatch(causal.cdcov(sim$Y, sim$Batch, data.frame(X=sim$X), match.form="X",
                                          match.args=list(method="nearest", exact=NULL, replace=FALSE, caliper=.1))$Test$p.value,
                             error=function(e) {
                               return(NaN)
                             })
      
      return(data.frame(Setting=name, Null=ifelse(null, "Null", "Signal"), i=i, p.value=c(testcd, testcausal),
                        Approach=c("Conditional", "Causal")))
    }, mc.cores=ncores))
  }))
}))

results$p.adj[!is.nan(results$p.value)] <- p.adjust(results$p.value[!is.nan(results$p.value)], method="BH")

saveRDS(results, '../data/dcorr/sim_results.rds')