
```{r}
meth <- "Raw"

R <- 1000
print(meth)
preproc.meth <- preproc[[meth]]
X = preproc.meth$graphs; Y = preproc.meth$covariates$Sex; Z = preproc.meth$covariates$Age
DY <- dist(Y); DZ <- dist(Z)
d <- dim(X)[2]
res <- do.call(rbind, mclapply(1:4486, function(i) {
  tryCatch({
    print(i)
    DX <- dist(X[,i])
    test <- pdcor.test(DX, DY, DZ, R=R)
    return(data.frame(Edge=i, Statistic=test$statistic,
                      p.value=(1 + sum(test$replicates >= test$statistic))/(1 + R)))
  }, error=function(e) {NULL})
}, mc.cores = ncores)) %>% mutate(Method=meth)

saveRDS(res, sprintf('../data/dcorr/outputs_edgewise_%s_1.rds', meth))
```

```{r}
meth <- "Raw"

R <- 1000
print(meth)
preproc.meth <- preproc[[meth]]
X = preproc.meth$graphs; Y = preproc.meth$covariates$Sex; Z = preproc.meth$covariates$Age
DY <- dist(Y); DZ <- dist(Z)
d <- dim(X)[2]
res <- do.call(rbind, mclapply(1:d, function(i) {
  tryCatch({
    print(i)
    DX <- dist(X[,i])
    test <- pdcor.test(DX, DY, DZ, R=R)
    return(data.frame(Edge=i, Statistic=test$statistic,
                      p.value=(1 + sum(test$replicates >= test$statistic))/(1 + R)))
  }, error=function(e) {NULL})
}, mc.cores = ncores)) %>% mutate(Method=meth)

saveRDS(res, sprintf('../data/dcorr/outputs_edgewise_%s_2.rds', meth))
```


```{r}
meth <- "Raw"

R <- 1000
print(meth)
preproc.meth <- preproc[[meth]]
X = preproc.meth$graphs; Y = preproc.meth$covariates$Sex; Z = preproc.meth$covariates$Age
DY <- dist(Y); DZ <- dist(Z)
d <- dim(X)[2]
res <- do.call(rbind, mclapply(8973:13586, function(i) {
  tryCatch({
    print(i)
    DX <- dist(X[,i])
    test <- pdcor.test(DX, DY, DZ, R=R)
    return(data.frame(Edge=i, Statistic=test$statistic,
                      p.value=(1 + sum(test$replicates >= test$statistic))/(1 + R)))
  }, error=function(e) {NULL})
}, mc.cores = ncores)) %>% mutate(Method=meth)

saveRDS(res, sprintf('../data/dcorr/outputs_edgewise_%s_3.rds', meth))
```

```{r}
meth <- "Raw"

R <- 1000
print(meth)
preproc.meth <- preproc[[meth]]
X = preproc.meth$graphs; Y = preproc.meth$covariates$Sex; Z = preproc.meth$covariates$Age
DY <- dist(Y); DZ <- dist(Z)
d <- dim(X)[2]
res <- do.call(rbind, mclapply(1:d, function(i) {
  print(i)
  DX <- dist(X[,i])
  test <- pdcor.test(DX, DY, DZ, R=R)
  return(data.frame(Edge=i, Statistic=test$statistic,
                    p.value=(1 + sum(test$replicates >= test$statistic))/(1 + R)))
}, mc.cores = ncores)) %>% mutate(Method=meth)

saveRDS(res, sprintf('../data/dcorr/outputs_edgewise_%s.rds', meth))
```


```{r}
meth <- "ComBat"

R <- 1000
print(meth)
preproc.meth <- preproc[[meth]]
X = preproc.meth$graphs; Y = preproc.meth$covariates$Sex; Z = preproc.meth$covariates$Age
DY <- dist(Y); DZ <- dist(Z)
d <- dim(X)[2]
res <- do.call(rbind, mclapply(1:d, function(i) {
  print(i)
  DX <- dist(X[,i])
  test <- pdcor.test(DX, DY, DZ, R=R)
  return(data.frame(Edge=i, Statistic=test$statistic,
                    p.value=(1 + sum(test$replicates >= test$statistic))/(1 + R)))
}, mc.cores = ncores)) %>% mutate(Method=meth)

saveRDS(res, '../data/dcorr/outputs_edgewise_combt.rds')
```


```{r}
meth <- "causal ComBat"

R <- 1000
print(meth)
preproc.meth <- preproc[[meth]]
X = preproc.meth$graphs; Y = preproc.meth$covariates$Sex; Z = preproc.meth$covariates$Age
DY <- dist(Y); DZ <- dist(Z)
d <- dim(X)[2]
res <- do.call(rbind, mclapply(1:d, function(i) {
  print(i)
  DX <- dist(X[,i])
  test <- pdcor.test(DX, DY, DZ, R=R)
  return(data.frame(Edge=i, Statistic=test$statistic,
                    p.value=(1 + sum(test$replicates >= test$statistic))/(1 + R)))
}, mc.cores = ncores)) %>% mutate(Method=meth)

saveRDS(res, '../data/dcorr/outputs_edgewise_causcombt.rds')
```