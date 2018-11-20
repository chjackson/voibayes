load("../data/mpesdata.rda")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Base case model
mpes.stan <- stan(file = "mpes.stan", data = mpesdat, init=mpesin, iter = 100000, chains=3, cores=3, seed=1)
summary(mpes.stan)
a <- rstan::extract(mpes.stan, permuted=FALSE)
sam <- matrix(a, nrow=prod(dim(a)[1:2]), ncol=dim(a)[3])
colnames(sam) <- dimnames(a)[[3]]
save(sam, file="sam.rda")

## Alternative assumption a) Undiagnosed prevalence from GUM Anon only, not GUMCAD
mpes.nogu.stan <- stan(file = "mpes-nogu.stan", data = mpesdat, init=mpesin, iter = 100000, chains=3, cores=3, seed=1)
a <- rstan::extract(mpes.nogu.stan, permuted=FALSE)
samnogu <- matrix(a, nrow=prod(dim(a)[1:2]), ncol=dim(a)[3])
colnames(samnogu) <- dimnames(a)[[3]]
save(samnogu, file="samnogu.rda")

## Alternative assumption b) GUMCAD informs both diagnosed and undiagnosed prevalence
mpes.gudnd.stan <- stan(file = "mpes-gudnd.stan", data = mpesdat, init=mpesin, iter = 100000, chains=3, cores=3, seed=1)
a <- rstan::extract(mpes.gudnd.stan, permuted=FALSE)
samgudnd <- matrix(a, nrow=prod(dim(a)[1:2]), ncol=dim(a)[3])
colnames(samgudnd) <- dimnames(a)[[3]]
save(samgudnd, file="samgudnd.rda")
