## Illustrate use of ENBS following EVSI to determine the optimal sample size
## Extra data from GUM Anon in "excluding GUMCAD data and strong priors" scenario

## Calculate EVSI for range of sample sizes

# load("samnogu.rda")
library(earth)
library(ggplot2)
library(mgcv)

out <- "munodelta[4]"
yvarnogu <- var(samnogu[,out])

evsi.gumanon <- function(n, sam=sam, out="munodelta[4]"){
    Y <- sam[,out,drop=FALSE]
    nsam <- nrow(Y)
    pu <- sam[,"piga"]
    nrep <- rbinom(nsam, n, pu)
    prep <- nrep / n
    var(fitted(earth(prep, Y)))
}

## sample sizes to calculate EVSI etc for 
ns <- seq(1,400,by=1)

recalc.evsi <- FALSE

if (recalc.evsi) { 
    evsi.ga.nogu <- numeric(length(ns))
    system.time({
        for (i in seq_along(ns)){
            print(ns)
            evsi.ga.nogu[i] <- evsi.gumanon(ns[i], sam=samnogu)
        }
    })
    save(evsi.ga.nogu, file="evsi.ga.nogu.rda")
} else 
    load(file="evsi.ga.nogu.rda")


## Function for cost of obtaining data from n people given up-front
## cost "up" and per-person cost "n"

## Values of 0 and 17 used here based on costs given by Louise Logan (in email) who was originally involved with GUM Anon
## Includes costs of collecting and analysing specimens, transport and consumables.
## Overall budget 35000
## 1551 specimens (note we just have 100 or so from MSM in London) 
## 4 per specimen to collect
## 8 per specimen to analyse 
## 8k transport + 500 for consumables, must increase with the number of people. 
## 1551*12 + 8000 + 500  = 27112. plausible underspend from 35000. 
## 27112 / 1551 = 17

survey_cost <- function(n, up=0, pp=17){ 
    up + pp*n
}

## Willingness to pay per unit of benefit
## Per reduction in variance by one unit.
## Of estimated number of people mu_U with undiagnosed HIV, as in Figure 6

## Show the WTP that would make the decision problem interesting
## note that current survey size is Gum-Anon 100, GMSHS 1000.
## Want optimal additional size between 10 and 10000

## Current info
## muU=804,  SD 323  with strong prior (base case) 
## muU=5164, SD 3271 with GUMCAD / strong priors excluded.

## (a) Excluding GUMCAD data and strong priors
## Benefit of more data from GUM Anon 

## Define willing to pay per unit of variance

wtpunit <- 3271^2 - 2771^2 ## variance reduction to reduce SD from e.g. 3000 to 2500
wtpperunit <- 5000 # WTP for this reduction
lam <- wtpperunit / wtpunit # WTP for variance reduction of 1

## expected benefit 
eb <- lam*evsi.ga.nogu
## expected cost
ec <- survey_cost(ns)
## expected net benefit of sampling 
enbs <- eb - ec

## smooth them all as functions of sample size
check.fit <- FALSE
if (check.fit){
    plot(ns, eb, type="l")
    lines(fitted(gam(eb ~ s(ns))), col="blue")
    plot(ns, ec, type="l")
    lines(fitted(gam(ec ~ s(ns))), col="blue")
    plot(ns, enbs, type="l")
    lines(fitted(gam(enbs ~ s(ns))), col="blue")
}
eb <- fitted(gam(eb ~ s(ns)))
ec <- fitted(gam(ec ~ s(ns)))
enbs <- fitted(gam(enbs ~ s(ns)))
dat <- data.frame(ns, eb, ec, enbs)
     
p <- ggplot(data=dat, aes(x=ns))
col2 <- "deepskyblue3"
maxbase <- ns[which.max(enbs)]
maxtwice <- ns[which.max(eb*2-ec)]
maxhalf <- ns[which.max(eb/2-ec)]
brks <- sort(c(setdiff(seq(0, 400, by=100), 300), maxbase, maxtwice))
