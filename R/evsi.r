load(file="~/scratch/uncertainty/evi/hiv/sam.rda")
load(file="~/scratch/uncertainty/evi/hiv/samnogu.rda")
library(earth)
library(tidyverse)

out <- "munodelta[4]"
yvar <- var(sam[,out])
yvarnogu <- var(samnogu[,out])
ns <- c(10,50,100,500,1000,10000,100000)

evsi.gumanon <- function(n, sam=sam, out="munodelta[4]"){
    Y <- sam[,out,drop=FALSE]
    nsam <- nrow(Y)
    pu <- sam[,"piga"]
    nrep <- rbinom(nsam, n, pu)
    prep <- nrep / n
    var(fitted(earth(prep, Y)))
}

evsi.ga <- numeric(length(ns))
for (i in seq_along(ns))
    evsi.ga[i] <- evsi.gumanon(ns[i], sam=sam)
evppi.ga <- var(fitted(earth(sam[,"piga"], sam[,out])))
(yvar - evsi.ga)/yvar

evsi.ga.nogu <- numeric(length(ns))
for (i in seq_along(ns))
    evsi.ga.nogu[i] <- evsi.gumanon(ns[i], sam=samnogu)
evppi.ga.nogu <- var(fitted(earth(samnogu[,"piga"], samnogu[,out])))
(yvarnogu - evsi.ga.nogu)/yvarnogu

evsi.gmshs <- function(n, sam=sam, out="munodelta[4]"){
    Y <- sam[,out,drop=FALSE]
    nsam <- nrow(Y)
    ping <- sam[,"prob_gmgum"]
    ngrep <- rbinom(nsam, n, ping)
    nnrep <- n - ngrep
    pg <- sam[,"prob_gmshs[1]"]
    pn <- sam[,"prob_gmshs[2]"]
    pgrep <- (rbinom(nsam, ngrep, pg)+0.5) / (ngrep+1) # posterior mean under jeffreys prior beta(0.5,0.5)
    pnrep <- (rbinom(nsam, nnrep, pn)+0.5) / (nnrep+1)
    orrep <- (pgrep/(1-pgrep)) / (pnrep/(1-pnrep))
    var(fitted(earth(orrep, Y)))
}

evsi.gm <- numeric(length(ns))
for (i in seq_along(ns))
    evsi.gm[i] <- evsi.gmshs(ns[i],sam=sam)
evppi.gm <- var(fitted(earth(sam[,"or_gmshs"], sam[,out])))
(yvar - evsi.gm)/yvar 

evsi.gm.nogu <- numeric(length(ns))
for (i in seq_along(ns))
    evsi.gm.nogu[i] <- evsi.gmshs(ns[i], sam=samnogu)
evppi.gm.nogu <- var(fitted(earth(samnogu[,"or_gmshs"], samnogu[,out])))
(yvarnogu - evsi.gm.nogu)/yvarnogu


## used for two graphs in paper, with / without GUMCAD

evsi.plot <- function(evsiga, evsigm, evppiga, evppigm, yvar, ys, labga, labgm, title){
    options(scipen=10000)
    labs <- eval(parse(text=paste0("expression(", paste0(ys,"^2",collapse=","), ")")))
    ggplot(data=NULL, aes(x=ns, y=yvar-evsiga)) +
      geom_point(col="red") +
      geom_line(col="red") +
      geom_point(data=NULL, aes(x=ns, y=yvar-evsigm), col="blue") +
      geom_line(data=NULL, aes(x=ns, y=yvar-evsigm), col="blue") +
      geom_hline(aes(yintercept=yvar), col="black") +
      geom_hline(aes(yintercept=yvar-evppiga), col="red", linetype=2) +
      geom_hline(aes(yintercept=yvar-evppigm), col="blue", linetype=2) +
      geom_text(data=NULL, aes(x=10, y=yvar+(diff(range(ys))/3)^2), label="Original posterior variance", col="blue", hjust=0, vjust=0) +
      geom_text(data=NULL, aes(x=10, y=yvar-evppiga), label="Expected~variance~knowing~pi^{(GA)}", parse=TRUE, col="red", hjust=0, vjust=2) +
      geom_text(data=NULL, aes(x=10, y=yvar-evppigm), label=bquote("Expected~variance~knowing~or^{(GM)}"), parse=TRUE, col="blue", hjust=0, vjust=2) +
      scale_x_continuous(trans="log", breaks=ns) +
      scale_y_continuous(breaks=ys^2, limits=range(ys^2), labels=labs) +
      annotate("text", x=1000, y=labga^2, label="GUM Anon data", col="red", hjust=0) +
      annotate("text", x=1000, y=labgm^2, label="GMSHS data", col="blue", hjust=0) +
      xlab(expression(paste("Planned additional size of survey ", bold(y)))) +
      ylab(expression(paste("Variance of ", mu[U]))) +
      ggtitle(title)
}


## standard errors 
if (0) { 
    
evsi.gumanon.se <- function(n, sam=sam, out="munodelta[4]"){
    Y <- sam[,out,drop=FALSE]
    nsam <- nrow(Y)
    pu <- sam[,"piga"]
    nrep <- rbinom(nsam, n, pu)
    prep <- nrep / n
    mod <- earth(prep, Y, nfold=10, ncross=30, varmod.method="const", Get.leverages=TRUE)
    se <- sqrt(mod$varmod$model.var)
    B <- 1000
    evi.rep <- numeric(B)
    for(i in 1:B)
        evi.rep[i] <- var(rnorm(length(fitted(mod)), fitted(mod), se)) 
    sd(evi.rep)
}

evsi.ga.se <- numeric(length(ns))
for (i in seq_along(ns))
    evsi.ga.se[i] <- evsi.gumanon.se(ns[i], sam=sam)

options(scipen=1e+07)
round(cbind(evsi.ga, evsi.ga.se), 2)
evsi.ga.se / evsi.ga ## SE 1% at most of evsi 
    
evsi.gmshs.se <- function(n, sam=sam, out="munodelta[4]"){
    Y <- sam[,out,drop=FALSE]
    nsam <- nrow(Y)
    ping <- sam[,"prob_gmgum"]
    ngrep <- rbinom(nsam, n, ping)
    nnrep <- n - ngrep
    pg <- sam[,"prob_gmshs[1]"]
    pn <- sam[,"prob_gmshs[2]"]
    pgrep <- (rbinom(nsam, ngrep, pg)+0.5) / (ngrep+1) # posterior mean under jeffreys prior beta(0.5,0.5)
    pnrep <- (rbinom(nsam, nnrep, pn)+0.5) / (nnrep+1)
    orrep <- (pgrep/(1-pgrep)) / (pnrep/(1-pnrep))
    mod <- earth(orrep, Y, nfold=10, ncross=30, varmod.method="const", Get.leverages=TRUE)
    se <- sqrt(mod$varmod$model.var)
    B <- 1000
    evi.rep <- numeric(B)
    for(i in 1:B)
        evi.rep[i] <- var(rnorm(length(fitted(mod)), fitted(mod), se)) 
    sd(evi.rep)
}

evsi.gm.se <- numeric(length(ns))
for (i in seq_along(ns))
    evsi.gm.se[i] <- evsi.gmshs.se(ns[i],sam=sam)
round(cbind(evsi.gm, evsi.gm.se), 2)
evsi.gm.se / evsi.gm # all < 1% of evsi

}
