library(rstan)
library(earth)
library(ggplot2)

phimat <- c("mean_obs_ons" = expression(paste("Male population size: ",mu[pop])),
            "rho[1]" = expression(paste("Prob. man is GMSM: ",rho[G])),
            "rho[2]" = expression(paste("Prob. man is NGMSM: ", rho[N])),
            "aux_handd"      = expression(paste("P(GMSM diag in GUM clinic): ",a[H])),
            "aux_sophid"     = expression(paste("SOPHID bias: ", a[S])),
            "aux_delta[1]"   = expression(paste("Diagnosed prop., GMSM: ", a[delta*G])),
            "gamma[1]"  = expression(paste("P(prev undiag): ", gamma[1])),
            "gamma[2]"  = expression(paste("P(test offered): ", gamma[2])),
            "gamma[3]"  = expression(paste("P(test accepted): ", gamma[3])),
            "gamma[4]"  = expression(paste("P(new diagnosis): ", gamma[4])),
            "aux_gumcad[1]"  = expression(paste("Excess prev for unoffered: ",a[UN])),
            "aux_gumcad[2]"  = expression(paste("Excess prev if opt out: ", a[OP])),
            "or_gmshs"       = expression(paste("OR between NGMSM/GMSM: ", or^{(GM)})),
            "aux_delta[2]"   = expression(paste("Diagnosed prop., NGMSM: ",a[delta*N])),
            "piga"  = expression(paste("Prevalence from GUM Anon: ", pi^{(GA)}))
            )
phi <- names(phimat)

phimat.nogu <- c("mean_obs_ons" = expression(paste("Male population size: ",mu[pop])),
            "rho[1]" = expression(paste("Prob. man is GMSM: ",rho[G])),
            "rho[2]" = expression(paste("Prob. man is NGMSM: ", rho[N])),
            "aux_handd"      = expression(paste("P(GMSM diag in GUM clinic): ",a[H])),
            "aux_sophid"     = expression(paste("SOPHID bias: ", a[S])),
            "aux_delta[1]"   = expression(paste("Diagnosed prop., GMSM: ", a[delta*G])),
            "pigd"  = expression(paste("Prev of newly-diag: ", pi^{(GD)})),
            "piun"  = expression(paste("Prev among unoffered: ", pi^{(UN)})),
            "piop"  = expression(paste("Prev among refusers: ", pi^{(OP)})),
            "or_gmshs"       = expression(paste("OR between NGMSM/GMSM: ", or^{(GM)})),
            "aux_delta[2]"   = expression(paste("Diagnosed prop., NGMSM: ",a[delta*N])),
            "piga"  = expression(paste("Prevalence from GUM Anon: ", pi^{(GA)}))
            )
phi.nogu <- names(phimat.nogu)

alphamat <- c("pidelta[1]"= expression((pi*delta)[G]),
              "pidelta[2]"= expression((pi*delta)[N]),
              "pinodelta[1]" = expression(bar((pi*delta))[G]),
              "pinodelta[2]"= expression(bar((pi*delta))[N]),
              "mudelta[1]"= expression(mu[DG]),
              "mudelta[2]"= expression(mu[DN]),
              "munodelta[1]"= expression(mu[UG]),
              "munodelta[2]"= expression(mu[UN]),
              "mu[4]"=expression(mu))
alpha <- names(alphamat)

evppi <- function(sam, phi) { 
    pe <- matrix(nrow=length(phi), ncol=length(alpha))
    rownames(pe) <- phi
    colnames(pe) <- alpha
    for (j in seq_along(alpha)) {
        Y <- sam[,alpha[j]]
        for (i in seq_along(phi)){
            X <- sam[,phi[i]]
            mod <- earth(X, Y)
            pe[i,j] <- var(fitted(mod)) / var(Y) ## EVPPI as prop of EVPI
        }
    }
    pe
}

evppi.plot <- function(pe, phi, phimat) { 
    pt <- data.frame(value=as.vector(pe),
                     input=factor(rep(rownames(pe), ncol(pe)), levels=rev(phi)),
                     output=factor(rep(colnames(pe), each=nrow(pe)), levels=names(alphamat)))
    ggplot(pt, aes(output, input)) +
      geom_raster(aes(fill=value)) +
      theme(axis.text.x = element_text(angle=0,vjust=1)) + # was hjust=0
      theme(strip.text.x = element_text(size = 7, angle = 90)) +
      xlab("Prevalences               Numbers of people") +
      ylab(expression(paste("Uncertain input parameter ",phi[r]))) +
      scale_fill_continuous(guide = "colorbar", name=expression(paste("EVPPI(", phi[r], ") / var(",alpha[s],")"))) +
      scale_y_discrete(labels=phimat) +
      scale_x_discrete(labels=alphamat)
}

pe <- evppi(sam, phi)
penogu <- evppi(samnogu, phi.nogu)
pegudnd <- evppi(samgudnd, phi)
save(pe, penogu, pegudnd, file="evppi.rda")
load(file="evppi.rda")

if (0) { 
pdf("../../write/evi/evppi.pdf")
evppi.plot(pe, phi, phimat) 
dev.off()

pdf("../../write/evi/evppi-nogu.pdf")
evppi.plot(penogu, phi.nogu, phimat.nogu)
dev.off()

pdf("../../write/evi/evppi-gudnd.pdf")
evppi.plot(pegudnd, phi, phimat)
dev.off()
}

## standard errors need an extra cross validation step: takes time.

evppi.se <- function(sam, phi) { 
    pe <- matrix(nrow=length(phi), ncol=length(alpha))
    rownames(pe) <- phi
    colnames(pe) <- alpha
    for (j in seq_along(alpha)) {
        Y <- sam[,alpha[j]]
        for (i in seq_along(phi)){
            X <- sam[,phi[i]]
            mod <- earth(X, Y, nfold=10, ncross=30, varmod.method="const", Get.leverages=TRUE)
            se <- sqrt(mod$varmod$model.var)
            B <- 1000
            evppi.rep <- numeric(B)
            for(i in 1:B)
                evppi.rep[i] <- var(rnorm(length(fitted(mod)), fitted(mod), se)) 
            pe[i,j] <- sd(evppi.rep) / var(Y)
        }
    }
    pe
}

if (0) {
    pese <- evppi.se(sam, phi) # negligible 
}
