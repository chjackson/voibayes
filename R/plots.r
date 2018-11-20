library(denstrip)

col1 <- "black"; col2="blue"; col3="red"; col4="purple"
colmin <- "gray92"
wd <- 0.3
ci <- c(0.025, 0.975)
base <- 1
dy <- 0.5
dg <- 0.7
base2 <- base + 3*dy + dg
base3 <- base + 6*dy + 2*dg

prev.plot <- function(sam) { 
    
    par(mar=c(3, 0, 0, 0), mgp=c(2,1,0))

    xmax <- 0.55
    plot(0, type="n", xlim=c(-0.15, xmax), ylim=c(1,7), axes=FALSE, xlab="Prevalence", ylab="")
    lim <- par("usr")
    rect(0, lim[3], lim[2], lim[4], col="gray92", border="gray92")
    axis(1, at=seq(0,1,by=0.05))
    abline(v=seq(0,1,by=0.05), col="white")

    denstrip(sam[,"pi[1]"],        at=base3+2*dy,  width=wd, colmax=col1, colmin=colmin, from=0, to=1, ticks=quantile(sam[,"pi[1]"], ci))
    denstrip(sam[,"pidelta[1]"],   at=base3+dy, width=wd, colmax=col2, colmin=colmin, from=0, to=1, ticks=quantile(sam[,"pidelta[1]"], ci))
    denstrip(sam[,"pinodelta[1]"], at=base3, width=wd, colmax=col3, colmin=colmin, from=0, to=1, ticks=quantile(sam[,"pinodelta[1]"], ci))

    denstrip(sam[,"pi[2]"],        at=base2+2*dy,  width=wd, colmax=col1, colmin=colmin, from=0, to=1, ticks=quantile(sam[,"pi[2]"], ci))
    denstrip(sam[,"pidelta[2]"],   at=base2+dy, width=wd, colmax=col2, colmin=colmin, from=0, to=1, ticks=quantile(sam[,"pidelta[2]"], ci))
    denstrip(sam[,"pinodelta[2]"], at=base2, width=wd, colmax=col3, colmin=colmin, from=0, to=1, ticks=quantile(sam[,"pinodelta[2]"], ci))

    denstrip(sam[,"pi[4]"],        at=base+2*dy, width=wd, colmax=col1, colmin=colmin, from=0, to=1, ticks=quantile(sam[,"pi[4]"], ci))
    denstrip(sam[,"pidelta[4]"],   at=base+dy,  width=wd, colmax=col2, colmin=colmin, from=0, to=1, ticks=quantile(sam[,"pidelta[4]"], ci))
    denstrip(sam[,"pinodelta[4]"], at=base,  width=wd, colmax=col3, colmin=colmin, from=0, to=1, ticks=quantile(sam[,"pinodelta[4]"], ci))

    text(-0.125, c(base3, base2, base)+3*dy, c("GMSM", "NGMSM", "All MSM"), pos=4)
    labs <- c("Overall","Diagnosed","Undiagnosed")
    labsn <- c(expression(paste("Overall ",pi[N])), expression(paste("Diagnosed ", (pi*delta)[N])), expression(paste("Undiagnosed ", bar((pi*delta))[N])))
    labsg <- c(expression(paste("Overall ",pi[G])), expression(paste("Diagnosed ", (pi*delta)[G])), expression(paste("Undiagnosed ", bar((pi*delta))[G])))
    text(-0.17, base+c(2*dy,dy,0), labs, pos=4, col=c(col1,col2,col3))
    text(-0.17, base2+c(2*dy,dy,0), labsn, pos=4, col=c(col1,col2,col3))
    text(-0.17, base3+c(2*dy,dy,0), labsg, pos=4, col=c(col1,col2,col3))

}



nums.plot <- function(sam){
    
    par(mar=c(3, 0, 0, 0), mgp=c(2,1,0))

    xmax <- 11000
    plot(0, type="n", xlim=c(-4000, xmax), ylim=c(1,7), axes=FALSE, xlab="Number of people living with HIV/AIDS", ylab="")
    lim <- par("usr")
    rect(0, lim[3], lim[2], lim[4], col="gray92", border="gray92")
    axis(1, at=seq(0,xmax,by=1000))
    abline(v=seq(0,xmax,by=500), col="white")
    denstrip(sam[,"mu[1]"],        at=base3+2*dy,  width=wd, colmax=col1, colmin=colmin, from=0, ticks=quantile(sam[,"mu[1]"], ci))
    denstrip(sam[,"mudelta[1]"],   at=base3+dy, width=wd, colmax=col2, colmin=colmin, from=0, ticks=quantile(sam[,"mudelta[1]"], ci))
    denstrip(sam[,"munodelta[1]"], at=base3, width=wd, colmax=col3, colmin=colmin, from=0, ticks=quantile(sam[,"munodelta[1]"], ci))
    denstrip(sam[,"mu[2]"],        at=base2+2*dy,  width=wd, colmax=col1, colmin=colmin, from=0, ticks=quantile(sam[,"mu[2]"], ci))
    denstrip(sam[,"mudelta[2]"],   at=base2+dy, width=wd, colmax=col2, colmin=colmin, from=0, ticks=quantile(sam[,"mudelta[2]"], ci))
    denstrip(sam[,"munodelta[2]"], at=base2, width=wd, colmax=col3, colmin=colmin, from=0, ticks=quantile(sam[,"munodelta[2]"], ci))
    denstrip(sam[,"mu[4]"],        at=base+2*dy, width=wd, colmax=col1, colmin=colmin, from=0, ticks=quantile(sam[,"mu[4]"], ci))
    denstrip(sam[,"mudelta[4]"],   at=base+dy,  width=wd, colmax=col2, colmin=colmin, from=0, ticks=quantile(sam[,"mudelta[4]"], ci))
    denstrip(sam[,"munodelta[4]"], at=base,  width=wd, colmax=col3, colmin=colmin, from=0, ticks=quantile(sam[,"munodelta[4]"], ci))
    text(-4250, c(base3, base2, base)+3*dy, c("GMSM", "NGMSM", "All MSM"), pos=4)

    labs <- c("Overall","Diagnosed","Undiagnosed")
    labs <- c(expression(paste("Overall ",mu)), expression(paste("Diagnosed ")), expression(paste("Undiagnosed ")))
    labsn <- c(expression(paste("Overall ",mu[N])), expression(paste("Diagnosed ", mu[DN])), expression(paste("Undiagnosed ", mu[UN])))
    labsg <- c(expression(paste("Overall ",mu[G])), expression(paste("Diagnosed ", mu[DG])), expression(paste("Undiagnosed ", mu[UG])))

    text(-4100, base+c(2*dy,dy,0), labs, pos=4, col=c(col1,col2,col3))
    text(-4100, base2+c(2*dy,dy,0), labsn, pos=4, col=c(col1,col2,col3))
    text(-4100, base3+c(2*dy,dy,0), labsg, pos=4, col=c(col1,col2,col3))

}
