## ----fig.cap="Lake of Menteith, Scotland", message=FALSE-----------------
library(bayess)
data("Menteith")
iter <- 800
burn <- iter/4 + 1
n <- prod(dim(Menteith))
k <- 6
image(as.matrix(Menteith),
      asp=1,xaxt='n',yaxt='n',col=gray(0:255/255))

## ------------------------------------------------------------------------
library(bayesImageS)
bcrit <- log(1 + sqrt(k))
beta <- sort(c(seq(0,1,by=0.1),seq(1.05,1.15,by=0.05),bcrit-0.05,bcrit-0.02,bcrit+0.02,
            seq(1.3,1.4,by=0.05),seq(1.5,2,by=0.1),2.5,3))
mask <- matrix(1, nrow=sqrt(n), ncol=sqrt(n))
neigh <- getNeighbors(mask, c(2,2,0,0))
block <- getBlocks(mask, 2)
edges <- getEdges(mask, c(2,2,0,0))
maxS <- nrow(edges)
E0 <- maxS/k
V0 <- maxS*(1/k)*(1 - 1/k)

## ----eval=FALSE----------------------------------------------------------
#  library(doParallel)
#  cores <- min(detectCores(), length(beta))
#  print(paste("Parallel computation using",cores,"CPU cores:",
#              iter,"iterations for",length(beta),"values of beta."))
#  cl <- makeCluster(cores)
#  print(cl)
#  clusterSetRNGStream(cl)
#  registerDoParallel(cl)

## ----eval=FALSE----------------------------------------------------------
#  tm <- system.time(matu <- foreach(i=1:length(beta),
#            .packages=c("bayesImageS"), .combine='cbind') %dopar% {
#    res <- swNoData(beta[i],k,neigh,block,iter)
#    res$sum
#  })
#  stopCluster(cl)

## ----echo=FALSE----------------------------------------------------------
data("simSW", package = "bayesImageS")
print(simSW$tm)

## ----fig.cap="Piecewise linear approximation"----------------------------
plot(x=beta, y=colMeans(simSW$matu),main="",
     xlab=expression(beta),ylab=expression(S(z)))
lrcst=approxfun(beta,colMeans(simSW$matu))
curve(lrcst,0,max(beta),add=T,col="blue")
abline(v=bcrit,col="red",lty=3)
abline(h=maxS,col=2,lty=2)
points(0,E0,col=2,pch=2)

## ----eval=FALSE----------------------------------------------------------
#  library(rstan)
#  options(mc.cores = min(4, parallel::detectCores()))
#  dat <- list(M=length(beta), N=iter-burn+1, maxY=maxS, e0=E0, v0=V0,
#              Vlim=2*maxS*log(maxS)/pi, tcrit=bcrit, y=t(simSW$matu[burn:iter,]), t=beta)
#  tm2 <- system.time(fit <- sampling(PFAB, data = dat, verbose=TRUE, iter=5000,
#                              control = list(adapt_delta = 0.9, max_treedepth=20)))

## ----echo=FALSE----------------------------------------------------------
data("fitStan", package = "bayesImageS")
print(fitStan$tm)

## ----results='asis'------------------------------------------------------
library(knitr)
kable(fitStan$summary, caption="95% highest posterior density intervals for pseudo-Voigt peaks",
      align = 'lrrrrrr', digits=3)

## ------------------------------------------------------------------------
ft <- function(t, tC, e0, ecrit, v0, vmax1, vmax2, phi1, phi2) {
  sqrtBcritPhi = sqrt(tC)*phi1
  fval <- numeric(length(t))
  for (i in 1:length(t)) {
    if (t[i] <= tC) {
      sqrtBdiffPhi = sqrt(tC - t[i])*phi1
      fval[i] <- e0 + t[i]*v0 - ((2*(vmax1-v0))/(phi1^2))*((sqrtBcritPhi + 1)/exp(sqrtBcritPhi) - (sqrtBdiffPhi + 1)/exp(sqrtBdiffPhi));
    } else {
      sqrtBdiff = sqrt(t[i] - tC)
      fval[i] <- ecrit - ((2*vmax2)/phi2)*(sqrtBdiff/exp(phi2*sqrtBdiff) + (exp(-phi2*sqrtBdiff) - 1)/phi2);
    }
  }
  return(fval)
}

plot(range(beta),range(simSW$matu),type='n',xlab=expression(beta),ylab=expression(S(z)))
idx <- burn+sample.int(iter-burn+1,size=20)
abline(v=bcrit,col="red",lty=3)
abline(h=maxS,col=2,lty=2)
points(rep(beta,each=20),simSW$matu[idx,],pch=20)
lines(beta, ft(beta, bcrit, E0, 14237, V0, 59019, 124668, 4.556, 6.691),
      col=4, lwd=2)

## ----error=FALSE, warning=FALSE------------------------------------------
residMx <- matrix(nrow=iter-burn+1, ncol=length(beta))
for (b in 1:length(beta)) {
  residMx[,b] <- simSW$matu[burn:iter,b] - ft(beta[b], bcrit, E0, 14237, V0, 59019, 124668, 4.556, 6.691)
}

dfdt <- function(t, tC, V0, Vmax1, Vmax2, r1, r2) {
  ifelse(t < tC,
         V0 + (Vmax1-V0)*exp(-r1*sqrt(tC - t)),
         Vmax2*exp(-r2*sqrt(t - tC)))
}

plot(range(beta),range(residMx),type='n',xlab=expression(beta),ylab="residuals")
abline(h=0,lty=2,col=4,lwd=2)
points(rep(beta,each=iter-burn+1),residMx,pch='.',cex=3)

x <- sort(c(seq(0,3,by=0.01),bcrit))
lines(x, 3*sqrt(dfdt(x, bcrit, V0, 59019, 124668, 4.556, 6.691)), col=2, lwd=2)
lines(x, -3*sqrt(dfdt(x, bcrit, V0, 59019, 124668, 4.556, 6.691)), col=2, lwd=2)

## ------------------------------------------------------------------------
mh <- list(algorithm="aux", bandwidth=0.02, Vmax1=59019,
           Vmax2=124668, E0=E0, Ecrit=14237, phi1=4.556,
           phi2=6.691, factor=1, bcrit=bcrit, V0=V0)
priors <- list()
priors$k <- k
priors$mu <- c(0, 50, 100, 150, 200, 250)
priors$mu.sd <- rep(10,k)
priors$sigma <- rep(20,k)
priors$sigma.nu <- rep(5, k)
priors$beta <- c(0,3)
iter <- 1e4
burn <- iter/2
y <- as.vector(as.matrix(Menteith))
print(date())
tm3 <- system.time(resPFAB <- 
             mcmcPotts(y,neigh,block,priors,mh,iter,burn))
print(tm3)

## ----eval=FALSE----------------------------------------------------------
#    mh <- list(algorithm="ex", bandwidth=0.02, auxiliary=200)
#    tm4 <- system.time(resAEA <-
#                         mcmcPotts(y,neigh,block,priors,mh,iter,burn))
#    resAEA$time <- tm4[3]

## ----echo=FALSE----------------------------------------------------------
data("aeaMenteith", package = "bayesImageS")
print(resAEA$time)

## ------------------------------------------------------------------------
densPFAB <- density(resPFAB$beta[burn:iter])
densAEA <- density(resAEA$beta[burn:iter])
plot(densAEA, col=4, lty=2, lwd=2, main="", xlab=expression(beta),
     xlim=range(resPFAB$beta[burn:iter],resAEA$beta[burn:iter]))
lines(densPFAB, col=2, lty=3, lwd=3)
abline(h=0,lty=2)
legend("topright",legend=c("AEA","PFAB"),col=c(4,2),lty=c(2,3),
       lwd=3)

