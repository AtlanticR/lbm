
library(spBayes)

data(meuse)
coordinates(meuse) = ~x+y

##Collect samples
y = log(meuse$zinc)
x = meuse$dist
fit <- glm( y~x-1, family="gaussian")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

coords = coordinates(meuse)

m.1 <- spLM( y~1, coords=coords, knots=c(6,6),
             starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n, size=weights, prob=p)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)


require(MBA)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(coords, label=paste("(",y,")",sep=""))

surf <- mba.surf(cbind(coords,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(coords, label=paste("(",weights,")",sep=""))

----
# dynamic space-time model

 ## Not run:
require(spBayes)

     data("NETemp.dat")
     ne.temp <- NETemp.dat
     
     set.seed(1)
     
     ##take a chunk of New England
     ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3000000,]
     
     ##subset first 2 years (Jan 2000 - Dec. 2002)
     y.t <- ne.temp[,4:27]
     N.t <- ncol(y.t) ##number of months
     n <- nrow(y.t) ##number of observation per months
     
     ##add some missing observations to illistrate prediction
     miss <- sample(1:N.t, 10)
     holdout.station.id <- 5
     y.t.holdout <- y.t[holdout.station.id, miss]
     y.t[holdout.station.id, miss] <- NA
     
     ##scale to km
     coords <- as.matrix(ne.temp[,c("UTMX", "UTMY")]/1000)
     max.d <- max(iDist(coords))
     
     ##set starting and priors
     p <- 2 #number of regression parameters in each month
     
     starting <- list("beta"=rep(0,N.t*p), "phi"=rep(3/(0.5*max.d), N.t),
                      "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                      "sigma.eta"=diag(rep(0.01, p)))
     
     tuning <- list("phi"=rep(5, N.t)) 
     
     priors <- list("beta.0.Norm"=list(rep(0,p), diag(1000,p)),
                    "phi.Unif"=list(rep(3/(0.9*max.d), N.t), rep(3/(0.05*max.d), N.t)),
                    "sigma.sq.IG"=list(rep(2,N.t), rep(10,N.t)),
                    "tau.sq.IG"=list(rep(2,N.t), rep(5,N.t)),
                    "sigma.eta.IW"=list(2, diag(0.001,p)))
     
     ##make symbolic model formula statement for each month
     mods <- lapply(paste(colnames(y.t),'elev',sep='~'), as.formula)
     
     n.samples <- 2000
     
     m.1 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), coords=coords,
                    starting=starting, tuning=tuning, priors=priors, get.fitted =TRUE,
                    cov.model="exponential", n.samples=n.samples, n.report=25) 
     
     burn.in <- floor(0.75*n.samples)


     quant <- function(x){quantile(x, prob=c(0.5, 0.025, 0.975))}
     
     beta <- apply(m.1$p.beta.samples[burn.in:n.samples,], 2, quant)
     beta.0 <- beta[,grep("Intercept", colnames(beta))]
     beta.1 <- beta[,grep("elev", colnames(beta))]
     
     plot(m.1$p.beta.0.samples)
     
     par(mfrow=c(2,1))
     plot(1:N.t, beta.0[1,], pch=19, cex=0.5, xlab="months", ylab="beta.0", ylim=range(beta.0))
     arrows(1:N.t, beta.0[1,], 1:N.t, beta.0[3,], length=0.02, angle=90)
     arrows(1:N.t, beta.0[1,], 1:N.t, beta.0[2,], length=0.02, angle=90)
     
     plot(1:N.t, beta.1[1,], pch=19, cex=0.5, xlab="months", ylab="beta.1", ylim=range(beta.1))
     arrows(1:N.t, beta.1[1,], 1:N.t, beta.1[3,], length=0.02, angle=90)
     arrows(1:N.t, beta.1[1,], 1:N.t, beta.1[2,], length=0.02, angle=90)
     
     theta <- apply(m.1$p.theta.samples[burn.in:n.samples,], 2, quant)
     sigma.sq <- theta[,grep("sigma.sq", colnames(theta))]
     tau.sq <- theta[,grep("tau.sq", colnames(theta))]
     phi <- theta[,grep("phi", colnames(theta))]
     
     par(mfrow=c(3,1))
     plot(1:N.t, sigma.sq[1,], pch=19, cex=0.5, xlab="months", ylab="sigma.sq", ylim=range(sigma.sq))
     arrows(1:N.t, sigma.sq[1,], 1:N.t, sigma.sq[3,], length=0.02, angle=90)
     arrows(1:N.t, sigma.sq[1,], 1:N.t, sigma.sq[2,], length=0.02, angle=90)
     
     plot(1:N.t, tau.sq[1,], pch=19, cex=0.5, xlab="months", ylab="tau.sq", ylim=range(tau.sq))
     arrows(1:N.t, tau.sq[1,], 1:N.t, tau.sq[3,], length=0.02, angle=90)
     arrows(1:N.t, tau.sq[1,], 1:N.t, tau.sq[2,], length=0.02, angle=90)
     
     plot(1:N.t, 3/phi[1,], pch=19, cex=0.5, xlab="months", ylab="eff. range (km)", ylim=range(3/phi))
     arrows(1:N.t, 3/phi[1,], 1:N.t, 3/phi[3,], length=0.02, angle=90)
     arrows(1:N.t, 3/phi[1,], 1:N.t, 3/phi[2,], length=0.02, angle=90)
     
     y.hat <- apply(m.1$p.y.samples[,burn.in:n.samples], 1, quant)
     y.hat.med <- matrix(y.hat[1,], ncol=N.t)
     y.hat.up <- matrix(y.hat[3,], ncol=N.t)
     y.hat.low <- matrix(y.hat[2,], ncol=N.t)
     
     y.obs <- as.vector(as.matrix(y.t[-holdout.station.id, -miss]))
     y.obs.hat.med <- as.vector(y.hat.med[-holdout.station.id, -miss])
     y.obs.hat.up <- as.vector(y.hat.up[-holdout.station.id, -miss])
     y.obs.hat.low <- as.vector(y.hat.low[-holdout.station.id, -miss])
     
     y.ho <- as.matrix(y.t.holdout)
     y.ho.hat.med <- as.vector(y.hat.med[holdout.station.id, miss])
     y.ho.hat.up <- as.vector(y.hat.up[holdout.station.id, miss])
     y.ho.hat.low <- as.vector(y.hat.low[holdout.station.id, miss])
     
     par(mfrow=c(2,1))
     plot(y.obs, y.obs.hat.med, pch=19, cex=0.5, xlab="observed",

     ylab="fitted", main="Observed vs. fitted")
     arrows(y.obs, y.obs.hat.med, y.obs, y.obs.hat.up, length=0.02, angle=90)
     arrows(y.obs, y.obs.hat.med, y.obs, y.obs.hat.low, length=0.02, angle=90)
     lines(-50:50, -50:50, col="blue")
     
     plot(y.ho, y.ho.hat.med, pch=19, cex=0.5, xlab="observed",
     ylab="predicted", main="Observed vs. predicted")
     arrows(y.ho, y.ho.hat.med, y.ho, y.ho.hat.up, length=0.02, angle=90)
     arrows(y.ho, y.ho.hat.med, y.ho, y.ho.hat.low, length=0.02, angle=90)
     lines(-50:50, -50:50, col="blue")
     ## End(Not run)
    



