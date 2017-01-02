
lbm_autoregressive_process = function( yt, lag=1, order=1 ) {

    require(LaplacesDemon)

    T <- length(yt)
    L <- 1 #Autoregressive lags
    mon.names <- "LP"
    parm.names <- as.parm.names(list(alpha=0, phi=0, sigma=0))
    pos.alpha <- grep("alpha", parm.names)
    pos.phi <- grep("phi", parm.names)
    pos.sigma <- grep("sigma", parm.names)
    PGF <- function(Data) {
      alpha <- rnorm(1)
      phi <- runif(1,-1,1)
      sigma <- rhalfcauchy(1,5)
      return(c(alpha, phi, sigma))
    }
    MyData <- list(L=L, PGF=PGF, T=T, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.phi=pos.phi,
    pos.sigma=pos.sigma, yt=yt)

    Model <- function(parm, Data) {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    phi <- parm[Data$pos.phi]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    phi.prior <- sum(dnormv(phi, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    i = (Data$L+1) : Data$T
    j = 1:(Data$T-Data$L)
    mu[i] <- rep(alpha, Data$T) + phi*Data$yt[j]
    LL <- sum( dnorm(Data$yt[-c(1:Data$L)], mu[-c(1:Data$L)], sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }

    Initial.Values <- c(rep(0,2), 1)

    parm=MyData$PGF()
    Data = MyData
    LaplacesDemon()
  }
