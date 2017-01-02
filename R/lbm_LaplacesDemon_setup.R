
lbm_LaplacesDemon_setup = function(DS="example.data", Data=NULL) {

  if (DS=="spatial.test") {  

    require(sp)
    data(meuse)
    meuse =meuse[ sample.int(nrow(meuse), 50),]
    nKs = nrow( meuse  )  # knots 
    xrange = range (c(meuse$x))
    yrange = range (c(meuse$y))
    dx = diff(xrange)
    dy = diff(yrange)
    dd = max( dx, dy )
    coordsK = data.frame( plon=(meuse$y-yrange[1])/dd, plat=(meuse$x-xrange[1])/dd) 
    y = log(meuse$zinc)  # LD likes to have a "y" as dep variable
    dKK = as.matrix(dist( coordsK, diag=TRUE, upper=TRUE)) # distance matrix between knots
    Data = list(
      n = nKs,  # required for LaplacesDemon
      nKs=nKs,
      dKK=dKK,
      y=y  
    )

    Data$mon.names = c( "LP" )
    Data$parm.names = as.parm.names(list(muKs=rep(0,Data$nKs), sigma=rep(0,2), phi=0 ))
    Data$pos = list(
      muKs = grep("muKs", Data$parm.names),
      sigma = grep("sigma", Data$parm.names),
      phi = grep("phi", Data$parm.names)
    )
    Data$PGF = function(Data) {
      sigma = runif(2,0.1,10)
      phi = runif(1,1,5)
      muKs = mvnfast::rmvn(1, rep(0,Data$nKs), sigma[2]*sigma[2]*exp(-phi*Data$dKK ) )
      return(c(muKs, sigma, phi))
    }
    Data$PGF  = compiler::cmpfun(Data$PGF)

    Data$Model = function(parm, Data){
      muKs = parm[Data$pos$muKs]
      # parm[Data$pos$kappa] = kappa = LaplacesDemonCpp::interval(parm[Data$pos$kappa], 1e-9, Inf)
      # parm[Data$pos$kappa] = kappa = 1
      parm[Data$pos$sigma] = sigma = LaplacesDemonCpp::interval(parm[Data$pos$sigma], 1e-9, Inf)
      parm[Data$pos$phi] = phi = LaplacesDemonCpp::interval(parm[Data$pos$phi], 1, 5)
      # rhoKs = exp(-phi * Data$dKK)^kappa   ## spatial correlation
      rhoKs = exp(-phi * Data$dKK )  ## spatial correlation
      
      covKs = sigma[2]*sigma[2] * rhoKs
      muKs.prior =  try(mvnfast::dmvn( muKs, rep(0, Data$nKs), sigma=covKs, log=TRUE  ))
      sigma.prior = sum(dgamma(sigma, 1, 100, log=TRUE))
      phi.prior = dunif(phi, 1, 5, log=TRUE)
      # kappa.prior = dgamma(kappa, 1, 100 log=TRUE)

      ### Interpolation
      errorSpatialK = rowSums(rhoKs / rowSums(rhoKs) * matrix(muKs, Data$nKs, Data$nKs, byrow=TRUE) ) # cov/cov --> sigmasq cancels leaving rho
      
      muK = errorSpatialK
      yK = rnorm(Data$nKs)*sigma[1] + muK 
      LL = sum(dnorm(Data$y, muK, sigma[1], log=TRUE))
      
      ### Log-Posterior
      LP = LL + muKs.prior + sigma.prior + phi.prior # + kappa.prior
      Modelout = list(LP=LP, Dev=-2*LL, Monitor=c(LP), yhat=yK, parm=parm)
      return(Modelout)
    }

    Data$Model.ML  = compiler::cmpfun( function(...) (Data$Model(...)$Dev / 2) )  # i.e. - log likelihood
    Data$Model.PML = compiler::cmpfun( function(...) (- Data$Model(...)$LP) ) #i.e., - log posterior 
    Data$Model = compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed

    print (Data$Model( parm=Data$PGF(Data), Data ) ) # test to see if return values are sensible

    return(Data)
  }


  # ----------------------------------------------
  
  if (DS=="spatial.predictive.process") {  
  
## based upon methods of: Viana, M., Jackson, A. L., Graham, N. and Parnell, A. C. 2012. 
## Disentangling spatio-temporal processes in a hierarchical system: a case study in fisheries discards. –
## Ecography 35:
## and Banerjee et al 2008 J. R. Statist. Soc. B (2008) 70, Part 4, pp. 825–848 Gaussian predictive process models for large spatial data sets ( sudipto.bol.ucla.edu/ResearchPapers/BGFS.pdf )
    require(sp)
    data(meuse)
    data(meuse.grid)
    nKs = nrow( meuse  )  # knots 
    nPs = nrow( meuse.grid )
    xrange = range (c(meuse$x, meuse.grid$x))
    yrange = range (c(meuse$y, meuse.grid$y))
    dx = diff(xrange)
    dy = diff(yrange)
    dd = max( dx, dy )
    coordsK = data.frame( plon=(meuse$y-yrange[1])/dd, plat=(meuse$x-xrange[1])/dd) 
    coordsP = data.frame( plon=(meuse.grid$y-yrange[1])/dd, plat=(meuse.grid$x-xrange[1])/dd) 
    # design matrix
    XK = as.matrix( data.frame( intercept=rep(1, nKs), dist=meuse$dist )) 
    XP = as.matrix( data.frame( intercept=rep(1, nPs), dist=meuse.grid$dist ))
    y = log(meuse$zinc)  # LD likes to have a "y" as dep variable
    dKK = as.matrix(dist( coordsK, diag=TRUE, upper=TRUE)) # distance matrix between knots
    dPK = matrix(0, nPs, nKs ) # distance from knot to prediction locations
    for (i in 1:nPs) {
      dPK[i,] = sqrt(( coordsK$plon - coordsP$plon[i] )^2 + (coordsK$plat - coordsP$plat[i] )^2)
    }
    Data = list(
      n = nKs,  # required for LaplacesDemon
      nKs=nKs,
      nPs=nPs,
      dKK=dKK,
      dPK=dPK,
      XP=XP,
      XK=XK,
      y=y  
    )
    Data$mon.names = c( "LP" )
    Data$parm.names = as.parm.names(list(muKs=rep(0,Data$nKs), beta=rep(0,2), sigma=rep(0,2), phi=0))
    Data$pos = list(
      muKs = grep("muKs", Data$parm.names),
      beta = grep("beta", Data$parm.names),
      sigma = grep("sigma", Data$parm.names),
      phi = grep("phi", Data$parm.names)
    )
    Data$PGF = function(Data) {
      beta = rnorm(2)
      sigma = runif(2,0.1,10)
      phi = runif(1,1,5)
      kappa = 1
      muKs = LaplacesDemonCpp::rmvn(1, rep(0,Data$nKs), sigma[2]*sigma[2]*exp(-phi*Data$dKK)^kappa )
      return(c(muKs, beta, sigma, phi))
    }
    Data$Model = function(parm, Data){
      beta = parm[Data$pos$beta]
      muKs = parm[Data$pos$muKs]
      kappa = 1
      parm[Data$pos$sigma] = sigma = LaplacesDemonCpp::interval(parm[Data$pos$sigma], 1, Inf)
      parm[Data$pos$phi] = phi = LaplacesDemonCpp::interval(parm[Data$pos$phi], 1, 5)
      
      # predictive process:
      rhoKs = exp(-phi * Data$dKK)^kappa   ## spatial correlation
      rhoPs = exp(-phi * Data$dPK)^kappa   ## spatial correlation
      covKs = sigma[2]*sigma[2] * rhoKs
      muKs.prior = mvtnorm::dmvnorm( muKs, rep(0, Data$nKs), sigma=covKs, log=TRUE  )
      beta.prior = sum(dnormv(beta, 0, 1000, log=TRUE))
      sigma.prior = sum(dhalfcauchy(sigma, 25, log=TRUE))
      phi.prior = dunif(phi, 1, 5, log=TRUE)
      
      ### Interpolation
      muK.beta = tcrossprod(Data$XK, t(beta))
      muP.beta = tcrossprod(Data$XP, t(beta))
      errorSpatialK = rowSums(rhoKs / rowSums(rhoKs) * matrix(muKs, Data$nKs, Data$nKs, byrow=TRUE) ) # cov/cov --> sigmasq cancels leaving rho
      errorSpatialP = rowSums(rhoPs / rowSums(rhoPs) * matrix(muKs, Data$nPs, Data$nKs, byrow=TRUE) )

      muK = muK.beta + errorSpatialK
      muP = muP.beta + errorSpatialP
      # y = rnorm(Data$nKs)*sigma[1] + muK 
      yP = rnorm(Data$nPs)*sigma[1] + muP 
      LL = sum(dnorm(Data$y, muK, sigma[1], log=TRUE))
      # LLP = sum(dnorm( yP, muP, sigma[1], log=TRUE)) # to include predictions into likelihood
      
    ## also possible: 
    #  if (0) {
    #      method1 = function() {
    #          ## This may be mathematically more correct but 10 X slower
    #          covPKs = sigma[2]*sigma[2] * rhoPs
    #          invcovKs = chol2inv(chol(covKs))
    #          errorSpatialP1 = covPKs %*% (invcovKs %*% errorSpatialK) # the nesting speeds up; 
    #          # also use: crossprod(X) instead of X'X ; also crossprod(X,Y) is faster than t(X) %*% Y , etc 
    #      }
    #      method2 = function() {
    #          errorSpatialP0 = rowSums(rhoPs / rowSums(rhoPs) * matrix(muKs, Data$nPs, Data$nKs, byrow=TRUE)  )
    #      }
    #      microbenchmark::microbenchmark( 
    #        method1(),
    #        method2(),
    #        times=10000 
    #      )
    #      expr     min       lq      mean   median       uq      max neval cld
    # method1() 400.130 474.5315 579.54565 492.8235 546.5975 93007.11 10000   b
    # method2()  30.673  33.5750  53.05837  39.9975  46.2095 87168.62 10000  a 
    # 10 times slower
    #  }


      ### Log-Posterior
      LP = LL + beta.prior + muKs.prior + sigma.prior + phi.prior
      Modelout = list(LP=LP, Dev=-2*LL, Monitor=c(LP), yhat=yP, parm=parm)
      return(Modelout)
    }

    Data$Model.ML  = compiler::cmpfun( function(...) (Data$Model(...)$Dev / 2) )  # i.e. - log likelihood
    Data$Model.PML = compiler::cmpfun( function(...) (- Data$Model(...)$LP) ) #i.e., - log posterior 
    Data$Model = compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed

    return(Data)
  }

}
