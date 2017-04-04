
lbm__spate = function( p, dat, pa, sloc, windowsize.half=NULL ) {
  #\\ SPDE solution via FFT using the spate library
  # require(spate)
  # based upon two-step process

  message( "This is not yet finished/debugged .. ")

  pa = lbm_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=windowsize.half, even=TRUE )
  sdTotal=sd(dat[,p$variable$Y], na.rm=T)
  ws = windowsize.half 
  nsq = 2*ws  # this needs to be an even matrix for spate (even=TRUE drops 1 from right edge)
  adims = c(p$nt, nsq, nsq ) 

  dat_id = array_map( "3->1", 
    coords = round( cbind( 
      ( (dat[,p$variables$TIME ] - p$prediction.ts[1] ) / p$tres) + 1 ),
      ( ws + (dat[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
      ( ws + (dat[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1), 
    dims=adims )
  
  xM = array( NA, dim=adims )
  xM[dat_id] = dat[,p$variables$Y]


  ts_gam = lbm__gam( p, dat, px ) # currently only a GAM is enabled for the TS component

  if (is.null( ts_gam)) return(NULL)
  if (ts_gam$lbm_stats$rsquared < p$lbm_rsquared_threshold ) return(NULL)

  # range checks
  # rY = range( dat[,p$variables$Y], na.rm=TRUE)
  # toosmall = which( ts_gam$predictions$mean < rY[1] )
  # toolarge = which( ts_gam$predictions$mean > rY[2] )
  # if (length(toosmall) > 0) ts_gam$predictions$mean[toosmall] = NA   # permit space modelling to fill this in
  # if (length(toolarge) > 0) ts_gam$predictions$mean[toolarge] = NA   
 
  pxts = ts_gam$predictions
  names(pxts)[which(names(pxts)=="mean")] = p$variables$Y
  names(pxts)[which(names(pxts)=="sd")] = paste(p$variables$Y, "sd", sep=".")
  
  ts_gam = NULL
  gc()

  pxts = pxtsts

  nsq = 2*ws + 2 # this needs to be an even matrix for spate
  adims = c(p$nt, nsq, nsq ) 

  pxts_id = array_map( "3->1", 
    coords = round( cbind( 
      ( ws + (pxts[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
      ( ws + (pxts[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1, 
      round( ( pxts[,p$variables$TIME ] - p$prediction.ts[1] ) / p$tres) + 1 )),
    dims=adims )
  
  xM = array( NA, dim=adims )
  xM[pxts_id] = pxts[,p$variables$Y]

  xM2 = matrix( xM, nrow=p$nt )
  xM2[mm] = xM[pxts_id] 

  g = spate.mcmc( y=xM2, n=nsq ) 

  # plot(g, postProcess=TRUE)

  predict <- spate.predict(y=mm, tPred=(1:p$nt), spateMCMC=g, Nsim=500, BurnIn=10, DataModel="Normal" )
  Pmean <- apply(predict, c(1,2), mean)
  Psd <- apply(predict, c(1,2), sd)

  # plot(mean ~ z , dat)
  ss = lm( dat$mean ~ dat[,p$variables$Y], na.action=na.omit )
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$lbm_rsquared_threshold ) return(NULL)

  lbm_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
 
  if (0) {
    par <- c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,
    gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
    spateSim <- spate.sim(par=par,n=20,T=20,seed=4)
    w <- spateSim$w
    ## This is an example to illustrate the use of the MCMC algorithm. 
    ## In practice, more samples (Nmc) are needed for a sufficiently 
    ## large effective sample size.
    spateMCMC <-spate.mcmc( y=w, x=NULL, 
      SV=c(rho0=0.2,sigma2=0.1,zeta=0.25,rho1=0.2,gamma=1, alpha=0.3,muX=0,muY=0,tau2=0.005),
      RWCov=diag(c(0.005,0.005,0.05,0.005, 0.005,0.001,0.0002,0.0002,0.0002)),
      DimRed=TRUE,NFour=29,
      Nmc=10000,BurnIn=2000,seed=4,NCovEst=500, BurnInCovEst=500,trace=FALSE,Padding=TRUE)
  
      ## Make predictions
    predict <- spate.predict(y=w, tPred=(21:23), 
                             spateMCMC=spateMCMC, Nsim = 100, 
                             BurnIn = 10, DataModel = "Normal",seed=4)
    Pmean <- apply(predict,c(1,2),mean)
    Psd <- apply(predict,c(1,2),sd)

    par(mfrow=c(2,3),mar=c(2,2,2,2))
    zlim=c(min(Pmean),max(Pmean))
    for(i in 1:3){
      image(1:20,1:20,matrix(Pmean[i,],nrow=20),zlim=zlim,
            main=paste("Mean predicted field at t=",i+20,sep=""),
            xlab="",ylab="",col=cols())
    }

    zlim=c(min(Psd),max(Psd))
    for(i in 1:3){
      image(1:20,1:20,matrix(Psd[i,],nrow=20),zlim=zlim,
            main=paste("Sd of predicted field at t=",i+20,sep=""),
            xlab="",ylab="",col=cols())
    }

  }

  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}

