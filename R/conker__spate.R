
conker__spate = function( p, x, pa, ws, sloc, smoothness=0.5 ) {
  #\\ SPDE solution via FFT using the spate library
  # require(spate)

  x_r = range(x[,p$variables$LOCS[1]])
  x_c = range(x[,p$variables$LOCS[2]])

  x_nr = diff(x_r)/p$pres + 1
  x_nc = diff(x_c)/p$pres + 1

  
  
  if ( exists("TIME", p$variables)) {
    xi = which( x[, p$variables$TIME]==p$ts[ti]) 
  } else {
    xi =1:nrow(x) 
  } 
 
  ws = trunc( conker_distance_cur / p$pres )
  nsq = 2*ws +1 

  x_id = cbind( ( ws + (x[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
                ( ws + (x[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1, 
                trunc( ( x[,p$variables$TIME ] - p$ts[1] ) / p$tres) + 1 
  )
  
  dups = which( duplicated(x_id) )
  if (length(dups) > 0 ) {
    ud = sort( unique( x_id[dups] ) ) 
    todrop = NULL 
    for ( uu in ud ) {
      ui = which( x_id == uu )
      todrop = c( todrop, ui[-1] )
      x[ui[1],p$variables$Y] = mean( x[ui,p$variables$Y], na.rm=TRUE )
    }
    x = x[-todrop, ]
    x_id =x_id[-todrop, ]
  }
  
  xM = array( NA, dim=c(nsq, nsq, p$nt) )
  xM[x_id] = x[,p$variables$Y]
  
  xM = xM[-nsq,-nsq,] # this needs to be an even matrix
  w = c(xM)

  ##Initial values for optim. This takes a couple of seconds.
  parI <- c(rho0=0.2, sigma2=0.1, zeta=0.25, rho1=0.01, 
    gamma=1, alpha=0.3, muX=0, muY=0, tau2=0.005)
  logInd=c(1,2,3,4,5,9)
  ##Transform to log-scale
  parI[logInd] <- log(parI[logInd])

  ### STOPPED HERE :: missing values in a time slice breaks this method  (when fft is computed by a time sline and there is no data to conduct the FFT .. ) :: revisit once the 2 -phase time-then space based interpolation is completed  .. where this is implements as the "spatial moethod"  


  ##Fourier transform needs to be done only once
  wFT <- real.fft.TS(w, n=(nsq-1), T=p$nt)
  ##ML estimation using optim, takes a couple of seconds
  ##Load the precomputed object a line below to save time
  spateMLE <- optim( par=parI, loglike, wFT=wFT,
      control=list(trace=TRUE,maxit=1000), method="L-BFGS-B",
      lower=c(-10,-10,-10,-10,-10,0,-0.5,-0.5,-10),
      upper=c(10,10,10,10,10,pi/2,0.5,0.5,10),negative=TRUE,
      logScale=TRUE,hessian=TRUE,n=nsq, T=p$nt)
  mle <- spateMLE$par
  mle[logInd] <- exp(mle[logInd])
  sd=sqrt(diag(solve(spateMLE$hessian)))


  x_plons = seq( x_r[1], x_r[2], length.out=x_nr )
  x_plats = seq( x_c[1], x_c[2], length.out=x_nc )

  x_locs = expand.grid( x_plons, x_plats ) # final output grid
  attr( x_locs , "out.attrs") = NULL
  names( x_locs ) = p$variables$LOCS

  x$mean = NA

  pa$mean = NA
  pa$sd = NA

  # locations of the new (output) coord system .. smaller than the data range of x
  pa_r = range(pa[,p$variables$LOCS[1]])
  pa_c = range(pa[,p$variables$LOCS[2]])
  
  pa_nr = diff(pa_r)/p$pres + 1
  pa_nc = diff(pa_c)/p$pres + 1

  pa_plons = seq( pa_r[1], pa_r[2], length.out=pa_nr )
  pa_plats = seq( pa_c[1], pa_c[2], length.out=pa_nc )

  pa_locs = expand.grid( pa_plons, pa_plats ) # final output grid
  attr( pa_locs , "out.attrs") = NULL
  names( pa_locs ) = p$variables$LOCS
  rm( pa_r, pa_c, pa_nr, pa_nc, pa_plons, pa_plats) 


  for ( ti in 1:p$nt ) {
     
    if ( exists("TIME", p$variables)) {
      xi = which( x[, p$variables$TIME]==p$ts[ti]) 
    } else {
      xi =1:nrow(x) 
    } 
    
    # map of row, col indices of input data in the new (output) coordinate system
    x_id = cbind( (x[xi,p$variables$LOCS[1]]-x_r[1])/p$pres + 1, 
                  (x[xi,p$variables$LOCS[2]]-x_c[1])/p$pres + 1 )

    # matrix representation of the output surface
    M = matrix( NA, nrow=x_nr, ncol=x_nc) 
    M[x_id] = x[xi,p$variables$Y] # fill with data in correct locations
    Z = try( fields::image.smooth( M, dx=p$pres, dy=p$pres, theta=p$conker_theta)$z )
  
    if (0) {
      # more control of covariance function .. but not behaving very well and slow .. better to copy internal and strip it down .. TODO
      Z = try( smooth.2d( Y=x[xi,p$variables$Y], x=x[xi,p$variables$LOCS], ncol=x_nc, nrow=x_nr, theta=p$conker_theta,
        cov.function=stationary.cov, Covariance="Exponential", p=smoothness, smoothness=smoothness ) )
      iZ = which( !is.finite( Z$z))
      if (length(iZ) > 0) Z$z[iZ] = NA
      rY = range( x[xi,p$variables$Y], na.rm=TRUE)
      nZ = which( Z$z < rY[1] )
      if (length(nZ) > 0) Z$z[nZ] = NA
      mZ = which( Z$z > rY[2] )
      if (length(mZ) > 0) Z$z[mZ] = NA
      
      x11(); image.plot(Z)
      Z = Z$z
    }
  
    if ( "try-error" %in% class(Z) ) next()
    # match prediction to input data 
    x$mean[xi] = Z[x_id]
    ss = lm( x$mean[xi] ~ x[xi,p$variables$Y], na.action=na.omit)
    if ( "try-error" %in% class( ss ) ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$conker_rsquared_threshold ) next()
    pa_i = ifelse( exists("TIME", p$variables), {which( pa[, p$variables$TIME]==p$ts[ti])}, {1:nrow(pa)} ) 
    Z_i = cbind( ( pa[pa_i,p$variables$LOCS[1]]-x_r[1])/p$pres + 1, 
                  (pa[pa_i,p$variables$LOCS[2]]-x_c[1])/p$pres + 1 )

    # make sure predictions exist .. kernel density can stop prediction beyond a given range if the xwidth/ywidth options are not used and/or the kernel distance (theta) is small 
    if ( any( Z_i<1) ) next()  
    if ( any( Z_i[,1] > x_nr) ) next()
    if ( any( Z_i[,2] > x_nc) ) next()
    pa$mean[pa_i] = Z[Z_i]
    pa$sd[pa_i] = 1
  }

  # plot(mean ~ z , x)
  ss = lm( x$mean ~ x[,p$variables$Y], na.action=na.omit )
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$conker_rsquared_threshold ) return(NULL)

  conker_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, conker_stats=conker_stats ) )  
}

