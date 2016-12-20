
hivemod__kerneldensity = function( p, x, pa, nu=NULL, phi=NULL ) {

  #\\ this is the core engine of hivemod .. localised space (no-time) modelling interpolation 
  #\\ note: time is not being modelled and treated independently 
  #\\      .. you had better have enough data in each time slice

  if (is.null(phi)) phi=p$hivemod_phi
  if (is.null(nu)) nu=p$hivemod_nu # nu=0.5 an exponential covariance

  rY = range( x[,p$variables$Y], na.rm=TRUE)

  x_r = range(x[,p$variables$LOCS[1]])
  x_c = range(x[,p$variables$LOCS[2]])

  nr = trunc( diff(x_r)/p$pres ) + 1
  nc = trunc( diff(x_c)/p$pres ) + 1

  x_plons = seq( x_r[1], x_r[2], length.out=nr )
  x_plats = seq( x_c[1], x_c[2], length.out=nc )

  x_locs = expand.grid( x_plons, x_plats ) # final output grid
  attr( x_locs , "out.attrs") = NULL
  names( x_locs ) = p$variables$LOCS

  x$mean = NA
  pa$mean = NA
  pa$sd = NA

  dx = dy = p$pres

  nr2 = 2 * nr
  nc2 = 2 * nc

  dgrid = make.surface.grid(list((1:nr2) * dx, (1:nc2) * dy))
  center = matrix(c((dx * nr2)/2, (dy * nc2)/2), nrow = 1, 
      ncol = 2)

  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1

  # first pass with the global params to get closest fit to data 
  AC = stationary.cov( dgrid, center, Covariance="Matern", range=p$hivemod_phi, nu=p$hivemod_nu )
  mAC = as.surface(dgrid, c(AC))$z
  fW = fft(mAC)/(fft(mC) * nr2 * nc2)

  rm(dgrid, AC,  mAC, mC); gc()

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$variables)) {
      xi = which( x[, p$variables$TIME]==p$ts[ti]) 
    } else {
      xi =1:nrow(x) 
    } 
    
    # map of row, col indices of input data in the new (output) coordinate system
    x_id = cbind( trunc(( x[xi,p$variables$LOCS[1]]-x_r[1])/p$pres) + 1, 
                  trunc(( x[xi,p$variables$LOCS[2]]-x_c[1])/p$pres) + 1 )
    xxii = array_map( "2->1", x_id, c(nr2, nc2) )
    
    # counts
    mN = matrix(0, nrow = nr2, ncol = nc2)
    # mN[xxii] = tapply( rep(1, length(xxii)), INDEX=xxii, FUN=sum, na.rm=TRUE )
    mN[xxii] = 1 # uniform weights .. more stable .. weights cause floating point over/underflow issues ..
    mN[!is.finite(mN)] = 0
    
    # density
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[xxii] = x[xi,p$variables$Y] # fill with data in correct locations
    mY[!is.finite(mY)] = 0
    
    # estimates based upon a global nu,phi .. they will fit to the immediate area near data and so retain their structure
    fN = Re(fft(fft(mN) * fW, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * fW, inverse = TRUE))[1:nr,1:nc]
    Z = fY/fN
    iZ = which( !is.finite( Z))
    if (length(iZ) > 0) Z[iZ] = NA
    lb = which( Z < rY[1] )
    if (length(lb) > 0) Z[lb] = NA
    ub = which( Z > rY[2] )
    if (length(ub) > 0) Z[ub] = NA
    # image(Z)

    # match prediction to input data 
    x$mean[xi] = Z[xxii]
    ss = try( lm( x$mean[xi] ~ x[xi,p$variables$Y], na.action=na.omit) )
    if ( "try-error" %in% class( ss ) ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$hivemod_rsquared_threshold ) next()
    
    if (exists("TIME", p$variables) ) {
      pa_i =  which( pa[, p$variables$TIME]==p$ts[ti])
    } else {
      pa_i = 1:nrow(pa)
    }

    Z_i = cbind( ( pa[pa_i,p$variables$LOCS[1]]-x_r[1])/p$pres + 1, 
                  (pa[pa_i,p$variables$LOCS[2]]-x_c[1])/p$pres + 1 )

    # make sure predictions exist .. kernel density can stop prediction beyond a given range if the xwidth/ywidth options are not used and/or the kernel distance (theta) is small 
    if ( any( Z_i<1) ) next()  
    if ( any( Z_i[,1] > nr) ) next()
    if ( any( Z_i[,2] > nc) ) next()
    pa$mean[pa_i] = Z[Z_i]
    pa$sd[pa_i] = 1

  }

  # plot(mean ~ z , x)
  ss = lm( x$mean ~ x[,p$variables$Y], na.action=na.omit )
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$hivemod_rsquared_threshold ) return(NULL)

  hivemod_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, hivemod_stats=hivemod_stats ) )  
}

