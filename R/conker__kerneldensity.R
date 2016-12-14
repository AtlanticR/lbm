
lstfilter__kerneldensity = function( p, x, pa, nu=NULL, phi=NULL ) {

  #\\ this is the core engine of lstfilter .. localised space (no-time) modelling interpolation 
  #\\ note: time is not being modelled and treated independently 
  #\\      .. you had better have enough data in each time slice


  if (is.null(phi)) phi=p$lstfilter_phi
  if (is.null(nu)) nu=p$lstfilter_nu # nu=0.5 an exponential covariance

  x_r = range(x[,p$variables$LOCS[1]])
  x_c = range(x[,p$variables$LOCS[2]])

  nr = diff(x_r)/p$pres + 1
  nc = diff(x_c)/p$pres + 1

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
  AC = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
    
  mAC = matrix(c(AC), nrow = nr2, ncol = nc2) # or .. mAC = as.surface(dgrid, c(AC))$z
  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1
  fW = fft(mAC)/(fft(mC) * nr2 * nc2)
  rm(dgrid, AC, mAC, mC); gc()

  rY = range( x[,p$variables$Y], na.rm=TRUE)


  for ( ti in 1:p$nt ) {
     
    if ( exists("TIME", p$variables)) {
      xi = which( x[, p$variables$TIME]==p$ts[ti]) 
    } else {
      xi =1:nrow(x) 
    } 
    
    # map of row, col indices of input data in the new (output) coordinate system
    x_id = cbind( (x[xi,p$variables$LOCS[1]]-x_r[1])/p$pres + 1, 
                  (x[xi,p$variables$LOCS[2]]-x_c[1])/p$pres + 1 )
    xxii = array_map( "2->1", x_id, c(nr2, nc2) )
    # counts
    mW = matrix(0, nrow = nr2, ncol = nc2)
    mW[xxii] = tapply( rep(1, length(xi)), INDEX=xxii, FUN=sum, na.rm=TRUE )
    mW[!is.finite(mW)] = 0
    fN = Re(fft(fft(mW) * fW, inverse = TRUE))[1:nr,1:nc]
   
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[xxii] = x[xi,p$variables$Y] # fill with data in correct locations
    mY[!is.finite(mY)] = 0
    fY = Re(fft(fft(mY) * fW, inverse = TRUE))[1:nr,1:nc]
    
    Z = fY/fN

    iZ = which( !is.finite( Z))
    if (length(iZ) > 0) Z[iZ] = NA
    lb = which( Z < rY[1] )
    if (length(lb) > 0) Z[lb] = NA
    ub = which( Z > rY[2] )
    if (length(ub) > 0) Z[ub] = NA
    
    # image(Z)

    # matrix representation of the output surface
     # Z = try( fields::image.smooth( mY, dx=p$pres, dy=p$pres, theta=p$lstfilter_phi)$z ) # phi==theta
  
    if (0) {
      # more control of covariance function .. but not behaving very well and slow .. better to copy internal and strip it down .. TODO
      Z = try( smooth.2d( Y=x[xi,p$variables$Y], x=x[xi,p$variables$LOCS], ncol=nc, nrow=nr, range=phi, nu=nu, cov.function=stationary.cov, Covariance="Exponential" ) )
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
    x$mean[xi] = Z[xxii]
    ss = lm( x$mean[xi] ~ x[xi,p$variables$Y], na.action=na.omit)
    if ( "try-error" %in% class( ss ) ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$lstfilter_rsquared_threshold ) next()
    
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
  if (rsquared < p$lstfilter_rsquared_threshold ) return(NULL)

  lstfilter_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, lstfilter_stats=lstfilter_stats ) )  
}

