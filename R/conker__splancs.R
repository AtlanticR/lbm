
conker__kerneldensity = function( p, x, pa, smoothness=0.5 ) {
  #\\ this is the core engine of conker .. localised space (no-time) modelling interpolation 
  #\\ note: time is not being modelled and treated independently 
  #\\      .. you had better have enough data in each time slice

  x_r = range(x[,p$variables$LOCS[1]])
  x_c = range(x[,p$variables$LOCS[2]])

  x_nr = diff(x_r)/p$pres + 1
  x_nc = diff(x_c)/p$pres + 1

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
    if ( exists("TIME", p$variables) ) {
      pa_i =  which( pa[, p$variables$TIME]==p$ts[ti] ) 
    } else { 
      pa_i = 1:nrow(pa) 
    }
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

