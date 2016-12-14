
lstfilter__gaussianprocess2Dt = function( p, x, pa ) {
  #\\ this is the core engine of lstfilter .. localised space (no-time) modelling interpolation 
  # \ as a 2D gaussian process (basically, simple krigimg or TPS -- time is treated as being independent)
  #\\ note: time is not being modelled and treated independently 
  #\\      .. you had better have enough data in each time slice ..  essentially this is kriging 

  if (!exists("fields.cov.function", p)) p$fields.cov.function="stationary.cov"
  if (!exists("fields.Covariance", p)) p$fields.Covariance="Exponential" # note that "Rad.cov" is TPS
  if (!exists("fields.cov.args", p) & p$fields.Covariance=="Matern") {
    if (!exists("fields.nu", p)) p$fields.nu=0.5  # note: this is the smoothness or shape parameter (fix at 0.5 if not calculated or given -- exponential)   
    p$fields.cov.args=list( Covariance=p$fields.Covariance, smoothness=p$fields.nu ) # this is exponential covariance 
  }
  
  x$mean = NA
  pa$mean = NA
  pa$sd = NA

  theta.grid = 10^seq( -6, 6, by=0.5) * p$lstfilter_distance_scale # aprox magnitude of the phi parameter
  lambda.grid = 10^seq( -9, 3, by=0.5) 

  for ( ti in 1:p$nt ) {
    
    if ( exists("TIME", p$variables) ) {
      xi = which( x[ , p$variables$TIME ] == p$ts[ti] )
    } else {
      xi = 1:nrow(x) # all data as p$nt==1
    }

    xy = x[xi, p$variables$LOCS]
    z = x[xi, p$variables$Y]
    
    fsp = try( MLESpatialProcess(xy, z, cov.function=p$fields.cov.function, cov.args=p$fields.cov.args ,
      theta.grid=theta.grid, lambda.grid=lambda.grid, ngrid = 10, niter = 15, tol = 0.01, 
      Distance = "rdist", nstep.cv = 50 ) )

    if (inherits(fsp, "try-error") )  next()
    if ( fsp$converge != 0 ) next()

    fspmodel <- try( Krig( xy, z, cov.function=p$fields.cov.function, cov.args=p$fields.cov.args, 
      theta=fsp$pars["theta"], lambda=fsp$pars["lambda"] ) )
    if (inherits(fspmodel, "try-error") )  next()

    x$mean[xi] = as.vector( predict(fspmodel, x=x[xi, p$variables$LOCS] ) )
    ss = lm( x$mean[xi] ~ x[xi,p$variables$Y], na.action=na.omit)
    if ( "try-error" %in% class( ss ) ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$lstfilter_rsquared_threshold ) next()

    if ( exists("TIME", p$variables) ) {
      pa_i = which( pa[, p$variables$TIME]==p$ts[ti])
    } else {
      pa_i = 1:nrow(pa)
    }

    pa$mean[pa_i] = predict(fspmodel, x=pa[pa_, p$variables$LOCS] )
    pa$sd[pa_i]   = predictSE(fspmodel, x=pa[pa_, p$variables$LOCS] )

    if ( 0 ){
      # debugging plots
      surface(fspmodel)
      fsp.p<- predictSurface(fspmodel, lambda=fsp$pars["lambda"], nx=200, ny=200, )
      surface(fsp.p, type="I")
      fsp.p2<- predictSurfaceSE(fspmodel)
      surface(fsp.p, type="C")
    }
 
  }

  # plot(pred ~ z , x)
  ss = lm( x$mean ~ x[,p$variables$Y], na.action=na.omit)
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$lstfilter_rsquared_threshold ) return(NULL)

  # TODO:: add some more stats: eg. range estimates, nugget/sill, etc..

  lstfilter_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, lstfilter_stats=lstfilter_stats ) )  
}

