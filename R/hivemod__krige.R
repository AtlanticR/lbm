
hivemod__krige = function( p, x, pa, nu, phi, method="default" ) {
  #\\ this is the core engine of hivemod .. localised space (no-time) modelling interpolation 
  # \ as a 2D gaussian process (basically, simple krigimg or TPS -- time is treated as being independent)
  #\\ note: time is not being modelled and treated independently 
  #\\      .. you had better have enough data in each time slice ..  essentially this is kriging 

  
  x$mean = NA
  pa$mean = NA
  pa$sd = NA

  for ( ti in 1:p$nt ) {
    
    if ( exists("TIME", p$variables) ) {
      xi = which( x[ , p$variables$TIME ] == p$ts[ti] )
    } else {
      xi = 1:nrow(x) # all data as p$nt==1
    }

    xy = x[xi, p$variables$LOCS]
    z = x[xi, p$variables$Y]

    if (method %in% c("default", "gstat") ) {
      xy$z = x[xi, p$variables$Y]
      names(xy) = c("plon", "plat", "z")
      vMod0 = vgm(psill=0.75, model="Mat", range=phi, nugget=0.25, kappa=nu ) # starting model parameters
      gs = gstat(id = "hmk", formula = z~1, locations=~plon+plat, data=xy, maxdist=maxdist, nmin=10, force=TRUE, model=vMod0 )

      preds <- predict(gs, newdata=pa )
      # spplot(preds)
    }

    if (method %in% c("fields") ) {
      fsp = try( MLESpatialProcess(xy, z, cov.function=p$fields.cov.function, cov.args=p$fields.cov.args ,
        theta.grid=p$phi.grid, lambda.grid=p$lambda.grid, ngrid = 10, niter = 15, tol = 0.01, 
        Distance = "rdist", nstep.cv = 50 ) )

      if (inherits(fsp, "try-error") )  next()
      if ( fsp$converge != 0 ) next()

      fspmodel <- try( Krig( xy, z, cov.function=p$fields.cov.function, cov.args=p$fields.cov.args, 
        theta=fsp$pars["theta"], lambda=fsp$pars["lambda"] ) )
      if (inherits(fspmodel, "try-error") )  next()
      x$mean[xi] = as.vector( predict(fspmodel, x=x[xi, p$variables$LOCS] ) )
      ss = lm( x$mean[xi] ~ x[xi,p$variables$Y], na.action=na.omit)
      if ( "try-error" %in% class( ss ) ) next()
    }

    rsquared = summary(ss)$r.squared
    if (rsquared < p$hivemod_rsquared_threshold ) next()

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
  if (rsquared < p$hivemod_rsquared_threshold ) return(NULL)

  # TODO:: add some more stats: eg. range estimates, nugget/sill, etc..

  hivemod_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, hivemod_stats=hivemod_stats ) )  
}

