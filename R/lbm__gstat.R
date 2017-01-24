
lbm__gstat = function( p, x, pa, nu, phi, varObs, varSpatial ) {
  #\\ this is the core engine of lbm .. localised space (no-time) modelling interpolation 
  #\\ note: time is not being modelled and treated independently 
  #\\      .. you had better have enough data in each time slice ..  essentially this is kriging 
  
  if (!exists( "lbm_gstat_formula", p)) p$lbm_gstat_formula = formula( paste( p$variables$Y, "~ 1 ")) 

  approx_range = phi*sqrt( 8*nu)

  sdTotal = sd(x[,p$variable$Y], na.rm=T)

  x$mean = NA
  pa$mean = NA
  pa$sd = sdTotal  # leave as this as sd estimation is too expensive

  for ( ti in 1:p$nt ) {
    
    if ( exists("TIME", p$variables) ) {
      xi = which( x[ , p$variables$TIME ] == p$prediction.ts[ti] )
      pa_i = which( pa[, p$variables$TIME]==p$prediction.ts[ti])
    } else {
      xi = 1:nrow(x) # all data as p$nt==1
      pa_i = 1:nrow(pa)
    }
    xy = x[xi, p$variables$LOCS]
    z = x[xi, p$variables$Y]

    vMod0 = vgm(psill=varSpatial, model="Mat", range=phi, nugget=varObs, kappa=nu ) # starting model parameters
    gs = gstat(id = "hmk", formula=p$lbm_gstat_formula, locations=~plon+plat, data=xy[xi,], maxdist=approx_range, nmin=p$n.min, nmax=p$n.max, force=TRUE, model=vMod0 )
    # this step adds a lot of time .. 
    preds = predict(gs, newdata=xy[xi,] )
    x$mean[xi] = as.vector( preds[,1] )
    ss = lm( x$mean[xi] ~ x[xi,p$variables$Y], na.action=na.omit)
    if ( "try-error" %in% class( ss ) ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$lbm_rsquared_threshold ) next()
    gsp = predict(gs, newdata=pa[pa_i,] ) # slow for large n
    pa$mean[pa_i] = as.vector(gsp[,1] )
    pa$sd[pa_i]   = as.vector(gsp[,2] )

  }

  # plot(pred ~ z , x)
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  ss = lm( x$mean ~ x[,p$variables$Y], na.action=na.omit)
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$lbm_rsquared_threshold ) return(NULL)

  lbm_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}

