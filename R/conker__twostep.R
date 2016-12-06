
conker__twostep = function( p, x, pa, smoothness=0.5 ) {
  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation

  # step 1 -- timeseries modelling
  # use all available data in 'x' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
  # currently only a GAM

  if ( exists("conker_local_model_distanceweighted", p) ) {
    if (p$conker_local_model_distanceweighted) {
      hmod = try( gam( p$conker_local_modelformula, data=x, weights=weights, optimizer=c("outer","optim")  ) )
    } else {
      hmod = try( gam( p$conker_local_modelformula, data=x, optimizer=c("outer","optim")  ) )
    }
  } else {
      hmod = try( gam( p$conker_local_modelformula, data=x ) )
  } 

  if ( "try-error" %in% class(hmod) ) return( NULL )

  ss = summary(hmod)
  if (ss$r.sq < p$conker_rsquared_threshold ) return(NULL)
    
  preds = try( predict( hmod, newdata=pa, type="response", se.fit=TRUE ) ) # should already be in the fit so just take the fitted values?

  reject = which( preds$se.fit > quantile( preds$se.fit, probs= p$conker_quantile_bounds[2], na.rm=TRUE ) 
                | preds$fit > p$qs[2] 
                | preds$fit < p$qs[1] )

  preds$fit[reject] = NA

  pa$mean = as.vector( preds$fit )
  pa$sd = as.vector( preds$se.fit )

  # locations of the new (output) coord system .. smaller than the data range of x
  pa_r = range(pa[,p$variables$LOCS[1]])
  pa_c = range(pa[,p$variables$LOCS[2]])
  
  pa_nr = diff(pa_r)/p$pres + 1
  pa_nc = diff(pa_c)/p$pres + 1
  
  # step 2 :: spatial modelling
  for ( ti in 1:p$nt ) {
     
    if ( exists("TIME", p$variables) ) {
      pa_i =  which( pa[, p$variables$TIME]==p$ts[ti] ) 
    } else { 
      pa_i = 1:nrow(pa) 
    }

    Z_i = cbind( ( pa[pa_i,p$variables$LOCS[1]]-pa_r[1])/p$pres + 1, 
                  (pa[pa_i,p$variables$LOCS[2]]-pa_c[1])/p$pres + 1 )

    if ( any( Z_i<1) ) next()  
    if ( any( Z_i[,1] > pa_nr) ) next()
    if ( any( Z_i[,2] > pa_nc) ) next()

    # matrix representation of the output surface
    M = matrix( NA, nrow=pa_nr, ncol=pa_nc) 
    M[Z_i] = pa[pa_i, "mean"] # fill with data in correct locations
    Z = try( fields::image.smooth( M, dx=p$pres, dy=p$pres, theta=p$conker_theta)$z )
  
    if (0) {
      # more control of covariance function .. but not behaving very well and slow .. better to copy internal and strip it down .. TODO
      Z = try( smooth.2d( Y=pa[pa_i,"mean"], x=pa[pa_i,p$variables$LOCS], ncol=pa_nc, nrow=pa_nr, theta=p$conker_theta,
        cov.function=stationary.cov, Covariance="Exponential", p=smoothness ) )
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
    # make sure predictions exist .. kernel density can stop prediction beyond a given range if the xwidth/ywidth options are not used and/or the kernel distance (theta) is small 
    pa$mean[pa_i] = Z[Z_i]
  }

  # plot(mean ~ z , x)
  rsquared = ss$r.sq 
  conker_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, conker_stats=conker_stats ) )  
}

