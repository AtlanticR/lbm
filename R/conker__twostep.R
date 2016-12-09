
conker__twostep = function( p, x, pa, mx=NULL ) {
  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation

  # step 1 -- timeseries modelling
  # use all available data in 'x' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
  # currently only a GAM is enable for the TS component

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

  if (is.null(mx)) mx=pa

  preds = try( predict( hmod, newdata=mx, type="response", se.fit=TRUE ) ) # should already be in the fit so just take the fitted values?

  reject = which( preds$se.fit > quantile( preds$se.fit, probs= p$conker_quantile_bounds[2], na.rm=TRUE ) 
                | preds$fit > p$qs[2] 
                | preds$fit < p$qs[1] )

  preds$fit[reject] = NA

  mx$mean = as.vector( preds$fit )
  mx$sd = as.vector( preds$se.fit )

  mx_r = range(mx[,p$variables$LOCS[1]])
  mx_c = range(mx[,p$variables$LOCS[2]])
  
  mx_nr = diff(mx_r)/p$pres + 1
  mx_nc = diff(mx_c)/p$pres + 1


  # step 2 :: spatial modelling
  Z_all = cbind( ( pa[,p$variables$LOCS[1]]-mx_r[1])/p$pres + 1, 
                (pa[,p$variables$LOCS[2]]-mx_c[1])/p$pres + 1 )

  M_all = cbind( ( mx[,p$variables$LOCS[1]]-mx_r[1])/p$pres + 1, 
                (mx[,p$variables$LOCS[2]]-mx_c[1])/p$pres + 1 )

  # default in case there is no time (a single time slice)
  pa_i = 1:nrow(pa)
  mx_i = 1:nrow(mx)

  for ( ti in 1:p$nt ) {
  
    if ( exists("TIME", p$variables) ) {
      pa_i =  which( pa[, p$variables$TIME]==p$ts[ti] ) 
      mx_i =  which( mx[, p$variables$TIME]==p$ts[ti] ) 
    } 

    if ( any( M_all[ mx_i,] < 1) ) next()  
    if ( any( M_all[ mx_i,1] > mx_nr) ) next()
    if ( any( M_all[ mx_i,2] > mx_nc) ) next()

    # matrix representation of the output surface
    Z = try( smooth.2d( Y=mx[mx_i,"mean"], x=mx[mx_i,p$variables$LOCS], ncol=mx_nc, nrow=mx_nr, theta=p$conker_theta, cov.function=stationary.cov, Covariance="Exponential" ) )
    if ( "try-error" %in% class(Z) ) next()

    iZ = which( !is.finite( Z$z))
    if (length(iZ) > 0) Z$z[iZ] = NA
    rY = range( x[ ,p$variables$Y ], na.rm=TRUE)
    nZ = which( Z$z < rY[1] )
    if (length(nZ) > 0) Z$z[nZ] = NA
    mZ = which( Z$z > rY[2] )
    if (length(mZ) > 0) Z$z[mZ] = NA
    # x11(); image.plot(Z)

    # make sure predictions exist .. kernel density can stop prediction beyond a given range if the xwidth/ywidth options are not used and/or the kernel distance (theta) is small 

    Zsd = try( smooth.2d( Y=mx[mx_i,"sd"], x=mx[mx_i,p$variables$LOCS], ncol=mx_nc, nrow=mx_nr, theta=p$conker_theta, cov.function=stationary.cov, Covariance="Exponential" ) )
    if ( "try-error" %in% class(Zsd) ) next()

    pa$mean[pa_i] = Z$z[Z_all[ pa_i, ]]
    pa$sd[pa_i] = Zsd$z[Z_all[ pa_i, ]]

  }
  
  rm (mx); gc()

  # plot(mean ~ z , x)
  rsquared = ss$r.sq 
  conker_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, conker_stats=conker_stats ) )  
}

