
lstfilter__twostep = function( p, x, pa, px=NULL, nu=NULL, phi=NULL ) {

  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation
  #\\ nu is the bessel smooth param

  # step 1 -- timeseries modelling
  # use all available data in 'x' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
  # currently only a GAM is enable for the TS component

  if (is.null(phi)) phi=p$lstfilter_phi # range parameter
  if (is.null(nu)) nu=p$lstfilter_nu  # this is an exponential covariance

  rY = quantile( x[,p$variables$Y], probs=p$lstfilter_quantile_bounds, na.rm=TRUE)

  if ( exists("lstfilter_local_model_distanceweighted", p) ) {
    if (p$lstfilter_local_model_distanceweighted) {
      hmod = try( gam( p$lstfilter_local_modelformula, data=x, weights=weights, optimizer=c("outer","optim")  ) )
    } else {
      hmod = try( gam( p$lstfilter_local_modelformula, data=x, optimizer=c("outer","optim")  ) )
    }
  } else {
      hmod = try( gam( p$lstfilter_local_modelformula, data=x ) )
  } 

  if ( "try-error" %in% class(hmod) ) return( NULL )

  ss = summary(hmod)
  if (ss$r.sq < p$lstfilter_rsquared_threshold ) return(NULL)

  if (is.null(px)) px=pa

  preds = try( predict( hmod, newdata=px, type="response", se.fit=TRUE ) ) # should already be in the fit so just take the fitted values?

  reject = which( preds$se.fit > quantile( preds$se.fit, probs= p$lstfilter_quantile_bounds[2], na.rm=TRUE ) 
                | preds$fit > rY[2] 
                | preds$fit < rY[1] )

  if (length(reject) > 0) preds$fit[reject] = NA

  px$mean = as.vector( preds$fit )
  px$sd = as.vector( preds$se.fit )

  px_r = range(px[,p$variables$LOCS[1]], na.rm=TRUE)
  px_c = range(px[,p$variables$LOCS[2]], na.rm=TRUE)
  
  nr = diff(px_r)/p$pres + 1
  nc = diff(px_c)/p$pres + 1

  # step 2 :: spatial modelling
  Z_all = cbind( ( pa[,p$variables$LOCS[1]]-px_r[1])/p$pres + 1, 
                (pa[,p$variables$LOCS[2]]-px_c[1])/p$pres + 1 )

  M_all = cbind( ( px[,p$variables$LOCS[1]]-px_r[1])/p$pres + 1, 
                (px[,p$variables$LOCS[2]]-px_c[1])/p$pres + 1 )

  # default in case there is no time (a single time slice)
  pa_i = 1:nrow(pa)
  px_i = 1:nrow(px)

  dx = dy = p$pres

  nr2 = 2 * nr
  nc2 = 2 * nc

  dgrid = make.surface.grid(list((1:nr2) * dx, (1:nc2) * dy))
  center = matrix(c((dx * nr2)/2, (dy * nc2)/2), nrow = 1, 
      ncol = 2)

  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1

  # first pass with the global params to get closest fit to data 
  AC_global = stationary.cov( dgrid, center, Covariance="Matern", range=p$lstfilter_phi, nu=p$lstfilter_nu )
  mAC_global = as.surface(dgrid, c(AC_global))$z
  fW_global = fft(mAC_global)/(fft(mC) * nr2 * nc2)

  # second pass with local fits to data to smooth what can be smoothed
  AC_local  = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
  mAC_local = as.surface(dgrid, c(AC_local))$z
  fW_local = fft(mAC_local)/(fft(mC) * nr2 * nc2)

  rm(dgrid, AC_global, AC_local, mAC_global, mAC_local, mC); gc()

  for ( ti in 1:p$nt ) {
  
    if ( exists("TIME", p$variables) ) {
      pa_i =  which( pa[, p$variables$TIME]==p$ts[ti] ) 
      px_i =  which( px[, p$variables$TIME]==p$ts[ti] ) 
    } 

    if ( any( M_all[ px_i,] < 1) ) next()  
    if ( any( M_all[ px_i,1] > nr) ) next()
    if ( any( M_all[ px_i,2] > nc) ) next()

    
    # matrix representation of the output surface
    # Z = try( smooth.2d( Y=px[px_i,"mean"], x=px[px_i,p$variables$LOCS], nrow=nr, ncol=nc, dx=p$pres, dy=p$pres, range=phi, cov.function=stationary.cov, Covariance="Matern", nu=nu ) )
    
    xi = cbind( (px[px_i,p$variables$LOCS[1]]-px_r[1])/p$pres + 1, 
                  (px[px_i,p$variables$LOCS[2]]-px_c[1])/p$pres + 1 )
    xxii = array_map( "2->1", xi, c(nr2, nc2) )
    

    # counts
    mN = matrix(0, nrow = nr2, ncol = nc2)
    mN[xxii] = tapply( rep(1, length(px_i)), INDEX=xxii, FUN=sum, na.rm=TRUE )
    mN[!is.finite(mN)] = 0
    
    # density
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[xxii] = px[px_i, "mean"] # fill with data in correct locations
    mY[!is.finite(mY)] = 0
    
    # estimates based upon a global nu,phi .. they will fit to the immediate area near data and so retain their structure
    fN = Re(fft(fft(mN) * fW_global, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * fW_global, inverse = TRUE))[1:nr,1:nc]
    Z = fY/fN
    iZ = which( !is.finite( Z))
    if (length(iZ) > 0) Z[iZ] = NA
    lb = which( Z < rY[1] )
    if (length(lb) > 0) Z[lb] = NA
    ub = which( Z > rY[2] )
    if (length(ub) > 0) Z[ub] = NA
    # image(Z)

    # estimates based upon local nu, phi .. this will over-smooth so if comes as a second pass 
    # to fill in areas with no data (e.g., far away from data locations)
    fN = Re(fft(fft(mN) * fW_local, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * fW_local, inverse = TRUE))[1:nr,1:nc]
    Z_local = fY/fN
    iZ = which( !is.finite( Z_local))
    if (length(iZ) > 0) Z_local[iZ] = NA
    lb = which( Z_local < rY[1] )
    if (length(lb) > 0) Z_local[lb] = NA
    ub = which( Z_local > rY[2] )
    if (length(ub) > 0) Z_local[ub] = NA

    toreplace = which(!is.finite(Z)) 
    if (length(toreplace) > 0 )  Z[toreplace] = Z_local[toreplace]
    pa$mean[pa_i] = Z[Z_all[ pa_i, ]]
    
    # Zsd = try( smooth.2d( Y=px[px_i,"sd"], x=px[px_i,p$variables$LOCS], nrow=nr, ncol=nc, dx=p$pres, dy=p$pres, range=phi, cov.function=stationary.cov, Covariance="Matern", nu=nu ) )
    
    # sd
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[xxii] = px[px_i,"sd"] # fill with data in correct locations
    mY[!is.finite(mY)] = 0
    
    # estimates based upon a global nu,phi .. they will fit to the immediate area near data and so retain their structure
    fN = Re(fft(fft(mN) * fW_global, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * fW_global, inverse = TRUE))[1:nr,1:nc]
    Z = fY/fN
    iZ = which( !is.finite( Z))
    if (length(iZ) > 0) Z[iZ] = NA
    lb = which( Z < 0 )
    if (length(lb) > 0) Z[lb] = NA
    ub = which( Z > rY[2] )
    if (length(ub) > 0) Z[ub] = NA
    # image(Z)

    # estimates based upon local nu, phi .. this will over-smooth so if comes as a second pass 
    # to fill in areas with no data (e.g., far away from data locations)
    fN = Re(fft(fft(mN) * fW_local, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * fW_local, inverse = TRUE))[1:nr,1:nc]
    Z_local = fY/fN
    iZ = which( !is.finite( Z_local))
    if (length(iZ) > 0) Z_local[iZ] = NA
    lb = which( Z_local < 0 )
    if (length(lb) > 0) Z_local[lb] = NA
    ub = which( Z_local > rY[2] )
    if (length(ub) > 0) Z_local[ub] = NA

    toreplace = which(!is.finite(Z)) 
    if (length(toreplace) > 0 )  Z[toreplace] = Z_local[toreplace]
    pa$sd[pa_i] = Z[Z_all[ pa_i, ]]
  }
  
  rm (px); gc()

  # plot(mean ~ z , x)
  rsquared = ss$r.sq 
  lstfilter_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, lstfilter_stats=lstfilter_stats ) )  
}

