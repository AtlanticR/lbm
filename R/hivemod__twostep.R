
hivemod__twostep = function( p, x, pa, px=NULL, nu=NULL, phi=NULL ) {

  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation
  #\\ nu is the bessel smooth param

  # step 1 -- timeseries modelling
  # use all available data in 'x' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
  # currently only a GAM is enable for the TS component

  if (is.null(phi)) phi=p$hivemod_lowpass_phi # range parameter
  if (is.null(nu)) nu=p$hivemod_lowpass_nu  # this is an exponential covariance

  rY = range( x[,p$variables$Y], na.rm=TRUE)

  if ( exists("hivemod_local_model_distanceweighted", p) ) {
    if (p$hivemod_local_model_distanceweighted) {
      hmod = try( gam( p$hivemod_local_modelformula, data=x, weights=weights, optimizer=c("outer","optim")  ) )
    } else {
      hmod = try( gam( p$hivemod_local_modelformula, data=x, optimizer=c("outer","optim")  ) )
    }
  } else {
      hmod = try( gam( p$hivemod_local_modelformula, data=x ) )
  } 

  if ( "try-error" %in% class(hmod) ) return( NULL )

  ss = summary(hmod)
  if (ss$r.sq < p$hivemod_rsquared_threshold ) return(NULL)

  if (is.null(px)) px=pa

  preds = try( predict( hmod, newdata=px, type="response", se.fit=TRUE ) ) # should already be in the fit so just take the fitted values?

  reject = which( preds$se.fit > quantile( preds$se.fit, probs= p$hivemod_quantile_bounds[2], na.rm=TRUE ) 
                | preds$fit > rY[2] 
                | preds$fit < rY[1] )

  if (length(reject) > 0) preds$fit[reject] = NA

  px$mean = as.vector( preds$fit )
  px$sd = as.vector( preds$se.fit )

  px_r = range(px[,p$variables$LOCS[1]], na.rm=TRUE)
  px_c = range(px[,p$variables$LOCS[2]], na.rm=TRUE)
  
  nr = trunc( diff(px_r)/p$pres) + 1
  nc = trunc( diff(px_c)/p$pres) + 1

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
  AC = stationary.cov( dgrid, center, Covariance="Matern", range=p$hivemod_lowpass_phi, nu=p$hivemod_lowpass_nu )
  mAC = as.surface(dgrid, c(AC))$z
  fW = fft(mAC)/(fft(mC) * nr2 * nc2)

  rm(dgrid, AC, mAC, mC); gc()

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
    xxii = array_map( "2->1", trunc(xi), c(nr2, nc2) )

    # counts
    mN = matrix(0, nrow = nr2, ncol = nc2)
    # mN[xxii] = tapply( rep(1, length(px_i)), INDEX=xxii, FUN=sum, na.rm=TRUE )
    mN[xxii] = 1  # counts cause numerical under/over flow
    mN[!is.finite(mN)] = 0
    
    # density
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[xxii] = px[px_i, "mean"] # fill with data in correct locations
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

     pa$mean[pa_i] = Z[Z_all[ pa_i, ]]
    
    # Zsd = try( smooth.2d( Y=px[px_i,"sd"], x=px[px_i,p$variables$LOCS], nrow=nr, ncol=nc, dx=p$pres, dy=p$pres, range=phi, cov.function=stationary.cov, Covariance="Matern", nu=nu ) )
    
    # sd
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[xxii] = px[px_i,"sd"] # fill with data in correct locations
    mY[!is.finite(mY)] = 0
    
    # estimates based upon a global nu,phi .. they will fit to the immediate area near data and so retain their structure
    fN = Re(fft(fft(mN) * fW, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * fW, inverse = TRUE))[1:nr,1:nc]
    Z = fY/fN
    iZ = which( !is.finite( Z))
    if (length(iZ) > 0) Z[iZ] = NA
    lb = which( Z < 0 )
    if (length(lb) > 0) Z[lb] = NA
    ub = which( Z > rY[2] )
    if (length(ub) > 0) Z[ub] = NA
    # image(Z)
  }
  
  rm (px); gc()

  # plot(mean ~ z , x)
  rsquared = ss$r.sq 
  hivemod_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, hivemod_stats=hivemod_stats ) )  
}

