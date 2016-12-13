
conker__twostep = function( p, x, pa, px=NULL, smoothness=0.5, phi=NULL ) {

  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation
  #\\ nu is the bessel smooth param

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

  if (is.null(px)) px=pa

  preds = try( predict( hmod, newdata=px, type="response", se.fit=TRUE ) ) # should already be in the fit so just take the fitted values?

  reject = which( preds$se.fit > quantile( preds$se.fit, probs= p$conker_quantile_bounds[2], na.rm=TRUE ) 
                | preds$fit > p$qs[2] 
                | preds$fit < p$qs[1] )

  preds$fit[reject] = NA

  px$mean = as.vector( preds$fit )
  px$sd = as.vector( preds$se.fit )

  px_r = range(px[,p$variables$LOCS[1]])
  px_c = range(px[,p$variables$LOCS[2]])
  
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
  AC = stationary.cov( dgrid, center, Covariance="Matern", theta=theta, smoothness=nu )
    
  mAC = matrix(c(AC), nrow = nr2, ncol = nc2) # or .. mAC = as.surface(dgrid, c(AC))$z
  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1
  fW = fft(mAC)/(fft(mC) * nr2 * nc2)
  rm(dgrid, AC, mAC, mC); gc()

  rY = range( x[,p$variables$Y], na.rm=TRUE)

  if (is.null(phi)) phi=p$conker_theta # range parameter
  if (is.null(smoothness)) smoothness=0.5 # this is an exponential covariance


  for ( ti in 1:p$nt ) {
  
    if ( exists("TIME", p$variables) ) {
      pa_i =  which( pa[, p$variables$TIME]==p$ts[ti] ) 
      px_i =  which( px[, p$variables$TIME]==p$ts[ti] ) 
    } 

    if ( any( M_all[ px_i,] < 1) ) next()  
    if ( any( M_all[ px_i,1] > nr) ) next()
    if ( any( M_all[ px_i,2] > nc) ) next()

    
    # matrix representation of the output surface
    # Z = try( smooth.2d( Y=px[px_i,"mean"], x=px[px_i,p$variables$LOCS], nrow=nr, ncol=nc, dx=p$pres, dy=p$pres, theta=theta, cov.function=stationary.cov, Covariance="Matern", smoothness=nu ) )
    
    x_id = cbind( (px[px_i,p$variables$LOCS[1]]-px_r[1])/p$pres + 1, 
                  (px[px_i,p$variables$LOCS[2]]-px_c[1])/p$pres + 1 )
    xxii = array_map( "2->1", x_id, c(nr2, nc2) )
    
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[xxii] = px[px_i,"mean"] # fill with data in correct locations
    mY[!is.finite(mY)] = 0
    fY = Re(fft(fft(mY) * fW, inverse = TRUE))[1:nr,1:nc]
    
    # counts
    mW = matrix(0, nrow = nr2, ncol = nc2)
    mW[xxii] = tapply( rep(1, length(xi)), INDEX=xxii, FUN=sum, na.rm=TRUE )
    mW[!is.finite(mW)] = 0
    fN = Re(fft(fft(mW) * fW, inverse = TRUE))[1:nr,1:nc]
    Z = fY/fN

    if ( "try-error" %in% class(Z) ) next()

    iZ = which( !is.finite( Z))
    if (length(iZ) > 0) Z[iZ] = NA
    lb = which( Z < rY[1] )
    if (length(lb) > 0) Z[lb] = NA
    ub = which( Z > rY[2] )
    if (length(ub) > 0) Z[ub] = NA

    # Zsd = try( smooth.2d( Y=px[px_i,"sd"], x=px[px_i,p$variables$LOCS], nrow=nr, ncol=nc, dx=p$pres, dy=p$pres, theta=theta, cov.function=stationary.cov, Covariance="Matern", smoothness=nu ) )
    
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[xxii] = px[px_i,"sd"] # fill with data in correct locations
    mY[!is.finite(mY)] = 0
    fY = Re(fft(fft(mY) * fW, inverse = TRUE))[1:nr,1:nc]
    # mW, fW already computed above
    fN = Re(fft(fft(mW) * fW, inverse = TRUE))[1:nr,1:nc]
    Z = fY/fN
  

    iZ = which( !is.finite( Z))
    if (length(iZ) > 0) Z[iZ] = NA
    lb = which( Z < 0 )
    if (length(lb) > 0) Z[lb] = NA
    ub = which( Z > rY[2] )
    if (length(ub) > 0) Z[ub] = NA

    # Z = try( smooth.2d( Y=px[px_i,"mean"], x=px[px_i,p$variables$LOCS], ncol=px_nc, nrow=px_nr, cov.function=stationary.cov, Covariance="Matern", smoothness=smoothness, range=phi ) )
    # if ( "try-error" %in% class(Z) ) next()

    # iZ = which( !is.finite( Z$z))
    # if (length(iZ) > 0) Z$z[iZ] = NA
    # rY = range( x[ ,p$variables$Y ], na.rm=TRUE)
    # nZ = which( Z$z < rY[1] )
    # if (length(nZ) > 0) Z$z[nZ] = NA
    # mZ = which( Z$z > rY[2] )
    # if (length(mZ) > 0) Z$z[mZ] = NA
    # # x11(); image.plot(Z)
    # Zsd = try( smooth.2d( Y=px[px_i,"sd"], x=px[px_i,p$variables$LOCS], ncol=px_nc, nrow=px_nr, cov.function=stationary.cov, Covariance="Matern", smoothness=smoothness, range=phi ) )

    if ( "try-error" %in% class(Zsd) ) next()
    pa$mean[pa_i] = Z$z[Z_all[ pa_i, ]]
    pa$sd[pa_i] = Zsd$z[Z_all[ pa_i, ]]
  }
  
  rm (px); gc()

  # plot(mean ~ z , x)
  rsquared = ss$r.sq 
  conker_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, conker_stats=conker_stats ) )  
}

