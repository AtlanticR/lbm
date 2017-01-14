
lbm__twostep = function( p, x, pa, px=NULL, nu=NULL, phi=NULL ) {

  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation
  #\\ nu is the bessel smooth param

  # step 1 -- timeseries modelling
  # use all available data in 'x' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
  # currently only a GAM is enable for the TS component

  rY = range( x[,p$variables$Y], na.rm=TRUE)

  if ( exists("lbm_local_model_distanceweighted", p) ) {
    if (p$lbm_local_model_distanceweighted) {
      hmod = try( gam( p$lbm_local_modelformula, data=x, weights=weights, optimizer=c("outer","optim")  ) )
    } else {
      hmod = try( gam( p$lbm_local_modelformula, data=x, optimizer=c("outer","optim")  ) )
    }
  } else {
      hmod = try( gam( p$lbm_local_modelformula, data=x ) )
  } 

  if ( "try-error" %in% class(hmod) ) return( NULL )

  ss = summary(hmod)
  if (ss$r.sq < p$lbm_rsquared_threshold ) return(NULL)

  if (is.null(px)) px=pa

  preds = try( predict( hmod, newdata=px, type="response", se.fit=TRUE ) ) # should already be in the fit so just take the fitted values?

  toosmall = which( preds$fit < rY[1] )
  toolarge = which( preds$fit > rY[2] )
  if (length(toosmall) > 0) preds$fit[toosmall] = rY[1]   
  if (length(toolarge) > 0) preds$fit[toolarge] = rY[2]   

  px$mean = as.vector( preds$fit )
  px$sd = as.vector( preds$se.fit )

  px_r = range(px[,p$variables$LOCS[1]], na.rm=TRUE)
  px_c = range(px[,p$variables$LOCS[2]], na.rm=TRUE)
  
  nr = round( diff(px_r)/p$pres) + 1
  nc = round( diff(px_c)/p$pres) + 1

  # step 2 :: spatial modelling

  Z_all = array_map( "xy->2", coords=pa[,p$variables$LOCS], 
    corner=c(px_r[1], px_c[1]), res=c(p$pres, p$pres) )

  M_all = array_map( "xy->2", coords=px[,p$variables$LOCS], 
    corner=c(px_r[1], px_c[1]), res=c(p$pres, p$pres) )

  # default in case there is no time (a single time slice)
  pa_i = 1:nrow(pa)
  px_i = 1:nrow(px)

  dx = dy = p$pres

  nr2 = 2 * nr
  nc2 = 2 * nc

  # constainer for spatial filters
  dgrid = make.surface.grid(list((1:nr2) * dx, (1:nc2) * dy))
  center = matrix(c((dx * nr2)/2, (dy * nc2)/2), nrow = 1, ncol = 2)
  
  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1
  

  if (p$lbm_fft_filter == "lowpass") {
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", range=p$lbm_lowpass_phi, nu=p$lbm_lowpass_nu )
    sp.covar.surf = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar.surf) / ( fft(mC) * nr2 * nc2 )
  }

  if (p$lbm_fft_filter == "spatial.process") {
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
    sp.covar.surf = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar.surf) / ( fft(mC) * nr2 * nc2 )
  }

  if (p$lbm_fft_filter == "lowpass_spatial.process") {
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", range=p$lbm_lowpass_phi, nu=p$lbm_lowpass_nu )
    sp.covar.surf = as.surface(dgrid, c(sp.covar))$z
    sp.covar2 = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
    sp.covar.surf2 = as.surface(dgrid, c(sp.covar2))$z
    sp.covar.kernel = {fft(sp.covar.surf)/ ( fft(mC) * nr2 * nc2 )} * {fft(sp.covar.surf2)/ ( fft(mC) * nr2 * nc2 )}
  }

  sp.covar = sp.covar2 = sp.covar.surf = sp.covar.surf2 = dgrid = center = mC = NULL
  gc()


  for ( ti in 1:p$nt ) {
  
    if ( exists("TIME", p$variables) ) {
      pa_i =  which( pa[, p$variables$TIME]==p$prediction.ts[ti] ) 
      px_i =  which( px[, p$variables$TIME]==p$prediction.ts[ti] ) 
    } 

    if ( any( M_all[ px_i,] < 1) ) next()  
    if ( any( M_all[ px_i,1] > nr) ) next()
    if ( any( M_all[ px_i,2] > nc) ) next()

    # matrix representation of the output surface
    # Z = try( smooth.2d( Y=px[px_i,"mean"], x=px[px_i,p$variables$LOCS], nrow=nr, ncol=nc, dx=p$pres, dy=p$pres, range=phi, cov.function=stationary.cov, Covariance="Matern", nu=nu ) )

    x_id = array_map( "xy->2", coords=px[px_i,p$variables$LOCS], 
      corner=c(px_r[1], px_c[1]), res=c(p$pres, p$pres) )

    u = as.image( px[px_i, "mean"], ind=as.matrix( x_id), na.rm=TRUE, nx=nr, ny=nc )
    
    mN = matrix(0, nrow = nr2, ncol = nc2)
    mN[1:nr,1:nc] = u$weights
    mN[!is.finite(mN)] = 0
    
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[1:nr,1:nc] = u$z
    mY[!is.finite(mY)] = 0
  
    u = NULL

    # low pass filter based upon a global nu,phi .. remove high freq variation
    fN = Re(fft(fft(mN) * sp.covar.kernel, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * sp.covar.kernel, inverse = TRUE))[1:nr,1:nc]
    
    mY = mN = NULL

    Z = fY/fN
    fY = fN = NULL
 
    # image.plot(Z)
    # lb = which( Z < rY[1] )
    # if (length(lb) > 0) Z[lb] = NA
    # ub = which( Z > rY[2] )
    # if (length(ub) > 0) Z[ub] = NA
   
    pa$mean[pa_i] = Z[Z_all[ pa_i, ]]
    
    
    # sd
    u = as.image( px[px_i, "sd"], ind=as.matrix( x_id), na.rm=TRUE, nx=nr, ny=nc )

    mN = matrix(0, nrow = nr2, ncol = nc2)
    mN[1:nr,1:nc] = u$weights
    mN[!is.finite(mN)] = 0

    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[1:nr,1:nc] = u$z #mean sd
    mY[!is.finite(mY)] = 0
    
    u = NULL

    # low pass filter based upon a global nu,phi .. remove high freq variation
    fN = Re(fft(fft(mN) * sp.covar.kernel, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * sp.covar.kernel, inverse = TRUE))[1:nr,1:nc]
    
    mY = mN = NULL

    Z = fY/fN
    fY = fN = NULL

    # lb = which( Z < rY[1] )
    # if (length(lb) > 0) Z[lb] = NA
    # ub = which( Z > rY[2] )
    # if (length(ub) > 0) Z[ub] = NA
    # # image.plot(Z)

    pa$sd[pa_i] = Z[Z_all[ pa_i, ]]
    

  }
  
  rm (px); gc()

  # plot(mean ~ z , x)
  rsquared = ss$r.sq 
  lbm_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}

