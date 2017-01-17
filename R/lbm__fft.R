
lbm__fft = function( p, x, pa, nu=NULL, phi=NULL ) {

  #\\ this is the core engine of lbm .. localised space (no-time) modelling interpolation 
  #\\ note: time is not being modelled and treated independently 
  #\\      .. you had better have enough data in each time slice
  #\\ first a low-pass filter as defined by p$lbm_lowpass_nu, p$lbm_lowpass_phi, then a simple covariance filter determined by nu,phi
  # varObs=varObs, varSpatial=varSpatial

  x_r = range(x[,p$variables$LOCS[1]])
  x_c = range(x[,p$variables$LOCS[2]])

  pa_r = range(pa[,p$variables$LOCS[1]])
  pa_c = range(pa[,p$variables$LOCS[2]])
  
  nr = round( diff(x_r)/p$pres ) + 1
  nc = round( diff(x_c)/p$pres ) + 1

  x_plons = seq( x_r[1], x_r[2], length.out=nr )
  x_plats = seq( x_c[1], x_c[2], length.out=nc )

  x_locs = expand.grid( x_plons, x_plats ) # final output grid
  attr( x_locs , "out.attrs") = NULL
  names( x_locs ) = p$variables$LOCS


  sdTotal = sd(x[,p$variable$Y], na.rm=T)

  x$mean = NA
  pa$mean = NA
  pa$sd = sdTotal  # this is ignored with fft

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
    # both ..
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", range=p$lbm_lowpass_phi, nu=p$lbm_lowpass_nu )
    sp.covar.surf = as.surface(dgrid, c(sp.covar))$z
    sp.covar2 = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
    sp.covar.surf2 = as.surface(dgrid, c(sp.covar2))$z
    sp.covar.kernel = {fft(sp.covar.surf)/ ( fft(mC) * nr2 * nc2 )} * {fft(sp.covar.surf2)/ ( fft(mC) * nr2 * nc2 )
  } }

  sp.covar = sp.covar2 = sp.covar.surf = sp.covar.surf2 = dgrid = center = mC = NULL
  gc()

  xi =1:nrow(x) 
  pa_i = 1:nrow(pa)
  origin=c(x_r[1], x_c[1])
  res=c(p$pres, p$pres)

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$variables)) xi = which( x[, p$variables$TIME]==p$prediction.ts[ti] ) 
    
    # map of row, col indices of input data in the new (output) coordinate system
    
    x_id = array_map( "xy->2", coords=x[xi,p$variables$LOCS], origin=origin, res=res )
    
    u = as.image( x[xi,p$variables$Y], ind=as.matrix( x_id), na.rm=TRUE, nx=nr, ny=nc )
    
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
    
    if (exists("TIME", p$variables) ) pa_i =  which( pa[, p$variables$TIME]==p$prediction.ts[ti])

    Z_i = array_map( "xy->2", coords=pa[pa_i,p$variables$LOCS], origin=origin, res=res )

    # make sure predictions exist
    if ( any( Z_i<1) ) next()  
    if ( any( Z_i[,1] > nr) ) next()
    if ( any( Z_i[,2] > nc) ) next()
    pa$mean[pa_i] = Z[Z_i]
    # pa$sd[pa_i] = NA  ## fix as NA
    Z = NULL
    
  }

  # plot(mean ~ z , x)
  ss = lm( x$mean ~ x[,p$variables$Y], na.action=na.omit )
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$lbm_rsquared_threshold ) return(NULL)

  lbm_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}

