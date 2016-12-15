
lstfilter_interpolate_fast = function( ip=NULL, p ) {
  #// designed to be called from lstfilter_interpolate
  #// for the sake of speed and parallelization, the kernel density method is written out again .. it is taken from fields::smooth.2d 

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  phi = p$lstfilter_phi # want phi to be small to retain local structure
  nu = p$lstfilter_nu

  P = lstfilter_attach( p$storage.backend, p$ptr$P )
  Psd = lstfilter_attach( p$storage.backend, p$ptr$Psd )
  Ploc = lstfilter_attach( p$storage.backend, p$ptr$Ploc )
  
  Z2P = as.matrix( cbind( Ploc[,1]-p$plons[1], Ploc[,2]-p$plats[1] ) /p$pres + 1 ) # row, col indices in matrix form Z

  dx = dy = p$pres
  nr = p$nplons
  nc = p$nplats
  nr2 = 2 * nr
  nc2 = 2 * nc
  zp = array_map( "2->1", Z2P, c(nr2, nc2) )

  dgrid = make.surface.grid(list((1:nr2) * dx, (1:nc2) * dy))
  center = matrix(c((dx * nr2)/2, (dy * nc2)/2), nrow = 1, 
      ncol = 2)

  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1

  # first pass with the global params to get closest fit to data 
  AC_global = stationary.cov( dgrid, center, Covariance="Matern", range=p$lstfilter_phi, nu=p$lstfilter_nu )
  mAC_global = matrix(c(AC_global), nrow = nr2, ncol = nc2) # or .. mAC = as.surface(dgrid, c(AC))$z
  fW_global = fft(mAC_global)/(fft(mC) * nr2 * nc2)

  # second pass with local fits to data to smooth what can be smoothed
  AC_local  = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
  mAC_local = matrix(c(AC_local), nrow = nr2, ncol = nc2) # or .. mAC = as.surface(dgrid, c(AC))$z
  fW_local = fft(mAC_local)/(fft(mC) * nr2 * nc2)

  rm(dgrid, AC_global, AC_local, mC, mAC_local, mAC_global); gc()

  for ( iip in ip ) {

    ww = p$runs[ iip, "time_index" ]

    tofill = which( ! is.finite( P[,ww] ) )
    if (length( tofill) > 0 ) {
      # counts
      mN = matrix(0, nrow = nr2, ncol = nc2)
      mN[zp] = tapply( rep(1, length(zp)), INDEX=zp, FUN=sum, na.rm=TRUE )
      mN[!is.finite(mN)] = 0
      
      # density
      mY = matrix(0, nrow = nr2, ncol = nc2)
      mY[zp] = P[,ww] # fill with data in correct locations
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

      # image(Z)
      Z[ Z>p$qs[2] ]=NA
      Z[ Z<p$qs[1] ]=NA
      P[,ww][tofill] = Z[zp][ tofill]
    }


    ## SD estimates
    tofill = which( ! is.finite( Psd[,ww] ) )
    if (length( tofill) > 0 ) {

      # counts
      mN = matrix(0, nrow = nr2, ncol = nc2)
      mN[zp] = tapply( rep(1, length(zp)), INDEX=zp, FUN=sum, na.rm=TRUE )
      mN[!is.finite(mN)] = 0
      
      # density
      mY = matrix(0, nrow = nr2, ncol = nc2)
      mY[zp] = Psd[,ww] # fill with data in correct locations
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

      # image(Z)

      Z[ Z>p$qs[2] ]=NA
      Z[ Z<0 ]=NA
      Psd[,ww][tofill] = Z[zp][ tofill]
    }


  }

  return( "complete" )

}

