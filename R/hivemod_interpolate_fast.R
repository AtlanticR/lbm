
hivemod_interpolate_fast = function( ip=NULL, p ) {
  #// designed to be called from hivemod_interpolate
  #// for the sake of speed and parallelization, the kernel density method is written out again .. it is taken from fields::smooth.2d 

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  S = hivemod_attach( p$storage.backend, p$ptr$S )
  
  # do this here as the Stats are from the most reliable estimates
  nu = median(  S[,which( p$statsvars=="nu" )], na.rm=TRUE )
  phi = median(  S[,which( p$statsvars=="phi" )], na.rm=TRUE )

  P = hivemod_attach( p$storage.backend, p$ptr$P )
  Psd = hivemod_attach( p$storage.backend, p$ptr$Psd )
  Ploc = hivemod_attach( p$storage.backend, p$ptr$Ploc )
  
  Z2P = as.matrix( round( cbind( Ploc[,1]-p$plons[1], Ploc[,2]-p$plats[1] ) /p$pres + 1) ) # row, col indices in matrix form Z

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

  # first a low pass filter 
  AC = stationary.cov( dgrid, center, Covariance="Matern", range=p$hivemod_lowpass_phi, nu=p$hivemod_lowpass_nu )
  mAC = as.surface(dgrid, c(AC))$z
  fW = fft(mAC)/(fft(mC) * nr2 * nc2)

  # spatial autocorrelation filter 
  AC = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
  mAC = as.surface(dgrid, c(AC))$z
  fW = fft(mAC)/(fft(mC) * nr2 * nc2)

  rm(dgrid, AC, mC, mAC); gc()

  for ( iip in ip ) {

    ww = p$runs[ iip, "time_index" ]

    tofill = which( ! is.finite( P[,ww] ) )
    if (length( tofill) > 0 ) {
      # counts
      rY = range( P[,ww], na.rm=TRUE )

      mN = matrix(0, nrow = nr2, ncol = nc2)
      mN[zp] = tapply( rep(1, length(zp)), INDEX=zp, FUN=sum, na.rm=TRUE )
      mN[!is.finite(mN)] = 0
      
      # density
      mY = matrix(0, nrow = nr2, ncol = nc2)
      mY[zp] = P[,ww] # fill with data in correct locations
      mY[!is.finite(mY)] = 0
      
      # estimates based upon a global nu,phi .. they will fit to the immediate area near data and so retain their structure
      fN = Re(fft(fft(mN) * fW, inverse = TRUE))[1:nr,1:nc]
      fY = Re(fft(fft(mY) * fW, inverse = TRUE))[1:nr,1:nc]
      Z = fY/fN
      lb = which( Z < rY[1] )
      if (length(lb) > 0) Z[lb] = NA
      ub = which( Z > rY[2] )
      if (length(ub) > 0) Z[ub] = NA

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

      Z[ Z>p$qs[2] ]=NA
      Z[ Z<0 ]=NA
      Psd[,ww][tofill] = Z[zp][ tofill]
    }


  }

  return( "complete" )

}

