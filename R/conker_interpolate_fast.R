
conker_interpolate_fast = function( ip=NULL, p ) {
  #// designed to be called from conker_interpolate
  #// for the sake of speed and parallelization, the kernel density method is written out again .. it is taken from fields::smooth.2d 

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  P = conker_attach( p$storage.backend, p$ptr$P )
  Psd = conker_attach( p$storage.backend, p$ptr$Psd )
  Ploc = conker_attach( p$storage.backend, p$ptr$Ploc )
  
  Z2P = as.matrix( cbind( Ploc[,1]-p$plons[1], Ploc[,2]-p$plats[1] ) /p$pres + 1 ) # row, col indices in matrix form Z

  dx = dy = p$pres
  nr = p$nplons
  nc = p$nplats
  nr2 = 2 * nr
  nc2 = 2 * nc
  zp = array_map( "2->1", Z2P, c(nr2, nc2) )

  theta = min( p$conker_theta, p$pres ) # want theta to be small to retain local structure
  dgrid = make.surface.grid(list((1:nr2) * dx, (1:nc2) * dy))
  center = matrix(c((dx * nr2)/2, (dy * nc2)/2), nrow = 1, 
      ncol = 2)
  AC = stationary.cov( dgrid, center, Covariance="Matern", theta=theta, smoothness=nu )
  mAC = matrix(c(AC), nrow = nr2, ncol = nc2) # or .. mAC = as.surface(dgrid, c(AC))$z
  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1
  fW = fft(mAC)/(fft(mC) * nr2 * nc2)
  rm(dgrid, AC, mAC, mC); gc()

  for ( iip in ip ) {

    ww = p$runs[ iip, "time_index" ]

    # counts
    mW = matrix(0, nrow = nr2, ncol = nc2)
    mW[zp] = tapply( rep(1, length(zp)), INDEX=zp, FUN=sum, na.rm=TRUE )
    mW[!is.finite(mW)] = 0
    fN = Re(fft(fft(mW) * fW, inverse = TRUE))[1:nr,1:nc]

    # mean estimates    
    tofill = which( ! is.finite( P[,ww] ) )
    if (length( tofill) > 0 ) {
      mY = matrix(0, nrow = nr2, ncol = nc2)
      mY[zp] = P[,ww] # fill with data in correct locations
      mY[!is.finite(mY)] = 0
      fY = Re(fft(fft(mY) * fW, inverse = TRUE))[1:nr,1:nc]
      Z = fY/fN
      Z[ Z>p$qs[2] ]=NA
      Z[ Z<p$qs[1] ]=NA
      P[,ww][tofill] = Z[zp][ tofill]
    }

    ## SD estimates
    tofill = which( ! is.finite( Psd[,ww] ) )
    if (length( tofill) > 0 ) {
      mY = matrix(0, nrow = nr2, ncol = nc2)
      mY[zp] = Psd[,ww] # fill with data in correct locations
      mY[!is.finite(mY)] = 0
      fY = Re(fft(fft(mY) * fW, inverse = TRUE))[1:nr,1:nc]
      Z = fY/fN
      Z[ Z>p$qs[2] ]=NA
      Z[ Z<0 ]=NA
      Psd[,ww][tofill] = Z[zp][ tofill]
    }


  }

  return( "complete" )

}

