
distance_matern = function( phi, nu=0.5, cor=0.95, dmax=phi*10, nX=1000, eps=1e-3 ) {   
  #\\ estimate distance at which spatial (matern) correlations drops to a given threshold cor
  phi = max( eps, phi, na.rm=TRUE )
  nu  = max( eps, nu, na.rm=TRUE )
  dmax = max( eps, dmax, na.rm=TRUE )
  u = matrix(0, ncol=3, nrow=nX )
  u[,1] = seq(0, dmax, length.out=nX )
  u[,2] = 1-fields::Matern( u[,1], range=phi, nu=nu )
  distance = approx( u[,c(2,1)], xout=cor )$y
  return(distance)
}

