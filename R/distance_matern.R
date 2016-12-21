
distance_matern = function( phi, nu=0.5, cor=0.95, dmax=phi*10, nX=1000 ) {   
  #\\ estimate distance at which spatial (matern) correlations drops to a given threshold cor
  u = matrix(0, ncol=3, nrow=nX )
  u[,1] = seq(0, dmax, length.out=nX )
  u[,2] = 1-fields::Matern( u[,1], range=phi, nu=nu )
  distance = approx( u[,c(2,1)], xout=cor )$y
  return(distance)
}

