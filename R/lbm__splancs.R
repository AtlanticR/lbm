
lbm__splancs = function( p, dat, pa, phi=NULL, nu=NULL ) {
  #\\ this is the core engine of lbm .. localised space (no-time) modelling interpolation 

  dat[, p$variables$Y] = p$lbm_local_family$linkfun ( dat[, p$variables$Y] ) 

  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}

