
lbm__splancs = function( p, dat, pa, phi=NULL, nu=NULL ) {
  #\\ this is the core engine of lbm .. localised space (no-time) modelling interpolation 

  dat[, p$variables$Y] = p$lbm_local_family$linkfun ( dat[, p$variables$Y] ) 


  pa$mean = p$lbm_local_family$linkinv( pa$mean )
  # pa$sd   = p$lbm_local_family$linkinv( pa$sd )

  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}

