
lbm__glm = function( p, dat, pa ) {
  #\\ this is the core engine of lbm .. localised space-time modelling interpolation and prediction
  #\\ simple GAM with spatial weights (inverse distance squared) and ts harmonics 

  sdTotal=sd(dat[,p$variable$Y], na.rm=T)

  if ( exists("lbm_local_model_distanceweighted", p) ) {
    if (p$lbm_local_model_distanceweighted) {
      hmod = try( glm( p$lbm_local_modelformula, data=dat, family=p$lbm_local_family, weights=Y_wgt  ) ) 
    } else {
      hmod = try( glm( p$lbm_local_modelformula, data=dat, family=p$lbm_local_family ) ) 
    }
  } else {
      hmod = try( glm( p$lbm_local_modelformula, data=dat, family=p$lbm_local_family  ) ) 
  } 

  if ( "try-error" %in% class(hmod) ) return( NULL )

  ss = summary(hmod)
  rsq = 1 - (ss$deviance/ss$null.deviance)
  if ( rsq < p$lbm_rsquared_threshold ) return(NULL)
  
  out = try( predict( hmod, newdata=pa, type="response", se.fit=TRUE ) ) 

  if ( "try-error" %in% class( out ) ) return( NULL )

  pa$mean = as.vector(out$fit)
  pa$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

  lbm_stats = list( sdTotal=sdTotal, rsquared=rsq, ndata=nrow(dat) ) # must be same order as p$statsvars .. pseudo rsquared for logistic .. for poisson {1- logLik(mod) / logLik(mod_saturated)} might be better
  
  # lattice::levelplot( mean ~ plon + plat, data=pa[pa$tiyr==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  
  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}
