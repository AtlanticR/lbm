
hivemod__glm = function( p, x, pa ) {
  #\\ this is the core engine of hivemod .. localised space-time modelling interpolation and prediction
  #\\ simple GAM with spatial weights (inverse distance squared) and ts harmonics 

  if ( exists("hivemod_local_model_distanceweighted", p) ) {
    if (p$hivemod_local_model_distanceweighted) {
      hmod = try( glm( p$hivemod_local_modelformula, data=x, weights=Y_wgt  ) ) 
    } else {
      hmod = try( glm( p$hivemod_local_modelformula, data=x  ) ) 
    }
  } else {
      hmod = try( glm( p$hivemod_local_modelformula, data=x  ) ) 
  } 

  if ( "try-error" %in% class(hmod) ) return( NULL )

  ss = summary(hmod)
  rsq = 1 - (ss$deviance/ss$null.deviance)
  if ( rsq < p$hivemod_rsquared_threshold ) return(NULL)
  
  out = try( predict( hmod, newdata=pa, type="response", se.fit=TRUE ) ) 

  if ( "try-error" %in% class( out ) ) return( NULL )

  pa$mean = as.vector(out$fit)
  pa$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

  hivemod_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsq, ndata=nrow(x) ) # must be same order as p$statsvars .. pseudo rsquared for logistic .. for poisson {1- logLik(mod) / logLik(mod_saturated)} might be better
  
  # lattice::levelplot( mean ~ plon + plat, data=pa[pa$tiyr==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  
  return( list( predictions=pa, hivemod_stats=hivemod_stats ) )  
}
