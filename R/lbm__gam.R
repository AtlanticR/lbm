
lbm__gam = function( p, x, pa ) {
  #\\ this is the core engine of lbm .. localised space-time modelling interpolation and prediction
  #\\ simple GAM with spatial weights (inverse distance squared) and ts harmonics 
  
  if ( exists("lbm_local_model_distanceweighted", p) ) {
    if (p$lbm_local_model_distanceweighted) {
      hmod = try( gam( p$lbm_local_modelformula, data=x, weights=weights, optimizer=p$lbm_gam_optimizer) )
    } else {
      hmod = try( gam( p$lbm_local_modelformula, data=x, optimizer=p$lbm_gam_optimizer  ) )
    }
  } else {
      hmod = try( gam( p$lbm_local_modelformula, data=x, optimizer=c("outer", "bfgs")  ) )
  } 


  if ( "try-error" %in% class(hmod) ) return( NULL )

  ss = summary(hmod)
  if (ss$r.sq < p$lbm_rsquared_threshold ) return(NULL)
    
  out = try( predict( hmod, newdata=pa, type="response", se.fit=T ) ) 

  if ( "try-error" %in% class( out ) ) return( NULL )

  pa$mean = as.vector(out$fit)
  pa$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

  lbm_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=ss$r.sq, ndata=ss$n ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa[pa$tiyr==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  
  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}
