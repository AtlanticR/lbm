
hivemod__bayesx = function( p, x, pa ) {
  #\\ this is the core engine of hivemod .. localised space-time modelling interpolation and prediction .. using bayesx 
   
  # EG: see: bayesx.term.options( bs="kr", method="REML" )
  
  #  logzinc ~  sx( x,y, nu=1.5, bs="kr")  # "kr" is perhaps overly smooth  ..  ie guassian process  .. kriging
  #  logzinc ~  sx( x,y, bs="te")  # more detail .. "te" is preferred

  if ( !exists( "hivemod_local_model_bayesxmethod", p) ) p$hivemod_local_model_bayesxmethod="MCMC"  # slightly more smoothing than the REML method
  if ( !exists( "hivemod_local_family_bayesx", p) ) p$hivemod_local_family_bayesx="gaussian"

  hmod = try( bayesx( p$hivemod_local_modelformula, data=x, method=p$hivemod_local_model_bayesxmethod, 
                     family=p$hivemod_local_family_bayesx ) )

  if ( "try-error" %in% class(hmod) ) return( NULL )

  px = predict(hmod)
  ss = summary(lm( px ~ x[, p$variables$Y ], na.action="na.omit" ))
  if (ss$r.squared < p$hivemod_rsquared_threshold ) return(NULL)
    
  out = try( predict( hmod, newdata=pa, type="response" ) ) 


  # plot( hmod, term = "sx(x,y)", image=TRUE)
  # summary(hmod)
  # lattice::levelplot( out ~ plon+plat, data=k, aspect="iso" )

  if (!inherits(out, "try-error")) return( NULL )

  pa$mean = as.vector( out )
  pa$sd = 1 # no option right now to estim posterior prediction errors .. may be possible with type="terms" but would be slow to simulate  and do not know how to do it yet .. please fix this ..
  varSpatial = hmod$smooth.hyp[,"Variance"]
  varObs = hmod$fixed.effects[1,"Std. Error"]  
  nu = 1.5
  phi=1/hmod$smooth.hyp[,"Smooth Par."] 
  # range = geoR::practicalRange("matern", phi=phi, kappa=nu  )
  range = distance_matern(phi=phi, nu=nu  )

  hivemod_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=ss$r.squared, ndata=nrow(x),
    sdSpatial=sqrt(varSpatial), sdObs=sqrt(varObs), phi=phi, nu=nu, range=range ) 

  # lattice::levelplot( mean ~ plon + plat, data=pa[pa$tiyr==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  
  return( list( predictions=pa, hivemod_stats=hivemod_stats ) )  
}
