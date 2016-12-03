
conker__gaussianprocess = function( p, x, pa ) {
  #\\ this is the core engine of conker .. localised space  and time modelling/ interpolation 
  # \ as a gaussian process
  # TODO 

  # plot(pred ~ z , x)
  ss = lm( x$mean ~ x[,p$variables$Y], na.action=na.omit)
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$conker_rsquared_threshold ) return(NULL)

  conker_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, conker_stats=conker_stats ) )  
}

