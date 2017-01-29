
lbm__spate = function( p, x, pa, sloc, px=NULL, ws=NULL ) {
  #\\ SPDE solution via FFT using the spate library
  # require(spate)
  # based upon two-step process

  # step 1 -- timeseries modelling
  # use all available data in 'x' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
  # currently only a GAM is enable for the TS component

  #  ws= round( lbm_distance_cur / p$pres )
  #  sloc=Sloc[Si,]

  if ( exists("lbm_local_model_distanceweighted", p) ) {
    if (p$lbm_local_model_distanceweighted) {
      hmod = try( gam( p$lbm_local_modelformula, data=x, weights=weights, optimizer=c("outer","optim")  ) )
    } else {
      hmod = try( gam( p$lbm_local_modelformula, data=x, optimizer=c("outer","optim")  ) )
    }
  } else {
      hmod = try( gam( p$lbm_local_modelformula, data=x ) )
  } 

  if ( "try-error" %in% class(hmod) ) return( NULL )

  ss = summary(hmod)
  if (ss$r.sq < p$lbm_rsquared_threshold ) return(NULL)

  if (is.null(px)) {
    px = pa
    ws = lbm_distance_cur
  }

  # ws =  p$lbm_distance_prediction

  preds = try( predict( hmod, newdata=px, type="response", se.fit=TRUE ) ) # should already be in the fit so just take the fitted values?

  reject = which( preds$se.fit > quantile( preds$se.fit, probs= p$lbm_quantile_bounds[2], na.rm=TRUE ) 
                | preds$fit > p$qs[2] 
                | preds$fit < p$qs[1] )

  preds$fit[reject] = NA

  px$mean = as.vector( preds$fit )
  px$sd = as.vector( preds$se.fit )

  nsq = 2*ws +1 
  adims = c(p$nt, nsq, nsq ) 

  xM = array( NA, dim=adims )
  px_id = cbind( 
    round( ( px[,p$variables$TIME ] - p$prediction.ts[1] ) / p$tres) + 1,
    round( ws + (px[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
    round( ws + (px[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1 
  )
  xM[px_id] = px[,"mean"]  
  xM2 = matrix( xM, nrow=p$nt )
  g = spate.mcmc( y=xM2, n=nsq-1 ) 
  # plot(g, postProcess=TRUE)


  px_id = array_map( "3->1", m=round( cbind( 
    ( ws + (px[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
    ( ws + (px[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1, 
    round( ( px[,p$variables$TIME ] - p$prediction.ts[1] ) / p$tres) + 1 )),
    n=adims )
  
  xM = array( NA, dim=adims )
  xM[px_id] = px[,"mean"]  # this needs to be an even matrix ???

  mm = array_map( "1->2" , px_id, n )
  xM2 = array( NA, dims=c(nsq^2, p$nt))
  xM2[mm] = xM[px_id] 

  g = spate.mcmc( y=xM2, n=nsq ) 

  # plot(g, postProcess=TRUE)

  predict <- spate.predict(y=mm, tPred=(1:p$nt), spateMCMC=g, Nsim=500, BurnIn=10, DataModel="Normal" )
  Pmean <- apply(predict, c(1,2), mean)
  Psd <- apply(predict, c(1,2), sd)

  # plot(mean ~ z , x)
  ss = lm( x$mean ~ x[,p$variables$Y], na.action=na.omit )
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$lbm_rsquared_threshold ) return(NULL)

  lbm_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}

