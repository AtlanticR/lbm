
lbm__spate = function( p, dat, pa, sloc, px=NULL, ws=NULL ) {
  #\\ SPDE solution via FFT using the spate library
  # require(spate)
  # based upon two-step process

  message( "This is not yet finished/debugged .. ")
  
  # step 1 -- timeseries modelling
  # use all available data in 'dat' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
  # currently only a GAM is enable for the TS component

  # ws =  p$lbm_distance_prediction
  #  ws= round( lbm_distance_cur / p$pres )
  #  sloc=Sloc[Si,]

  sdTotal=sd(dat[,p$variable$Y], na.rm=T)

  if (is.null(px)) {
    px = pa
    ws = lbm_distance_cur
  }

  ts_gam = lbm__gam( p, dat, px ) # currently only a GAM is enabled for the TS component

  if (is.null( ts_gam)) return(NULL)
  if (ts_gam$lbm_stats$rsquared < p$lbm_rsquared_threshold ) return(NULL)

  # range checks
  # rY = range( dat[,p$variables$Y], na.rm=TRUE)
  # toosmall = which( ts_gam$predictions$mean < rY[1] )
  # toolarge = which( ts_gam$predictions$mean > rY[2] )
  # if (length(toosmall) > 0) ts_gam$predictions$mean[toosmall] = NA   # permit space modelling to fill this in
  # if (length(toolarge) > 0) ts_gam$predictions$mean[toolarge] = NA   
 
  pxts = ts_gam$predictions
  names(pxts)[which(names(pxts)=="mean")] = p$variables$Y
  names(pxts)[which(names(pxts)=="sd")] = paste(p$variables$Y, "sd", sep=".")
  
  ts_gam = NULL
  gc()

  pxts = pxtsts
  nsq = 2*ws +1 
  adims = c(p$nt, nsq, nsq ) 

  xM = array( NA, dim=adims )
  pxts_id = cbind( 
    round( ( pxts[,p$variables$TIME ] - p$prediction.ts[1] ) / p$tres) + 1,
    round( ws + (pxts[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
    round( ws + (pxts[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1 
  )
  xM[pxts_id] = pxts[,"mean"]  
  xM2 = matrix( xM, nrow=p$nt )
  g = spate.mcmc( y=xM2, n=nsq-1 ) 
  # plot(g, postProcess=TRUE)


  pxts_id = array_map( "3->1", m=round( cbind( 
    ( ws + (pxts[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
    ( ws + (pxts[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1, 
    round( ( pxts[,p$variables$TIME ] - p$prediction.ts[1] ) / p$tres) + 1 )),
    n=adims )
  
  xM = array( NA, dim=adims )
  xM[pxts_id] = pxts[,"mean"]  # this needs to be an even matrix ???

  mm = array_map( "1->2" , pxts_id, n )
  xM2 = array( NA, dims=c(nsq^2, p$nt))
  xM2[mm] = xM[pxts_id] 

  g = spate.mcmc( y=xM2, n=nsq ) 

  # plot(g, postProcess=TRUE)

  predict <- spate.predict(y=mm, tPred=(1:p$nt), spateMCMC=g, Nsim=500, BurnIn=10, DataModel="Normal" )
  Pmean <- apply(predict, c(1,2), mean)
  Psd <- apply(predict, c(1,2), sd)

  # plot(mean ~ z , dat)
  ss = lm( dat$mean ~ dat[,p$variables$Y], na.action=na.omit )
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$lbm_rsquared_threshold ) return(NULL)

  lbm_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}

