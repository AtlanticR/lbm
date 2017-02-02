
lbm__twostep = function( p, dat, pa, px=NULL, nu=NULL, phi=NULL, varObs=varObs, varSpatial=varSpatial ) {

  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation
  #\\ nu is the bessel smooth param

  # step 1 -- timeseries modelling
  # use all available data in 'dat' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
  if (is.null(px)) px=pa

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

  out = NULL

  # step 2 :: spatial modelling .. essentially a time-space separable solution
  if (!exists( "lbm_twostep_space", p)) p$lbm_twostep_space="krige" # default
  
  if ( p$lbm_twostep_space == "krige" ) {
    out = lbm__krige( p, dat=pxts, pa=pa, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial ) 
    if (is.null( out)) return(NULL)
  }

  if ( p$lbm_twostep_space == "gstat" ) {
    out = lbm__gstat( p, dat=pxts, pa=pa, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial ) 
    if (is.null( out)) return(NULL)
  }

  if (p$lbm_twostep_space %in% c("tps") ) {
    out = lbm__tps( p, dat=pxts, pa=pa, phi=phi, lambda=varObs/varSpatial  )  
    if (is.null( out)) return(NULL)
  }

  if (p$lbm_twostep_space %in% c("fft", "lowpass", "spatial.process", "lowpass_spatial.process") ) {
    out = lbm__fft( p, dat=pxts, pa=pa, nu=nu, phi=phi )  
    if (is.null( out)) return(NULL)
  }


  return( out )  
}

