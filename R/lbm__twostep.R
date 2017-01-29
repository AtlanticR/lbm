
lbm__twostep = function( p, x, pa, px=NULL, nu=NULL, phi=NULL, varObs=varObs, varSpatial=varSpatial ) {

  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation
  #\\ nu is the bessel smooth param

  # step 1 -- timeseries modelling
  # use all available data in 'x' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
  if (is.null(px)) px=pa

  ts_gam = lbm__gam( p, x, px ) # currently only a GAM is enabled for the TS component

  if (is.null( ts_gam)) return(NULL)
  if (ts_gam$lbm_stats$rsquared < p$lbm_rsquared_threshold ) return(NULL)


  # range checks
  rY = range( x[,p$variables$Y], na.rm=TRUE)
  toosmall = which( ts_gam$predictions$mean < rY[1] )
  toolarge = which( ts_gam$predictions$mean > rY[2] )
  if (length(toosmall) > 0) ts_gam$predictions$mean[toosmall] = rY[1]   
  if (length(toolarge) > 0) ts_gam$predictions$mean[toolarge] = rY[2]   

  px = ts_gam$predictions
  names(px)[which(names(px)=="mean")] = p$variables$Y
  names(px)[which(names(px)=="sd")] = paste(p$variables$Y, "sd", sep=".")
  
  ts_gam = NULL
  gc()

  out = NULL

  # step 2 :: spatial modelling
  if (!exists( "lbm_fft_filter", p)) p$lbm_fft_filter="krige" # default
  
  if ( p$lbm_fft_filter == "krige" ) {
    out = lbm__krige( p, x=px, pa=pa, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial ) 
    if (is.null( out)) return(NULL)
  }

  if (p$lbm_fft_filter %in% c("lowpass", "spatial.process", "lowpass_spatial.process") ) {
    out = lbm__fft( p, x=px, pa=pa, nu=nu, phi=phi )  ## px vs pa ... fix this
    if (is.null( out)) return(NULL)
  }

  return( out )  
}

