
lbm__twostep = function( p, dat, pa, nu=NULL, phi=NULL, varObs=varObs, varSpatial=varSpatial ) {

  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation
  #\\ nu is the bessel smooth param

  # step 1 -- timeseries modelling
  # use all available data in 'dat' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
     # some methods require a uniform (temporal with associated covariates) prediction grid based upon all dat locations 

  px = dat # only the static parts .. time has to be a uniform grid so reconstruct below

  ids = array_map( "xy->1", px[, c("plon", "plat")], gridparams=p$gridparams ) # 100X faster than paste / merge
  todrop = which(duplicated( ids) )
  if (length(todrop>0)) px = px[-todrop,]
  rm(ids, todrop)

  # static vars .. don't need to look up
  tokeep = c(p$variables$LOCS )
  if (exists("weights", dat) ) tokeep = c(tokeep, "weights")
  if (p$nloccov > 0) {
    for (ci in 1:p$nloccov) {
      vn = p$variables$local_cov[ci]
      pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
      nts = ncol(pu)
      if ( nts==1 ) tokeep = c(tokeep, vn ) 
    }
  }
  px = px[ , tokeep ]
  px_n = nrow(px)
  nts = vn = NULL

  # add temporal grid
  if ( exists("TIME", p$variables) ) {
    px = cbind( px[ rep.int(1:px_n, p$nt), ], 
                    rep.int(p$prediction.ts, rep(px_n, p$nt )) )
    names(px)[ ncol(px) ] = p$variables$TIME 
    px = cbind( px, lbm_timecovars ( vars=p$variables$local_all, ti=px[,p$variables$TIME]  ) )
  }

  if (p$nloccov > 0) {
    # add time-varying covars .. not necessary except when covars are modelled locally
    for (ci in 1:p$nloccov) {
      vn = p$variables$local_cov[ci]
      pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
      nts = ncol(pu)
      if ( nts== 1) {
        # static vars are retained in the previous step
      } else if ( nts == p$ny )  {
        px$iy = px$yr - p$yrs[1] + 1 #yr index
        px[,vn] = pu[ cbind(px$i, px$iy) ]  
       } else if ( nts == p$nt) {
        px$it = p$nw*(px$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
        px[,vn] = pu[ cbind(px$i, px$it) ]  
      }
    } # end for loop
    nts = vn = NULL
  } # end if

  ts_gam = lbm__gam( p, dat, px ) # currently only a GAM is enabled for the TS component

  if (is.null( ts_gam)) return(NULL)
  if (ts_gam$lbm_stats$rsquared < p$lbm_rsquared_threshold ) return(NULL)

  # range checks
  rY = range( dat[,p$variables$Y], na.rm=TRUE)
  toosmall = which( ts_gam$predictions$mean < rY[1] )
  toolarge = which( ts_gam$predictions$mean > rY[2] )
  if (length(toosmall) > 0) ts_gam$predictions$mean[toosmall] =  rY[1]  
  if (length(toolarge) > 0) ts_gam$predictions$mean[toolarge] =  rY[2]
 
  pxts = ts_gam$predictions
  # revert to response scale as the following expects this:
  pxts$mean = p$lbm_local_family$linkinv( pxts$mean )
  pxts$sd = p$lbm_local_family$linkinv( pxts$sd )

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
    out = lbm__tps( p, dat=pxts, pa=pa, lambda=varObs/varSpatial  )  
    if (is.null( out)) return(NULL)
  }

  if (p$lbm_twostep_space %in% c("fft", "lowpass", "spatial.process", "lowpass_spatial.process") ) {
    out = lbm__fft( p, dat=pxts, pa=pa, nu=nu, phi=phi )  
    if (is.null( out)) return(NULL)
  }


  return( out )  
}

