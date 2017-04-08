

lbm_interpolate = function( ip=NULL, p, debug=FALSE ) {
  #\\ core function to intepolate (model and predict) in parllel

  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) ) 
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  #---------------------
  # data for modelling
  S = lbm_attach( p$storage.backend, p$ptr$S )
  Sflag = lbm_attach( p$storage.backend, p$ptr$Sflag )
  
  Sloc = lbm_attach( p$storage.backend, p$ptr$Sloc )
  Ploc = lbm_attach( p$storage.backend, p$ptr$Ploc )
  Yloc = lbm_attach( p$storage.backend, p$ptr$Yloc )

  Y = lbm_attach( p$storage.backend, p$ptr$Y )

  P = lbm_attach( p$storage.backend, p$ptr$P )
  Pn = lbm_attach( p$storage.backend, p$ptr$Pn )
  Psd = lbm_attach( p$storage.backend, p$ptr$Psd )


  if (p$nloccov > 0) {
    Ycov = lbm_attach( p$storage.backend, p$ptr$Ycov )
  }
  
  if ( exists("TIME", p$variables) ) {
    Ytime = lbm_attach( p$storage.backend, p$ptr$Ytime )
  }

  if ( p$storage.backend != "bigmemory.ram" ) {
    # force copy into RAM to reduce thrashing ?
    # Sloc = Sloc[]
    # Yloc = Yloc[]
    # Y = Y[]
  }

  Yi = lbm_attach( p$storage.backend, p$ptr$Yi )
  # Yi = as.vector(Yi[])  #force copy to RAM as a vector

  # misc intermediate calcs to be done outside of parallel loops
  upsampling = sort( p$sampling[ which( p$sampling > 1 ) ] )
  upsampling = upsampling[ which(upsampling*p$lbm_distance_scale <= p$lbm_distance_max )]
  downsampling = sort( p$sampling[ which( p$sampling < 1) ] , decreasing=TRUE )
  downsampling = downsampling[ which(downsampling*p$lbm_distance_scale >= p$lbm_distance_min )]

  localcount = -1 

  stime = Sys.time()

# main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    localcount = localcount + 1 
    if (( localcount %% 37 )== 0) {
      varstoout = c("n.total", "n.shallow", "n.todo", "n.skipped", "n.outsidepreddomain", "n.outside", "n.complete", "prop_incomp" )
      header = paste( c( varstoout) )
      currentstatus = lbm_db( p=p, DS="statistics.status" )
      currentstatus = c( unlist( currentstatus[ varstoout ] ) )
      dtime = difftime( Sys.time(), stime )
      dtimehr = difftime( Sys.time(), stime, units="hours" )
      nrate = currentstatus["n.complete"]/ as.numeric(dtimehr)
      tmore = currentstatus["n.todo"] / nrate
      tall = (currentstatus["n.todo"]+currentstatus["n.complete"]) / nrate
      cat( paste( "---", p$project.root, p$variables$Y, p$spatial.domain, "--- \n\n"), file=p$lbm_current_status, append=FALSE )
      cat( paste( "Start time :  ", stime, "\n"), file=p$lbm_current_status, append=TRUE )
      cat( paste( "Current time :", Sys.time(), "\n"), file=p$lbm_current_status, append=TRUE )
      cat( paste( "Elapsed time :", format(dtime), "\n" ), file=p$lbm_current_status, append=TRUE)
      cat( paste( "Est. time remaining (hrs) :", round( tmore,3), "\n" ), file=p$lbm_current_status, append=TRUE)
      cat( paste( "Est. time total (hrs) :", round( tall,3), "\n" ), file=p$lbm_current_status, append=TRUE)
      for ( hd in varstoout ){
        cat( paste( hd, ":", currentstatus[hd], "\n" ), file=p$lbm_current_status, append=TRUE)
      }
    }

    Si = p$runs[ iip, "locs" ]

    # Sflag: 
    #   0=TODO, 1=complete, 9=problem, 2=oustide bounds(if any), 3=shallow(if z is a covariate) 
    if ( Sflag[Si] != 0L ) next() 
    Sflag[Si] = 9L   # mark as skipped here. if not it is over-written below 
    print( iip )

    # find data nearest S[Si,] and with sufficient data
    dlon = abs( Sloc[Si,1] - Yloc[Yi[],1] ) 
    dlat = abs( Sloc[Si,2] - Yloc[Yi[],2] ) 
    U =  which( (dlon  <= p$lbm_distance_scale)  & (dlat <= p$lbm_distance_scale) )
    lbm_distance_cur = p$lbm_distance_scale
    ndata = length(U)

    if (0) {
      plot( Sloc[,], pch=20, cex=0.5, col="gray")
      points( Yloc[,], pch=20, cex=0.2, col="green")
      points( Yloc[U,], pch=20, cex=1, col="yellow" )
      points( Sloc[Si,2] ~ Sloc[Si,1], pch=20, cex=5, col="blue" )
    }

    o = ores = NULL

    if (ndata > p$n.min ) {
      if (ndata < p$n.max ) {
        # nothing to do
      } else {
        if ( ndata <= p$n.max * 1.5 ) { 
          # if close to p$n.max, subsample quickly 
          U = U[ .Internal( sample( length(U), p$n.max, replace=FALSE, prob=NULL)) ] 
          ndata = p$n.max
        } else {
          # need to downsample
          for ( dsamp in downsampling )  { # lots of data .. downsample
            lbm_distance_cur = p$lbm_distance_scale * dsamp
            U = which( dlon < lbm_distance_cur & dlat < lbm_distance_cur )# faster to take a block 
            ndata = length(U)
            if ( ndata <= p$n.max ) break()
            if ( lbm_distance_cur <= p$lbm_distance_min ) {
              # reached lower limit in distance, taking a subsample instead
              U = which( dlon < p$lbm_distance_min & dlat < p$lbm_distance_min ) # faster to take a block 
              U = U[ .Internal( sample( length(U), p$n.max, replace=FALSE, prob=NULL)) ]
              ndata = length(U)
              break()
            }
          }
        }
      }
    } else {
      # need to upsample
      for ( usamp in upsampling )  {
        lbm_distance_cur = p$lbm_distance_scale * usamp
        U = which( dlon < lbm_distance_cur & dlat < lbm_distance_cur ) # faster to take a block 
        ndata = length(U)
        if ( ndata >= p$n.min ) {
          if (ndata >= p$n.max) {
            U = U[ .Internal( sample( length(U), p$n.max, replace=FALSE, prob=NULL)) ] 
            ndata = p$n.max
            break()
          }
        }
      }
    }

   if (ndata < p$n.min)  next() # check in case a fault in logic, above

   # crude mean variogram across all time slices
    o = try( lbm_variogram( xy=Yloc[U,], z=Y[U], methods=p$lbm_variogram_method, family=p$lbm_local_family ) )
      if ( !is.null(o)) {
        if (!inherits(o, "try-error")) {
          if (exists(p$lbm_variogram_method, o)) {
            ores = o[[p$lbm_variogram_method]] # store current best estimate of variogram characteristics
            if (exists("range", ores) ) {
              if (is.finite(ores[["range"]] )) {
                if ( (ores[["range"]] > p$lbm_distance_min) & (ores[["range"]] <= p$lbm_distance_max) ) {
                  lbm_distance_cur = ores[["range"]]
                  vario_U  = which( dlon  <= ores[["range"]]  & dlat <= ores[["range"]] )
                  vario_ndata =length(vario_U)                
                  if ((vario_ndata > p$n.min) & (vario_ndata < p$n.max) ) { 
                    U  = vario_U
                    ndata = vario_ndata
                  }  
                }
              }
            } 
          }
        }   
      }

    if (is.null(ores)) {
      ndata = U = o = ores = NULL
      next()
    }

    if (ndata > p$n.max & ndata <= p$n.max * 1.5) {
      U = U[ .Internal( sample( length(U), p$n.max, replace=FALSE, prob=NULL)) ] 
      ndata = p$n.max
    } 

    if (ndata < p$n.min)  next() # check in case a fault in logic, above
    if (ndata > p$n.max)  next() # check in case a fault in logic, above

    dlon=dlat=o=NULL; gc()

    if (debug) {
      print( ores )
    }

    YiU = Yi[U]  
    # So, YiU and p$lbm_distance_prediction determine the data entering into local model construction
    # dist_model = lbm_distance_cur

    pa = lbm_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=p$windowsize.half )
    if (is.null(pa)) next()

      if (debug) {
        # check that position indices are working properly
        Sloc = lbm_attach( p$storage.backend, p$ptr$Sloc )
        Yloc = lbm_attach( p$storage.backend, p$ptr$Yloc )
        plot( Yloc[U,2]~ Yloc[U,1], col="red", pch=".", 
          ylim=range(c(Yloc[U,2], Sloc[Si,2], Ploc[pa$i,2]) ), 
          xlim=range(c(Yloc[U,1], Sloc[Si,1], Ploc[pa$i,1]) ) ) # all data
        points( Yloc[YiU,2] ~ Yloc[YiU,1], col="green" )  # with covars and no other data issues
        points( Sloc[Si,2] ~ Sloc[Si,1], col="blue" ) # statistical locations
        # statistical output locations
        grids= spatial_grid(p, DS="planar.coords" )
        points( grids$plat[round( (Sloc[Si,2]-p$origin[2])/p$pres) + 1] 
              ~ grids$plon[round( (Sloc[Si,1]-p$origin[1])/p$pres) + 1] , col="purple", pch=25, cex=5 ) 

        points( grids$plat[pa$iplat] ~ grids$plon[ pa$iplon] , col="cyan", pch=20, cex=0.01 ) # check on Proc iplat indexing
        points( Ploc[pa$i,2] ~ Ploc[ pa$i, 1] , col="black", pch=20, cex=0.7 ) # check on pa$i indexing -- prediction locations
      }


    # prep dependent data 
    # reconstruct data for modelling (dat) and data for prediction purposes (pa)
    dat = data.frame( Y[YiU] ) # these are residuals if there is a global model
    names(dat) = p$variables$Y
    dat$plon = Yloc[YiU,1]
    dat$plat = Yloc[YiU,2]
    dat$weights = 1
    # dat$weights = 1 / (( Sloc[Si,2] - dat$plat)**2 + (Sloc[Si,1] - dat$plon)**2 )# weight data in space: inverse distance squared
    # dat$weights[ which( dat$weights < 1e-4 ) ] = 1e-4
    # dat$weights[ which( dat$weights > 1 ) ] = 1
    
    if (p$nloccov > 0) {
      for (i in 1:p$nloccov) dat[, p$variables$local_cov[i] ] = Ycov[YiU,i] # no need for other dim checks as this is user provided 
    }
     
    if (exists("TIME", p$variables)) {
      dat[, p$variables$TIME ] = Ytime[YiU,] 
      dat = cbind( dat, lbm_timecovars ( vars=p$variables$local_all, ti=dat[,p$variables$TIME]  ) )
    }

    nu = phi = varSpatial = varObs = NULL
    if ( exists("nu", ores) ) nu = ores$nu
    if ( exists("phi", ores) && ores$phi > (p$pres/2)) phi = ores$phi 
    if ( exists("varSpatial", ores) ) varSpatial = ores$varSpatial
    if ( exists("varObs", ores) ) varObs = ores$varObs
  
    if (is.null(nu)) nu = p$lbm_lowpass_nu
    if (is.null(phi)) phi = lbm_distance_cur/sqrt( 8*nu) # crude estimate of phi based upon current scaling  distance approximates the range at 90% autocorrelation(e.g., see Lindgren et al. 2011)
    if (is.null(varSpatial)) varSpatial =0.5 * var( p$lbm_local_family$linkfun(dat[, p$variables$Y]), na.rm=TRUE)
    if (is.null(varObs)) varObs = varSpatial
    
    # model and prediction .. outputs are in scale of the link (and not response)
    # the following permits user-defined models (might want to use compiler::cmpfun )
    gc()
    res =NULL
    res = try( switch( p$lbm_local_modelengine, 
      bayesx = lbm__bayesx( p, dat, pa ),
      habitat = lbm__habitat( p, dat, pa ), # TODO 
      inla = lbm__inla( p, dat, pa ),
      gam = lbm__gam( p, dat, pa ), 
      gaussianprocess2Dt = lbm__gaussianprocess2Dt( p, dat, pa ), 
      gaussianprocess = lbm__gaussianprocess( p, dat, pa ),  # TODO
      glm = lbm__glm( p, dat, pa), 
      gstat = lbm__gstat( p, dat, pa, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial ), 
      krige = lbm__krige( p, dat, pa, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial ), 
      LaplacesDemon = lbm__LaplacesDemon( p, dat, pa ),
      splancs = lbm__splancs( p, dat, pa ), # TODO
      spate = lbm__spate( p, dat, pa, sloc=Sloc[Si,], distance=lbm_distance_cur, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial), 
      fft = lbm__fft( p, dat, pa, nu=nu, phi=phi ), 
      tps = lbm__tps( p, dat, pa, lambda=varObs/varSpatial ), 
      twostep = lbm__twostep( p, dat, pa, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial )
    ) )

    ### from here on predictions and sd are scaled by family p$lbm_local_family$linkfun ( ) 
    ## must return to response scale with p$lbm_local_family$linkinv( ) in parent call


    if (debug) {
      print( str(res))
    }

    if (0) {
      lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
       
      lattice::levelplot( mean ~ plon + plat, data=res$predictions, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
   
      for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
    }

    rm(dat); gc()
    if ( inherits(res, "try-error") ) {
      dat = pa = NULL
      next()
    }
    
    if ( is.null(res)) {
      dat = pa = NULL
      next()
    }
    if ( all( !is.finite(res$predictions$mean ))) {
      dat = pa = res = NULL
      next()
    }

    if (exists( "lbm_quantile_bounds", p)) {
      tq = quantile( p$lbm_local_family$linkfun(Y[YiU]), probs=p$lbm_quantile_bounds, na.rm=TRUE  )
      toolow  = which( res$predictions$mean < tq[1] )
      toohigh = which( res$predictions$mean > tq[2] )
      if (length( toolow) > 0)  res$predictions$mean[ toolow] = tq[1]
      if (length( toohigh) > 0) res$predictions$mean[ toohigh] = tq[2]
    }
    
    ii = which( is.finite(res$predictions$mean ))
    if (length(ii) < 5) {
      dat = pa = res = NULL
      next()  # looks to be a faulty solution
    }

    # stats collator
    if (!exists("lbm_stats",  res) ) res$lbm_stats = list()
    
    if (!exists("sdSpatial", res$lbm_stats)) {
      # some methods can generate spatial stats simultaneously .. 
      # it is faster to keep them all together instead of repeating here
      # field and RandomFields gaussian processes seem most promising ... 
      # default to fields for speed:
      res$lbm_stats["sdSpatial"] = NA 
      res$lbm_stats["sdObs"] = NA 
      res$lbm_stats["range"] = NA
      res$lbm_stats["phi"] = NA
      res$lbm_stats["nu"] = NA
      if ( !is.null(ores)) {
        res$lbm_stats["sdSpatial"] = sqrt( ores[["varSpatial"]] ) 
        res$lbm_stats["sdObs"] = sqrt(ores[["varObs"]]) 
        res$lbm_stats["range"] = ores[["range"]]
        res$lbm_stats["phi"] = ores[["phi"]]
        res$lbm_stats["nu"] = ores[["nu"]]
      } 
    }
    
    if ( exists("TIME", p$variables) ){
      # annual ts, seasonally centered and spatially 
      # pa_i = which( Sloc[Si,1]==Ploc[,1] & Sloc[Si,2]==Ploc[,2] )
      pac_i = which( res$predictions$plon==Sloc[Si,1] & res$predictions$plat==Sloc[Si,2] )
      # plot( mean~tiyr, res$predictions[pac_i,])
      # plot( mean~tiyr, res$predictions, pch="." )
      res$lbm_stats["ar_timerange"] = NA 
      res$lbm_stats["ar_1"] = NA
            
      if (length(pac_i) > 5) {
        pac = res$predictions[ pac_i, ]
        pac$dyr = pac[, p$variables$TIME] - trunc(pac[, p$variables$TIME] )
        piid = which( zapsmall( pac$dyr - p$dyear_centre) == 0 )
        pac = pac[ piid, c(p$variables$TIME, "mean")]
        pac = pac[ order(pac[,p$variables$TIME]),]
        if (length(piid) > 5 ) {
          ts.stat = NULL
          ts.stat = try( lbm_timeseries( pac$mean, method="fft" ) )
          if (!is.null(ts.stat) && !inherits(ts.stat, "try-error") ) {
            res$lbm_stats["ar_timerange"] = ts.stat$quantilePeriod 
            if (all( is.finite(pac$mean))) {
              afin = which (is.finite(pac$mean) )
              if (length(afin) > 5 && var( pac$mean, na.rm=TRUE) > p$eps ) {
                ar1 = NULL
                ar1 = try( ar( pac$mean, order.max=1 ) )
                if (!inherits(ar1, "try-error")) {
                  if ( length(ar1$ar) == 1 ) {
                    res$lbm_stats["ar_1"] = ar1$ar
                  }  
                } 
              }
            }
            if ( !is.finite(res$lbm_stats[["ar_1"]]) ) {
              ar1 = try( cor( pac$mean[1:(length(piid) - 1)], pac$mean[2:(length(piid))], use="pairwise.complete.obs" ) )
              if (!inherits(ar1, "try-error")) res$lbm_stats["ar_1"] = ar1 
            }
          } 

          ### Do the logistic model here ! -- if not already done ..
          if (!exists("ts_K", res$lbm_stats)) {
            # model as a logistic with ts_r, ts_K, etc .. as stats outputs
            
          } 

        } 
        rm ( pac, piid )
      } 
      rm(pac_i)
    }


    # update SD estimates of predictions with those from other locations via the
    # incremental  method ("online algorithm") of mean estimation after Knuth ;
    # see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    # update means: inverse-variance weighting   
    # see https://en.wikipedia.org/wiki/Inverse-variance_weighting
   
    npred = nrow(res$predictions)

    if ( ! exists("TIME", p$variables) ) {

      u = which( is.finite( P[res$predictions$i] ) )  # these have data already .. update
      if ( length( u ) > 1 ) {
        ui = res$predictions$i[u]  # locations of P to modify
        Pn[ui] = Pn[ui] + 1 # update counts
        stdev_update =  Psd[ui] + ( res$predictions$sd[u] -  Psd[ui] ) / Pn[ui]
        means_update = ( P[ui] / Psd[ui]^2 + res$predictions$mean[u] / res$predictions$sd[u]^2 ) / ( Psd[ui]^(-2) + res$predictions$sd[u]^(-2) )
        mm = which(is.finite( means_update + stdev_update ))
        if( length(mm)> 0) {
          iumm = ui[mm]
          Psd[iumm] = stdev_update[mm]
          P  [iumm] = means_update[mm]
        } 
        stdev_update = NULL
        means_update = NULL

        rm(ui, mm, iumm)
      }

      # first time # no data yet
      v = setdiff(1:npred, u)         
      if ( length(v) > 0 ) {
        vi = res$predictions$i[v]
        Pn [vi] = 1
        P  [vi] = res$predictions$mean[v]
        Psd[vi] = res$predictions$sd[v]
      }
    }

    if ( exists("TIME", p$variables) ) {
      u = which( is.finite( P[res$predictions$i,1] ) )  # these have data already .. update
      u_n = length( u ) 
      if ( u_n > 1 ) {  # ignore if only one point .. mostly because it can cause issues with matrix form .. 
        # locations of P to modify
        ui = sort(unique(res$predictions$i[u]))
        nc = ncol(P)
        if (p$storage.backend == "ff" ) {
          add.ff(Pn, 1, ui, 1:nc ) # same as Pn[ui,] = Pn[ui]+1 but 2X faster
        } else {
          Pn[ui,] = Pn[ui,] + 1
        }
        stdev_update =  Psd[ui,] + ( res$predictions$sd[u] -  Psd[ui,] ) / Pn[ui,]
        means_update = ( P[ui,] / Psd[ui,]^2 + res$predictions$mean[u] / res$predictions$sd[u]^2 ) / 
          ( Psd[ui,]^(-2) + res$predictions$sd[u]^(-2) )
        
        updates = means_update + stdev_update 
        if (!is.matrix(updates)) next()

        mm = which( is.finite( rowSums(updates)))  # created when preds go outside quantile bounds .. this removes all data from a given location rather than the space-time .. severe but likely due to a poor prediction and so remove all (it is also faster this way as few manipulations)
        if( length(mm)> 0) {
          iumm = ui[mm] 
          Psd[iumm,] = stdev_update[mm,]
          P  [iumm,] = means_update[mm,]
          iumm = NULL
        } 
        stdev_update = NULL
        means_update = NULL
        rm(ui, mm)

      }

      # do this as a second pass in case NA's were introduced by the update .. unlikely , but just in case
      v = which( !is.finite( P[res$predictions$i,1] ) )  # these have data already .. update
      nv = length(v)          # no data yet
      if ( nv > 0 ) {
        vi = sort(unique(res$predictions$i[v]))
        Pn [vi,] = 1
        P  [vi,] = res$predictions$mean[v]
        Psd[vi,] = res$predictions$sd[v]
        rm(vi)
      } 
    }

    # save stats
    for ( k in 1: length(p$statsvars) ) {
      if (exists( p$statsvars[k], res$lbm_stats )) {
        S[Si,k] = res$lbm_stats[[ p$statsvars[k] ]]
      }
    }
    
    if (debug) {
        v = res$predictions
        if ( exists("TIME", p$variables) ){
          v = v[which( v[,p$variables$TIME]==2000.55),]
        }
        require(lattice)
        print(
          levelplot( mean ~ plon+plat, v, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
      )
    }
      
      if (0) {
     
        lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
        
        for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      
        for (i in 1:p$nt) {
          print( lattice::levelplot( P[pa$i,i] ~ Ploc[pa$i,1] + Ploc[ pa$i, 2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
        }

      }
   
    res = NULL
    pa = NULL

    # ----------------------
    # do last. it is an indicator of completion of all tasks 
    # restarts would be broken otherwise
    Sflag[Si] = 1L  # mark as done 

  }  # end for loop
  
  invisible()

}

