
lbm_interpolate = function( ip=NULL, p ) {
  #\\ core function to intepolate (model and predict) in parllel

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  #---------------------
  # data for modelling
  # dependent vars # already link-transformed in lbm_db("dependent")
  
  S = lbm_attach( p$storage.backend, p$ptr$S )
  Sflag = lbm_attach( p$storage.backend, p$ptr$Sflag )
  
  Sloc = lbm_attach( p$storage.backend, p$ptr$Sloc )
  Ploc = lbm_attach( p$storage.backend, p$ptr$Ploc )
  Yloc = lbm_attach( p$storage.backend, p$ptr$Yloc )

  Y = lbm_attach( p$storage.backend, p$ptr$Y )

  P = lbm_attach( p$storage.backend, p$ptr$P )
  Pn = lbm_attach( p$storage.backend, p$ptr$Pn )
  Psd = lbm_attach( p$storage.backend, p$ptr$Psd )

  if (exists("local_cov", p$variables)) {
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

  if (p$lbm_local_modelengine=="habitat") {
    Ylogit = lbm_attach( p$storage.backend, p$ptr$Ylogit )
    Plogit = lbm_attach( p$storage.backend, p$ptr$Plogit )
    Plogitsd = lbm_attach( p$storage.backend, p$ptr$Plogitsd )
  }

  Yi = lbm_attach( p$storage.backend, p$ptr$Yi )
  # Yi = as.vector(Yi[])  #force copy to RAM as a vector

  # misc intermediate calcs to be done outside of parallel loops
  upsampling = sort( p$sampling[ which( p$sampling > 1 ) ] )
  upsampling = upsampling[ which(upsampling*p$lbm_distance_scale <= p$lbm_distance_max )]
  downsampling = sort( p$sampling[ which( p$sampling < 1) ] , decreasing=TRUE )
  downsampling = downsampling[ which(downsampling*p$lbm_distance_scale >= p$lbm_distance_min )]

  am = c(p$nplons, p$nplats)
  ploc_ids = array_map( "2->1", round(cbind(Ploc[,1]-p$plons[1], Ploc[,2]-p$plats[1])/p$pres+1), am )

  localcount = -1 

# main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    localcount = localcount + 1 
    if (( localcount %% 100 )== 0) {
      varstoout = c("n.total", "n.land", "n.todo", "n.problematic", "n.outside", "n.complete", "prop_incomp" )
      header = paste( c( "system.time", varstoout) )
      currentstatus = lbm_db( p=p, DS="statistics.status" )
      currentstatus = c( Sys.time(), unlist( currentstatus[ varstoout ] ) )
      cat( header, file=p$lbm_current_status, append=FALSE)
      cat( paste("\n"), file=p$lbm_current_status, append=TRUE)
      cat( currentstatus, file=p$lbm_current_status, append=TRUE )
    }

    Si = p$runs[ iip, "locs" ]

    # Sflag: 
    #   0=TODO, 1=complete, 9=problem, 2=oustide bounds(if any), 3=land(if z is a covariate) 
    if ( Sflag[Si] != 0L ) next() 
    Sflag[Si] = 9L   # mark as problematic here. if not it is over-written below 
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

    o = try( lbm_variogram( xy=Yloc[U,], z=p$lbm_local_family$linkfun(Y[U]), 
      methods=p$lbm_variogram_method ) )
      if ( !is.null(o)) {
        if (!inherits(o, "try-error")) {
          if (exists(p$lbm_variogram_method, o)) {
            ores = o[[p$lbm_variogram_method]] # store current best estimate of variogram characteristics
            if (!is.null(ores[["range"]] ) ) {
              if ( (ores[["range"]] > p$pres) & (ores[["range"]] <= p$lbm_distance_max) ) {
                vario_U  = which( dlon  <= ores[["range"]]  & dlat <= ores[["range"]] )
                vario_ndata =length(vario_U)                
                if ((vario_ndata > p$n.min) & (vario_ndata < p$n.max) ) { 
                  U  = vario_U
                  ndata = vario_ndata
                  lbm_distance_cur = ores[["range"]]
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

    if (ndata < p$n.min)  next() # check in case a fault in logic, above
    if (ndata > p$n.max) {
      U = U[ .Internal( sample( length(U), p$n.max, replace=FALSE, prob=NULL)) ] 
      ndata = p$n.max
    }

    dlon=dlat=o=NULL; gc()

    YiU = Yi[U]  
    # So, YiU and dist_prediction determine the data entering into local model construction
    # dist_model = lbm_distance_cur

    if ( lbm_distance_cur < p$lbm_distance_prediction ) dist_prediction = p$lbm_distance_prediction # do not predict greater than p$lbm_distance_prediction

    # construct prediction/output grid area ('pa')
    windowsize.half = floor(dist_prediction/p$pres) # convert distance to discretized increments of row/col indices

    pa_w = -windowsize.half : windowsize.half
    pa_w_n = length(pa_w)
    
    iwplon = round( (Sloc[Si,1]-p$plons[1])/p$pres + 1 + pa_w )
    iwplat = round( (Sloc[Si,2]-p$plats[1])/p$pres + 1 + pa_w )
    
    pa = NULL
    pa = data.frame( iplon = rep.int(iwplon, pa_w_n) , 
                     iplat = rep.int(iwplat, rep.int(pa_w_n, pa_w_n)) )
    rm(iwplon, iwplat, pa_w, pa_w_n, windowsize.half)

    bad = which( (pa$iplon < 1 & pa$iplon > p$nplons) | (pa$iplat < 1 & pa$iplat > p$nplats) )
    if (length(bad) > 0 ) pa = pa[-bad,]
    if (nrow(pa)< 5)  {
      bad = pa = YiU = o = U = NULL
      next()
    }

    pa$i = match( array_map( "2->1", cbind(pa$iplon, pa$iplat), am ), ploc_ids )
        
    bad = which( !is.finite(pa$i))
    if (length(bad) > 0 ) pa = pa[-bad,]
    pa_n = nrow(pa)
    if ( pa_n < 5) {
      bad = pa = YiU = o = U = NULL
      next()
    }

      if (0) {
        # check that position indices are working properly
        Sloc = lbm_attach( p$storage.backend, p$ptr$Sloc )
        Yloc = lbm_attach( p$storage.backend, p$ptr$Yloc )
        plot( Yloc[U,2]~ Yloc[U,1], col="red", pch=".", 
          ylim=range(c(Yloc[U,2], Sloc[Si,2], Ploc[pa$i,2]) ), 
          xlim=range(c(Yloc[U,1], Sloc[Si,1], Ploc[pa$i,1]) ) ) # all data
        points( Yloc[YiU,2] ~ Yloc[YiU,1], col="green" )  # with covars and no other data issues
        points( Sloc[Si,2] ~ Sloc[Si,1], col="blue" ) # statistical locations
        # statistical output locations
        
        points( p$plats[round( (Sloc[Si,2]-p$plats[1])/p$pres) + 1] ~ p$plons[round((Sloc[Si,1]-p$plons[1])/p$pres) + 1] , col="purple", pch=25, cex=5 ) 

        points( p$plats[pa$iplat] ~ p$plons[ pa$iplon] , col="cyan", pch=20, cex=0.01 ) # check on Proc iplat indexing
        points( Ploc[pa$i,2] ~ Ploc[ pa$i, 1] , col="black", pch=20, cex=0.7 ) # check on pa$i indexing -- prediction locations
      }
   
    pa$plon = Ploc[ pa$i, 1]
    pa$plat = Ploc[ pa$i, 2]

 
    # prediction covariates i.e., independent variables/ covariates
    pvars = c("plon", "plat", "i")
    if (exists("local_cov", p$variables)) {
      # .. not necessary except when covars are modelled locally
      for (ci in 1:length(p$variables$local_cov)) {
        vn = p$variables$local_cov[ci]
        pu = NULL
        pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
        nts = ncol(pu)
        if ( nts== 1 ) {
          pvars = c( pvars, vn )
          pa[,vn] = pu[pa$i]  # ie. a static variable
        }
      }
    }
    pa = pa[, pvars]

    if ( exists("TIME", p$variables) ) {
      pa = cbind( pa[ rep.int(1:pa_n, p$nt), ], 
                      rep.int(p$ts, rep(pa_n, p$nt )) )
      names(pa) = c( pvars, p$variables$TIME )
      if ( p$variables$TIME != "yr" ) pa$yr = trunc( pa[,p$variables$TIME] )
      # where time exists and there are seasonal components, 
      # additional variables are created/needed here: cos.w, sin.w, etc.. 
      # for harmonic analysis: to add an offset to a trig function (b) must add cos to a sin function
      # y ~ a + c*sin(x+b)
      # y ~ a + c*sin(b)*cos(x) + c*cos(b)*sin(x)  
      #   .. as C*sin(x+b) = C*( cos(b) * sin(x) + sin(b) * cos(x) )
      # y ~ b0 + b1*x1 + b2*x2
      # where: 
      #   a = b0
      #   c^2 = b1^2 + b2^2 = c^2*(sin^2(b) + cos^2(b))
      #   c = sqrt(b1^2 + b2^2)
      #   b1/b2 = tan(b)  
      #   b = arctan(b1/b2)
      if ("dyear" %in% p$variables$local_all)  pa$dyear = pa[, p$variables$TIME] - pa$yr  # fractional year
      if ("cos.w" %in% p$variables$local_all)  pa$cos.w  = cos( pa[,p$variables$TIME] )
      if ("sin.w" %in% p$variables$local_all)  pa$sin.w  = sin( pa[,p$variables$TIME] )
      if ("cos.w2" %in% p$variables$local_all) pa$cos.w2 = cos( 2*pa[,p$variables$TIME] )
      if ("sin.w2" %in% p$variables$local_all) pa$sin.w2 = sin( 2*pa[,p$variables$TIME] )
      if ("cos.w3" %in% p$variables$local_all) pa$cos.w3 = cos( 3*pa[,p$variables$TIME] )
      if ("sin.w3" %in% p$variables$local_all) pa$sin.w3 = sin( 3*pa[,p$variables$TIME] )
      # more than 3 harmonics would not be advisable .. but you would add them here..
      
      if (exists("local_cov", p$variables)) {
        # add time-varying covars .. not necessary except when covars are modelled locally
        pvars2 = names(pa)
        pa$iy = pa$yr - p$yrs[1] + 1 #yr index
        pa$it = p$nw*(pa$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
        for (ci in 1:length(p$variables$local_cov)) {
          vn = p$variables$local_cov[ci]
          pu = NULL
          pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
          nts = ncol(pu)
          if ( nts == p$ny )  {
            pvars2 = c( pvars2, vn )
            pa[,vn] = pu[pa$i, pa$iy ]  
            message("Need to check that data order is correct")
           } else if ( nts == p$nt) {
            pvars2 = c( pvars2, vn )
            pa[,vn] = pu[pa$i, pa$it ]  
            message("Need to check that data order is correct")
          }
        }
      }
    }
    
    # prep dependent data 
    # reconstruct data for modelling (dat) and data for prediction purposes (pa)
    dat = data.frame( Y[YiU] )
    names(dat) = p$variables$Y
    dat[, p$variables$Y] = p$lbm_local_family$linkfun ( dat[, p$variables$Y] ) 
    if (p$lbm_local_modelengine=="habitat") {
      dat[, p$variables$Ylogit ] = p$lbm_local_family_logit$linkfun ( dat[, p$variables$Ylogit] ) ### -- need to conform with data structure ... check once ready
    }
    dat$plon = Yloc[YiU,1]
    dat$plat = Yloc[YiU,2]
    dat$weights = 1 / (( Sloc[Si,1] - dat$plat)**2 + (Sloc[Si,2] - dat$plon)**2 )# weight data in space: inverse distance squared
    dat$weights[ which( dat$weights < 1e-3 ) ] = 1e-3
    dat$weights[ which( dat$weights > 1 ) ] = 1
    
    if (exists("local_cov", p$variables)) {
      for (i in 1:length(p$variables$local_cov )) dat[, p$variables$local_cov[i] ] = Ycov[YiU,i]
    }
     
    if (exists("TIME", p$variables)) {
      dat[, p$variables$TIME ] = Ytime[YiU,] 
      if ( p$variables$TIME != "yr" ) dat$yr = trunc( dat[, p$variables$TIME]) 
      if ("dyear" %in% p$variables$local_all)  dat$dyear = dat[, p$variables$TIME] - dat$yr
      if ("cos.w" %in% p$variables$local_all)  dat$cos.w  = cos( 2*pi*dat[,p$variables$TIME] )
      if ("sin.w" %in% p$variables$local_all)  dat$sin.w  = sin( 2*pi*dat[,p$variables$TIME] )
      if ("cos.w2" %in% p$variables$local_all) dat$cos.w2 = cos( 2*dat[,p$variables$TIME] )
      if ("sin.w2" %in% p$variables$local_all) dat$sin.w2 = sin( 2*dat[,p$variables$TIME] )
      if ("cos.w3" %in% p$variables$local_all) dat$cos.w3 = cos( 3*dat[,p$variables$TIME] )
      if ("sin.w3" %in% p$variables$local_all) dat$sin.w3 = sin( 3*dat[,p$variables$TIME] )
    }

  # use a larger data grid to interpolate.. right now too slow to use so skip this step
    if (p$lbm_local_modelengine %in% c( "spate", "twostep", "lbm_local_modelengine_userdefined", "fft" ) ) {
      # some methods require a uniform prediction grid based upon all dat locations (and time) 
      # begin with "dat"    
      px = dat # only the static parts .. time has to be a uniform grid so reconstruct below
      px$plat = grid.internal( px$plat, p$plats )
      px$plon = grid.internal( px$plon, p$plons )

      # ids = paste(px[,p$variables$LOCS[1] ], px[,p$variables$LOCS[2] ] ) 
      # test which is faster ... remove non-unique?
      ids = array_map( "2->1", round(cbind(px$plon, px$plat)/p$pres+1), c(p$nplons, p$nplats) ) # 100X faster than paste / merge
      todrop = which(duplicated( ids) )
      if (length(todrop>0)) px = px[-todrop,]
      rm(ids, todrop)

      # static vars .. don't need to look up
      tokeep = c(p$variables$LOCS )
      if (exists("weights", dat) ) tokeep = c(tokeep, "weights")
      if (exists("local_cov", p$variables)) {
        for (ci in 1:length(p$variables$local_cov)) {
          vn = p$variables$local_cov[ci]
          pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
          nts = ncol(pu)
          if ( nts==1 ) tokeep = c(tokeep, vn ) 
        }
      }
      px = px[ , tokeep ]
      px_n = nrow(px)
      nts = vn = pvars2 = NULL

      # add temporal grid
      if ( exists("TIME", p$variables) ) {
        px = cbind( px[ rep.int(1:px_n, p$nt), ], 
                        rep.int(p$ts, rep(px_n, p$nt )) )
        names(px)[ ncol(px) ] = p$variables$TIME 
        if ( p$variables$TIME != "yr" ) px$yr = trunc( px[,p$variables$TIME] )
        if ("dyear" %in% p$variables$local_all)  px$dyear = px[, p$variables$TIME] - px$yr  # fractional year
        if ("cos.w" %in% p$variables$local_all)  px$cos.w  = cos( px[,p$variables$TIME] )
        if ("sin.w" %in% p$variables$local_all)  px$sin.w  = sin( px[,p$variables$TIME] )
        if ("cos.w2" %in% p$variables$local_all) px$cos.w2 = cos( 2*px[,p$variables$TIME] )
        if ("sin.w2" %in% p$variables$local_all) px$sin.w2 = sin( 2*px[,p$variables$TIME] )
        if ("cos.w3" %in% p$variables$local_all) px$cos.w3 = cos( 3*px[,p$variables$TIME] )
        if ("sin.w3" %in% p$variables$local_all) px$sin.w3 = sin( 3*px[,p$variables$TIME] )
        # more than 3 harmonics would not be advisable .. but you would add them here..
      }

      if (exists("local_cov", p$variables)) {
        # add time-varying covars .. not necessary except when covars are modelled locally
        pvars2 = names(px)
        px$iy = px$yr - p$yrs[1] + 1 #yr index
        px$it = p$nw*(px$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
        for (ci in 1:length(p$variables$local_cov)) {
          vn = p$variables$local_cov[ci]
          pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
          nts = ncol(pu)
          if ( nts== 1) {
            # static vars are retained in the previous step
          } else if ( nts == p$ny )  {
            pvars2 = c( pvars2, vn )
            px[,vn] = pu[px$i, px$iy ]  
           } else if ( nts == p$nt) {
            pvars2 = c( pvars2, vn )
            px[,vn] = pu[px$i, px$it ]  
          }
        } # end for loop
        nts = vn = pvars2 = NULL
      } # end if
    }

    nu = phi = varSpatial = varObs = NULL
    if ( exists("nu", ores) ) nu = ores$nu
    if ( exists("phi", ores) && ores$phi > (p$pres/2)) phi = ores$phi 
    if ( exists("varSpatial", ores) ) varSpatial = ores$varSpatial
    if ( exists("varObs", ores) ) varObs = ores$varObs
  
    if (is.null(nu)) nu = p$lbm_lowpass_nu
    if (is.null(phi)) phi = lbm_distance_cur/sqrt( 8*nu) # crude estimate of phi based upon current scaling  distance approximates the range at 90% autocorrelation(e.g., see Lindgren et al. 2011)
    if (is.null(varSpatial)) varSpatial =0.5 * var(dat[, p$variables$Y], na.rm=TRUE)
    if (is.null(varObs)) varObs = varSpatial
    
    # model and prediction 
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
      glm = lbm__glm( p, dat, pa, nu=nu, phi=phi ), 
      krige = lbm__krige( p, dat, pa, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial ), # TODO
      LaplacesDemon = lbm__LaplacesDemon( p, dat, pa ),
      splancs = lbm__spate( p, dat, pa ), # TODO
      spate = lbm__spate( p, dat, pa, sloc=Sloc[Si,], px=px ), 
      fft = lbm__fft( p, dat, pa, nu=nu, phi=phi ), 
      twostep = lbm__twostep( p, dat, pa, px=px, nu=nu, phi=phi  ), # slow ...!
      lbm_local_modelengine_userdefined = p$lbm_local_modelengine_userdefined( p, dat, pa)
    ) )


    if (0) {
      lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
       
      lattice::levelplot( mean ~ plon + plat, data=res$predictions, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
   
      for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
    }

    rm(dat); gc()
    if ( inherits(res, "try-error") ) {
      dat = pa = px = NULL
      next()
    }
    
    if ( is.null(res)) {
      dat = pa = px = NULL
      next()
    }
    if ( all( !is.finite(res$predictions$mean ))) {
      dat = pa = px = res = NULL
      next()
    }


    res$predictions$mean = p$lbm_local_family$linkinv( res$predictions$mean )
    # res$predictions$sd   = p$lbm_local_family$linkinv( res$predictions$sd )
    if (p$lbm_local_modelengine=="habitat") {
      res$predictions$logitmean = p$lbm_local_family_logit$linkinv( res$predictions$logitmean )
      # res$predictions$logitsd   = p$lbm_local_family_logit$linkinv( res$predictions$logitsd )
    }
 
    if (exists( "lbm_quantile_bounds", p)) {
      tq = quantile( Y[YiU], probs=p$lbm_quantile_bounds, na.rm=TRUE  )
      toolow  = which( res$predictions$mean < tq[1] )
      toohigh = which( res$predictions$mean > tq[2] )
      if (length( toolow) > 0)  res$predictions$mean[ toolow] = tq[1]
      if (length( toohigh) > 0) res$predictions$mean[ toohigh] = tq[2]
    }
    
    ii = which( is.finite(res$predictions$mean ))
    if (length(ii) < 5) {
      dat = pa = px = res = NULL
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

        if (p$lbm_local_modelengine=="habitat") {
          logit_stdev_update =  Plogitsd[ui] + ( res$predictions$logitsd[u] -  Plogitsd[ui] ) / Pn[ui]
          logit_means_update = ( Plogit[ui] / Plogitsd[ui]^2 + res$predictions$logitmean[u] / res$predictions$logitsd[u]^2 ) / ( Plogitsd[ui]^(-2) + res$predictions$logitsd[u]^(-2) )
          mm = which(is.finite( logit_means_update + logit_stdev_update ))
          if( length(mm)> 0) {
            iumm = ui[mm]
            Plogitsd[iumm] = logit_stdev_update[mm]
            Plogit  [iumm] = logit_means_update[mm]
          }
          logit_stdev_update = NULL
          logit_means_update = NULL
        }
        rm(ui, mm, iumm)
      }

      # first time # no data yet
      v = setdiff(1:npred, u)         
      if ( length(v) > 0 ) {
        vi = res$predictions$i[v]
        Pn [vi] = 1
        P  [vi] = res$predictions$mean[v]
        Psd[vi] = res$predictions$sd[v]
        if (p$lbm_local_modelengine=="habitat") {
          Plogit  [vi] = res$predictions$logitmean[v]
          Plogitsd[vi] = res$predictions$logitsd[v]
        }
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
        if (p$lbm_local_modelengine=="habitat") {
          logit_stdev_update =  Plogitsd[ui,] + ( res$predictions$logitsd[u] -  Plogitsd[ui,] ) / Pn[ui]
          logit_means_update = ( Plogit[ui,] / Plogitsd[ui,]^2 + res$predictions$logitmean[u] / res$predictions$logitsd[u]^2 ) / ( Plogitsd[ui,]^(-2) + res$predictions$logitsd[u]^(-2) )
          updates = logit_means_update + logit_stdev_update
          if (!is.matrix(updates)) next()
          mm = which( is.finite( rowSums(updates)))  # created when preds go outside quantile bounds .. this removes 
          if( length(mm)> 0) {
            iumm = ui[mm]
            Plogitsd[iumm,] = logit_stdev_update[mm,]
            Plogit  [iumm,] = logit_means_update[mm,]
            iumm = NULL
          } 
          logit_stdev_update = NULL
          logit_means_update = NULL
        }
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
        if (p$lbm_local_modelengine=="habitat") {
          Plogit  [vi,] = res$predictions$logitmean[v]
          Plogitsd[vi,] = res$predictions$logitsd[v]
        }
        rm(vi)
      } 
    }


    # save stats
    for ( k in 1: length(p$statsvars) ) {
      if (exists( p$statsvars[k], res$lbm_stats )) {
        S[Si,k] = res$lbm_stats[[ p$statsvars[k] ]]
      }
    }
    
      if (0) {
     
        lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
        
        for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      
        for (i in 1:p$nt) {
          print( lattice::levelplot( P[pa$i,i] ~ Ploc[pa$i,1] + Ploc[ pa$i, 2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
        }

        v = res$predictions
        if ( exists("TIME", p$variables) ){
          v = v[which( v[,p$variables$TIME]==1990.55),]
        }
        require(lattice)
        levelplot( mean ~ plon+plat, v, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
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

