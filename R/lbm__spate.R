
lbm__spate = function( p, dat, pa, sloc, distance, nu, phi, varObs, varSpatial ) {

  # require(spate) #\\ SPDE solution via FFT using the spate library

  sdTotal=sd(dat[,p$variable$Y], na.rm=T)
  ndata = nrow(dat)
  
  datgridded = dat # only the static parts .. time has to be a uniform grid so reconstruct below

  ids = array_map( "xy->1", datgridded[, c("plon", "plat")], gridparams=p$gridparams ) # 100X faster than paste / merge
  todrop = which(duplicated( ids) )
  if (length(todrop>0)) datgridded = datgridded[-todrop,]
  rm(ids, todrop)

  if (p$lbm_spate_boost_timeseries ) {
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

    datgridded = datgridded[ , tokeep ]
    datgridded_n = nrow(datgridded)
    nts = vn = NULL

    # add temporal grid
    if ( exists("TIME", p$variables) ) {
      datgridded = cbind( datgridded[ rep.int(1:datgridded_n, p$nt), ], 
                      rep.int(p$prediction.ts, rep(datgridded_n, p$nt )) )
      names(datgridded)[ ncol(datgridded) ] = p$variables$TIME 
      datgridded = cbind( datgridded, lbm_timecovars ( vars=p$variables$local_all, ti=datgridded[,p$variables$TIME]  ) )
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
          datgridded$iy = datgridded$yr - p$yrs[1] + 1 #yr index
          datgridded[,vn] = pu[ cbind(datgridded$i, datgridded$iy) ]  
         } else if ( nts == p$nt) {
          datgridded$it = p$nw*(datgridded$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
          datgridded[,vn] = pu[ cbind(datgridded$i, datgridded$it) ]  
        }
      } # end for loop
      nts = vn = NULL
    } # end if

    ts_gam = lbm__gam( p, dat, datgridded ) # currently only a GAM is enabled for the TS component

    if (is.null( ts_gam)) return(NULL)
    if (ts_gam$lbm_stats$rsquared < p$lbm_rsquared_threshold ) return(NULL)

    # range checks
    rY = range( dat[,p$variables$Y], na.rm=TRUE)
    toosmall = which( ts_gam$predictions$mean < rY[1] )
    toolarge = which( ts_gam$predictions$mean > rY[2] )
    if (length(toosmall) > 0) ts_gam$predictions$mean[toosmall] =  NA 
    if (length(toolarge) > 0) ts_gam$predictions$mean[toolarge] =  NA
   
    # overwrite dat with interpolated predictions
    datgridded = ts_gam$predictions
    ts_gam = NULL
    
    # revert to response scale as the following expects this:
    datgridded$mean = p$lbm_local_family$linkinv( datgridded$mean )
    datgridded$sd = p$lbm_local_family$linkinv( datgridded$sd )

    names(datgridded)[which(names(datgridded)=="mean")] = p$variables$Y
    names(datgridded)[which(names(datgridded)=="sd")] = paste(p$variables$Y, "sd", sep=".")


    gc()
  }

  windowsize.half = floor(distance/p$pres)
  pa_w = -windowsize.half : windowsize.half # default window size 
  # pa_w = pa_w[ -length(pa_w)] # must be even 
  pa_w_n = length(pa_w)
  adims = c(p$nt, pa_w_n, pa_w_n ) 
  datgridded$id = array_map( "3->1", 
    coords = round( cbind( 
      ( (datgridded[,p$variables$TIME ] - p$prediction.ts[1] ) / p$tres) + 1 ,
      ( windowsize.half + (datgridded[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
      ( windowsize.half + (datgridded[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1)), 
    dims=adims )
  o = which( datgridded$id > prod(adims)  | datgridded$id <= 0 ) # remove the area outside the aoi
  if (length(o) > 0 ) datgridded = datgridded[-o,]
  
  ddup = which(duplicated(datgridded$id))
  
  if ( length( ddup) > 0 ) {
    dups = unique(datgridded$id[ddup])
    for ( i in dups ) {
      j = which( datgridded$id== i)
      meanvalue = mean( datgridded[j, p$variable$Y], na.rm=TRUE)  
      datgridded[j, p$variable$Y] = NA
      datgridded[j[1], p$variable$Y] = meanvalue
    }
    datgridded = datgridded[ which(is.finite(datgridded[, p$variable$Y]) ) , ]
  } 

  xM = array( NA, dim=adims )
  xM[datgridded$id] = datgridded[,p$variables$Y]


  if (0) {
  # ignore this for now ..
  # prediction covariates i.e., independent variables/ covariates
    pvars = c("plon", "plat", "i")
    if (p$nloccov > 0) {
      # .. not necessary except when covars are modelled locally
      for (ci in 1:p$nloccov) {
        vn = p$variables$local_cov[ci]
        pu = NULL
        pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
        nts = ncol(pu)
        if ( nts== 1 ) {
          pvars = c( pvars, vn )
          datgridded[,vn] = pu[datgridded$i]  # ie. a static variable
        }
      }
    }
    datgridded = datgridded[, pvars]
    datgridded = cbind( datgridded, lbm_timecovars ( vars=p$variables$local_all, ti=datgridded[,p$variables$TIME]  ) )

    if (p$nloccov > 0) {
      # add time-varying covars .. not necessary except when covars are modelled locally
      for (ci in 1:p$nloccov) {
        vn = p$variables$local_cov[ci]
        pu = NULL
        pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
        nts = ncol(pu)
        if ( nts == p$ny )  {
          datgridded$iy = datgridded$yr - p$yrs[1] + 1 #yr index
          datgridded[,vn] = pu[ cbind(datgridded$i, datgridded$iy) ]  
          message("Need to check that datgriddeda order is correct")
        } else if ( nts == p$nt ) {
          datgridded$it = p$nw*(datgridded$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
          datgridded[,vn] = pu[ cbind(datgridded$i, datgridded$it) ]  
          message("Need to check that data order is correct")
        } else if (nts==1) { } #nothing to do .. already processed above }
      }
    }
  }

  # shorten by 1 row and column to make it even
  nsq = pa_w_n-1
  w = matrix( xM[,1:nsq, 1:nsq ], nrow=p$nt )

  SV = c(rho0 = 0.25, 
        sigma2 = max(0.1, varSpatial, na.rm=TRUE), zeta = 0.25, rho1 = 0.25, gamma = 1, alpha = 0.1, 
        muX = 0, muY = 0, tau2 = max(0.005, varObs, na.rm=TRUE) )

  g = try( spate.mcmc( y=w, n=nsq, Padding=FALSE, trace=FALSE, Nmc=5000, SV=SV,
    adaptive=TRUE, Separable=FALSE, Drift=TRUE, Diffusion=TRUE, plotTrace=FALSE, BurnIn=1500,
    BurnInCovEst=1000, NCovEst=1000 ), silent=TRUE)

  if ("try-error" %in% class(g)) {
    g = try( spate.mcmc( y=w, n=nsq, Padding=FALSE, trace=FALSE, Nmc=5000, 
      adaptive=TRUE, Separable=FALSE, Drift=TRUE, Diffusion=TRUE, plotTrace=FALSE , BurnIn=2000 ), silent=TRUE) # longer burn-in (1000 is default) and alternate rnd seed
  }

  if ("try-error" %in% class(g)) {
    g = try( spate.mcmc( y=w, n=nsq, Padding=FALSE, trace=FALSE, seed=1001, Nmc=5000,
      adaptive=TRUE, Separable=FALSE, Drift=TRUE, Diffusion=TRUE, plotTrace=FALSE , BurnIn=2000 ), silent=TRUE) # longer burn-in (1000 is default) and alternate rnd seed
  }

  if ("try-error" %in% class(g)) return(NULL)
 
  spp <- spate.predict(y=w, tPred=(1:p$nt), 
    spateMCMC=g, Nsim=1010, BurnIn=10, DataModel="Normal", trace=FALSE ) # nu=nu, defulat is to assume nu =1 
  #  DimRed=TRUE, NFour=101

# must be same order as p$statsvars
  pmean_vars = c( "rho_0", "zeta", "rho_1", "gamma", "alpha", "mu_x", "mu_y", "sigma^2", "tau^2" ) #reorder
  psd_vars   = c( "rho_0.sd", "zeta.sd", "rho_1.sd", "gamma.sd", "alpha.sd", "mu_x.sd", "mu_y.sd" ) # drop sd of variance terms

  pmean =  apply(g$Post, 1, mean)[pmean_vars] 
  psd =  apply(g$Post, 1, sd)[pmean_vars] 

  names(pmean) = pmean_vars
  names(psd) = paste( pmean_vars, "sd", sep=".")
  psd = psd[psd_vars]

  rm(g, w); gc()

  # determine prediction locations and time slices
  iwplon = round( (sloc[1]-p$origin[1])/p$pres + 1 + pa_w )
  iwplat = round( (sloc[2]-p$origin[2])/p$pres + 1 + pa_w )

  oo = NULL
  oo = data.frame( iplon = rep.int(iwplon, pa_w_n) , 
                   iplat = rep.int(iwplat, rep.int(pa_w_n, pa_w_n)) )

  Ploc = lbm_attach( p$storage.backend, p$ptr$Ploc )
  ploc_ids = array_map( "xy->1", Ploc[], gridparams=p$gridparams )

  oo$i = match( array_map( "2->1", oo[, c("iplon", "iplat")], gridparams=p$gridparams ), ploc_ids )
  oo = cbind( oo[ rep.int(1:nrow(oo), p$nt), ], 
                  rep.int(p$prediction.ts, rep(nrow(oo), p$nt )) )
  names(oo)[4] = p$variables$TIME 

  bad = which( (oo$iplon < 1 & oo$iplon > p$nplons) | (oo$iplat < 1 & oo$iplat > p$nplats) )
  if (length(bad) > 0 ) oo = oo[-bad,]
  if (nrow(oo)< 5) return(NULL)
  
  pa$id = array_map( "3->1", 
    coords = round( cbind( 
      ( (pa[,p$variables$TIME ] - p$prediction.ts[1] ) / p$tres) + 1 ,
      ( windowsize.half + (pa[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
      ( windowsize.half + (pa[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1)), 
    dims=adims )

  # means 
  xM[,1:nsq,1:nsq] = apply(spp, c(1,2), mean)
  pa$mean = NA
  pa$mean = xM[pa$id ]

  datgridded$mean = NA
  datgridded$mean = xM[datgridded$id]

  # sd
  xM[,1:nsq,1:nsq] = apply(spp, c(1,2), sd)
  pa$sd = NA
  pa$sd = xM[pa$id]


  if (0) {
    
    Pmean = apply(spp, c(1,2), mean)
    zlim=range(Pmean)
    x11()
    for( i in 1:p$nt) {  
      # i = 1
      pause(0.2)
  #    image(1:nsq,1:nsq, matrix(Pmean[i,], nrow=nsq), zlim=zlim,
  #            main=paste("Mean predicted field at t=",i,sep=""), xlab="",ylab="",col=cols())
      oo = expand.grid( 1:nsq, 1:nsq)
      oo = as.data.frame( cbind(oo, Pmean[i,] ) )
      oom = MBA::mba.surf( oo, nsq, nsq, extend=TRUE)$xyz.est
      surface(oom, zlim=zlim, type="I") # I = no countour, C=with couontours
    }

    Psd = apply(spp, c(1,2), sd)
    zlim=range(Psd)
    for( i in 1:p$nt) {  
    # i = 1
      pause(0.25)
      # image(1:nsq,1:nsq,matrix(Psd[i,],nrow=nsq), zlim=zlim,
      #       main=paste("Sd of predicted field at t=",i,sep=""), xlab="",ylab="",col=cols())
      oo = expand.grid( 1:nsq, 1:nsq)
      oo = as.data.frame( cbind(oo, Psd[i,] ) )
      oom = MBA::mba.surf( oo, 300, 300, extend=TRUE)$xyz.est
      surface(oom, zlim=zlim, contour=FALSE)
    } 

    # just raw data running average
    x11()
    zlim=range(dat[, p$variables$Y ])
    ww = 3
    for( i in 1:p$nt) {  
      pause(0.2)
      print(p$prediction.ts[i] )
      ii = i+c(-floor(p$nw/3):floor(p$nw/3))
      ii = ii[ ii >= 1 & ii <= p$nt ]
      trange = range( p$prediction.ts[ii] )
      di = which( dat[ , p$variables$TIME ] > trange[1] & dat[ , p$variables$TIME ] < trange[2]  )
      if( length(di) < 10) next() 
      dmbas = MBA::mba.surf( dat[di, c( p$variables$LOCS, p$variables$Y) ], 300, 300, extend=TRUE)$xyz.est
      surface(dmbas, zlim=zlim)
    }
   
     plot( plat~plon, dat[ dat$yr==2016 ,], pch=20 )
     text( plat~plon, labels=dat[dat$yr==2016,"t"] , dat[ dat$yr==2016 ,] )

        # debugging plots
        # boosted "data"
        i = 658
        zlim=range(dat[, p$variables$Y ])
 
        x11()
        for (i in 1:p$nt) {
         pause(0.2)
        xi = which( datgridded[ , p$variables$TIME ] == p$prediction.ts[i] )
        mbas = MBA::mba.surf( datgridded[xi, c( p$variables$LOCS, p$variables$Y) ], 300, 300, extend=TRUE)$xyz.est
        surface(mbas, zlim=zlim)
        }

        #data
        x11()
        for (i in 1:p$nt) {
          pause(0.2)
          di = which( floor(dat[ , p$variables$TIME ]) == floor(p$prediction.ts[i]) )
          dmbas = MBA::mba.surf( dat[di, c( p$variables$LOCS, p$variables$Y) ], 300, 300, extend=TRUE)$xyz.est
          surface(dmbas, zlim=zlim)
        }

        # spate
        x11()
        for (i in 1:p$nt) {
          pause(0.2)
          oo = expand.grid( 1:nsq, 1:nsq)
          oo = as.data.frame( cbind(oo, Pmean[i,] ) )
          oom = MBA::mba.surf( oo, 300, 300, extend=TRUE)$xyz.est
          surface(oom, zlim=zlim)
        }

        zlim=range(dat[, p$variables$Y ])
        x11()
        par( mfrow=c(1,3))
        dmbas0 = MBA::mba.surf( dat[, c( p$variables$LOCS, p$variables$Y) ], 300, 300, extend=TRUE)$xyz.est
        for (i in 1:p$nt) {
           pause(1)
          # data
          di = which( floor(dat[ , p$variables$TIME ]) == floor(p$prediction.ts[i]) )
          if (length(di) > 5) {
            dmbas = MBA::mba.surf( dat[di, c( p$variables$LOCS, p$variables$Y) ], 300, 300, extend=TRUE)$xyz.est
           } else { 
            dmbas = dmbas0 
          }
          surface(dmbas, zlim=zlim)

          # grid boosted
          xi = which( datgridded[ , p$variables$TIME ] == p$prediction.ts[i] )
          mbas = MBA::mba.surf( na.omit(datgridded[xi, c( p$variables$LOCS, p$variables$Y) ]), 300, 300, extend=TRUE)$xyz.est
          surface(mbas, zlim=zlim)

         # spate
          oo = expand.grid( 1:nsq, 1:nsq)
          oo = as.data.frame( cbind(oo, Pmean[i,] ) )
          oom = MBA::mba.surf( oo, 300, 300, extend=TRUE)$xyz.est
          surface(oom, zlim=zlim)
      }  # end for  

  } # end debug

  rm(oo, xM, spp); gc()

  # plot(mean ~ z , datgridded)

  ss = lm( datgridded$mean ~ datgridded[,p$variables$Y], na.action=na.omit )
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rm(datgridded); gc()

  rsquared = summary(ss)$r.squared
  if (rsquared < p$lbm_rsquared_threshold ) return(NULL)


  lbm_stats = c(list( sdTotal=sdTotal, rsquared=rsquared, ndata=ndata), pmean, psd)  
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, lbm_stats=lbm_stats ) )  
}

