

lstfilter_interpolate_xy_simple = function( interp.method, data, locsout, datagrid=NULL,
  trimquants=TRUE, trimprobs=c(0.025, 0.975), 
  nr=NULL, nc=NULL, phi=1, xwidth=phi*10, ywidth=phi*10, nu=0.5 ) {
  #\\ reshape after interpolating to fit the output resolution 
  #\\ designed for interpolating statistics  ...

  # trim quaniles in case of extreme values
  if (trimquants) {
    vq = quantile( data$z, probs=trimprobs, na.rm=TRUE )
    ii = which( data$z < vq[1] )
    if ( length(ii)>0) data$z[ii] = vq[1]
    jj = which( data$z > vq[2] )
    if ( length(jj)>0) data$z[jj] = vq[2]
  }

  out = NULL

  # ------

  if (interp.method == "multilevel.b.splines") {
    library(MBA)
    out = mba.surf(data, no.X=nr, no.Y=nc, extend=TRUE)
    if (0) {
      image(out, xaxs = "r", yaxs = "r", main="Observed response")
      locs= cbind(data$x, data$y)
      points(locs)
      contour(out, add=T)
    }
    return(out$xyz.est)
  }

  # ------

  if (interp.method == "fft") {
    require(fields)

    rY = quantile( x[,p$variables$Y], probs=p$lstfilter_quantile_bounds, na.rm=TRUE)

    im = as.image(data, x=locsout, nx=nr, ny=nc, na.rm=TRUE)
    grid <- list(x = im$x, y = im$y)
    dx <- grid$x[2] - grid$x[1]
    dy <- grid$y[2] - grid$y[1]
    nr2 = 2 * nr
    nc2 = 2 * nc

    dgrid = make.surface.grid(list((1:nr2) * dx, (1:nc2) * dy))
    center = matrix(c((dx * nr2)/2, (dy * nc2)/2), nrow = 1, 
        ncol = 2)

    mC = matrix(0, nrow = nr2, ncol = nc2)
    mC[nr, nc] = 1

    
    # first pass with the global params to get closest fit to data 
    AC_global = stationary.cov( dgrid, center, Covariance="Matern", range=p$lstfilter_phi, nu=p$lstfilter_nu )
    mAC_global = as.surface(dgrid, c(AC_global))$z
    fW_global = fft(mAC_global)/(fft(mC) * nr2 * nc2)

    # second pass with local fits to data to smooth what can be smoothed
    AC_local  = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
    mAC_local = as.surface(dgrid, c(AC_local))$z
    fW_local = fft(mAC_local)/(fft(mC) * nr2 * nc2)

    rm(dgrid, AC_global, AC_local, mC, mAC_local, mAC_global); gc()

    
    # not finished, see lstfilter__kerneldensity
    x_id = cbind( (x[xi,p$variables$LOCS[1]]-x_r[1])/p$pres + 1, 
                  (x[xi,p$variables$LOCS[2]]-x_c[1])/p$pres + 1 )
        
    # counts
    mN = matrix(0, nrow = nr2, ncol = nc2)
    mN[x_id] = tapply( rep(1, length(xi)), INDEX=x_id, FUN=sum, na.rm=TRUE )
    mN[!is.finite(mN)] = 0
    
    # density
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[x_id] = x[xi,p$variables$Y] # fill with data in correct locations
    mY[!is.finite(mY)] = 0
    
    # estimates based upon a global nu,phi .. they will fit to the immediate area near data and so retain their structure
    fN = Re(fft(fft(mN) * fW_global, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * fW_global, inverse = TRUE))[1:nr,1:nc]
    Z = fY/fN
    iZ = which( !is.finite( Z))
    if (length(iZ) > 0) Z[iZ] = NA
    lb = which( Z < rY[1] )
    if (length(lb) > 0) Z[lb] = NA
    ub = which( Z > rY[2] )
    if (length(ub) > 0) Z[ub] = NA
    # image(Z)

    # estimates based upon local nu, phi .. this will over-smooth so if comes as a second pass 
    # to fill in areas with no data (e.g., far away from data locations)
    fN = Re(fft(fft(mN) * fW_local, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * fW_local, inverse = TRUE))[1:nr,1:nc]
    Z_local = fY/fN
    iZ = which( !is.finite( Z_local))
    if (length(iZ) > 0) Z_local[iZ] = NA
    lb = which( Z_local < rY[1] )
    if (length(lb) > 0) Z_local[lb] = NA
    ub = which( Z_local > rY[2] )
    if (length(ub) > 0) Z_local[ub] = NA

    toreplace = which(!is.finite(Z)) 
    if (length(toreplace) > 0 )  Z[toreplace] = Z_local[toreplace]

    # image(Z)
    return (Z)
  }

  # ------

  if (interp.method == "kernel.density") {
    # default :: create a "surface" and reshape to a grid using (gaussian) kernel-based smooth via FFT
    require(fields)

    if(0) {
      data = RMelevation
      image( data )
      locsout = expand.grid( data$x, data$y)
      nr = length( data$x)
      nc = length(data$y)
      phi = 1
      xwidth = 2
      ywidth = 2 
    }

    isurf = fields::interp.surface( data, loc=locsout  )
    # lattice::levelplot( isurf ~ locsout[,1] + locsout[,2], aspect="iso",  col=topo.colors(256) )
    zout = matrix( isurf, nrow=nr, ncol=nc )
    # image(zout)
    out = fields::image.smooth( zout, theta=phi, xwidth=xwidth, ywidth=ywidth ) 
    image(out)
    return (out$z)
  }

  # ------------------------

  if (interp.method == "fft") {
    # default :: create a "surface" and reshape to a grid using (gaussian) kernel-based smooth via FFT
    require(fields)

    if(0) {
      data = RMelevation
      image( data )
      datagrid = data[c("x", "y")] 
      locsout = expand.grid( x=datagrid$x, y=datagrid$y)
      x = locsout[,1]
      y = locsout[,2]
      z = c(data$z)
      nr = length( datagrid$x)
      nc = length( datagrid$y)
      dx = min(diff(datagrid$x))
      dy = min(diff(datagrid$y))
      nu = 0.5
      phi = min(dx, dy) / 10

      o = lstfilter::lstfilter_variogram( xy=cbind(x,y), z=z, methods="gstat" ) 
      # suggest: nu=0.3; phi =1.723903

      keep = sample.int( nrow(locsout), 1000 ) 
      x = x[keep]
      y = y[keep]
      z = z[keep]
      
      o = lstfilter::lstfilter_variogram( xy=cbind(x,y), z=z, methods="gstat" ) 
      o = lstfilter::lstfilter_variogram( xy=cbind(x,y), z=z, methods="fast" ) 
 
    }
    qz = range(z, na.rm=TRUE)

    nr2 = 2 * nr
    nc2 = 2 * nc
    
    Z2P = as.matrix( cbind( (x-datagrid$x[1])/dx + 1 , (y-datagrid$y[1] )/dy+1) ) # row, col indices in matrix form Z
    zp = lstfilter::array_map( "2->1", Z2P, c(nr2, nc2) ) # map to larger grid

    dgrid = make.surface.grid(list((1:nr2) * dx, (1:nc2) * dy))
    center = matrix(c((dx * nr2)/2, (dy * nc2)/2), nrow = 1, ncol = 2)
    AC = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
    mAC = matrix(c(AC), nrow = nr2, ncol = nc2) # or .. mAC = as.surface(dgrid, c(AC))$z
    mC = matrix(0, nrow = nr2, ncol = nc2)
    mC[nr, nc] = 1
    fW = fft(mAC)/(fft(mC) * nr2 * nc2)
    rm(dgrid, AC, mAC, mC); gc()

    # counts
    mW = matrix(0, nrow = nr2, ncol = nc2)
    mW[zp] = tapply( ifelse( is.finite(z), 1, 0 ), INDEX=zp, FUN=sum, na.rm=TRUE )
    mY[!is.finite(mW)] = 0
    fN = Re(fft(fft(mW) * fW, inverse = TRUE))[1:nr,1:nc]

    # data
    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[zp] = tapply( z, INDEX=zp, FUN=mean, na.rm=TRUE ) # fill with data in correct locations
    mY[!is.finite(mY)] = 0
    fY = Re(fft(fft(mY) * fW, inverse = TRUE))[1:nr,1:nc]
    Z = fY/fN
    Z[ Z>qz[2] ]=NA
    Z[ Z<qz[1] ]=NA
    x11(); image(Z)


    return (out)
  }

  # ------

  if (interp.method == "inla.spde" ) {
    require(INLA)
    FM = formula( "Y ~ -1 + intercept + f( spatial.field, model=SPDE )" )
    Y = data$z
    locs= cbind(data$x, data$y)
    rm (data); gc()
    method="fast"  # "direct" is slower and more accurate
    nsamples=5000 
  # identity links by default .. add more if needed here
    locs = as.matrix( locs)
    lengthscale = max( diff(range( locs[,1])), diff(range( locs[,2]) )) / 10  # in absence of range estimate take 1/10 of domain size
    ndata = length(Y)
    noise = lengthscale * 1e-9
    locs = locs + runif( ndata*2, min=-noise, max=noise ) # add  noise  to prevent a race condition .. inla does not like uniform grids
    MESH = lstfilter_mesh_inla( locs, lengthscale=lengthscale )
    if ( is.null( MESH) ) return( "Mesh Error" )
    SPDE = inla.spde2.matern( MESH,  alpha=2 ) # alpha is 2*nu (Bessel smoothness factor)
    varY = as.character( FM[2] )
    obs_index = inla.spde.make.index(name="spatial.field", SPDE$n.spde)
    obs_eff = list()
    obs_eff[["spde"]] = c( obs_index, list(intercept=1) )
    obs_A = list( inla.spde.make.A( mesh=MESH, loc=locs[,] ) ) # no effects
    obs_ydata = list()
    obs_ydata[[ varY ]] =  Y 
    DATA = inla.stack( tag="obs", data=obs_ydata, A=obs_A, effects=obs_eff, remove.unused=FALSE )
    if ( method=="direct") {
      # direct method
      preds_index = inla.spde.make.index( name="spatial.field", SPDE$n.spde)
      preds_eff = list()
      preds_eff[["spde"]] = c( list( intercept=rep(1,MESH$n )),
           inla.spde.make.index(name="spatial.field", n.spde=SPDE$n.spde) )
      ydata = list()
      ydata[[ varY ]] = NA
      Apreds = inla.spde.make.A(MESH, loc=as.matrix( locsout ) )
      PREDS = inla.stack( tag="preds", data=list( Y=NA), A=list(Apreds),
        effects=list( c( list(intercept=rep(1, MESH$n)), inla.spde.make.index( name="spatial.field", MESH$n))) )
      DATA = inla.stack(DATA, PREDS)
      i_data = inla.stack.index( DATA, "preds")$data
    }
    RES = NULL
    RES = lstfilter_inla_call( FM=FM, DATA=DATA, SPDE=SPDE, FAMILY="gaussian" )
    # extract summary statistics from a spatial (SPDE) analysis and update the output file
    # inla.summary = lstfilter_summary_inla_spde2 = ( RES, SPDE )
    # inla.spde2.matern creates files to disk that are not cleaned up:
    spdetmpfn = SPDE$f$spde2.prefix
    fns = list.files( dirname( spdetmpfn ), all.files=TRUE, full.names=TRUE, recursive=TRUE, include.dirs=TRUE )
    oo = grep( basename(spdetmpfn), fns )
    if(length(oo)>0) file.remove( sort(fns[oo], decreasing=TRUE) )
    rm( SPDE, DATA ); gc()
    # predict upon grid
    if ( method=="direct" ) {
      # direct method ... way too slow to use for production runs
      preds = as.data.frame( locsout )
      preds$xmean =  RES$summary.fitted.values[ i_data, "mean"] 
      preds$xsd   =  RES$summary.fitted.values[ i_data, "sd"] 
      rm(RES, MESH); gc()
    }
    if (method=="fast") {
      posterior.extract = function(s, rnm) {
        # rnm are the rownames that will contain info about the indices ..
        # optimally the grep search should only be done once but doing so would
        # make it difficult to implement in a simple structure/manner ...
        # the overhead is minimal relative to the speed of modelling and posterior sampling
        i_intercept = grep("intercept", rnm, fixed=TRUE ) # matching the model index "intercept" above .. etc
        i_spatial.field = grep("spatial.field", rnm, fixed=TRUE )
        return(  s$latent[i_intercept,1] + s$latent[ i_spatial.field,1] )
      }
      # note: locsout seems to be treated as token locations and really its range and dims controls output
      pG = inla.mesh.projector( MESH, loc=as.matrix( locsout ) )
      posterior.samples = inla.posterior.sample(n=nsamples, RES)
      rm(RES, MESH); gc()
      rnm = rownames(posterior.samples[[1]]$latent )
      posterior = sapply( posterior.samples, posterior.extract, rnm=rnm )
      posterior =  posterior    # return to original scale
      rm(posterior.samples); gc()
          # robustify the predictions by trimming extreme values .. will have minimal effect upon mean
          # but variance estimates should be useful/more stable as the tails are sometimes quite long
          for (ii in 1:nrow(posterior )) {
            qnt = quantile( posterior[ii,], probs=c(0.025, 0.975), na.rm=TRUE )
            toolow = which( posterior[ii,] < qnt[1] )
            toohigh = which (posterior[ii,] > qnt[2] )
            if (length( toolow) > 0 ) posterior[ii,toolow] = qnt[1]
            if (length( toohigh) > 0 ) posterior[ii,toohigh] = qnt[2]
          }
      # posterior projection is imperfect right now .. not matching the actual requested locations
      preds = data.frame( plon=pG$loc[,1], plat = pG$loc[,2])
      preds$xmean = c( inla.mesh.project( pG, field=apply( posterior, 1, mean, na.rm=TRUE )  ))
      preds$xsd   = c( inla.mesh.project( pG, field=apply( posterior, 1, sd, na.rm=TRUE )  ))
      rm (pG)
    }
    if (0) {
      require(lattice)
      levelplot( log( xmean)  ~ plon+plat, preds, aspect="iso",
                labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      levelplot( log (xsd )  ~ plon+plat, preds, aspect="iso",
                labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }
    return(preds$xmean) 
  }

  if (interp.method=="spBayes"){
    #// interpolate via spBayes ~ random field approx via matern
    #// method determines starting parameter est method to seed spBayes
    #// NOTE:: TOO SLOW TO USE FOR OPERATIONAL WORK
    require(spBayes)
    library(MBA)
    require( geoR )
    require( stat )
    require( sp )
    require( rgdal )

    xy = as.data.frame( cbind( x,y) )
    z = z / sd( z, na.rm=TRUE)

    pCovars = as.matrix( rep(1, nrow(locsout)))  # in the simplest model, 1 col matrix for the intercept

    stv = lstfilter_variogram( xy, z, methods=method )
    rbounds = stv[[method]]$range * c( 0.01, 1.5 )
    phibounds = range( -log(0.05) / rbounds ) ## approximate
    nubounds = c(1e-3, stv[[method]]$kappa * 1.5 )# Finley et al 2007 suggest limiting this to (0,2)
    # Finley, Banerjee Carlin suggest that kappa_geoR ( =nu_spBayes ) > 2 are indistinguishable .. identifiability problems cause slow solutions

    starting = list( phi=median(phibounds), sigma.sq=0.5, tau.sq=0.5, nu=1  ) # generic start
    tuning   = list( phi=starting$phi/10, sigma.sq=starting$sigma.sq/10, tau.sq=starting$tau.sq/10, nu=starting$nu/10 ) # MH variance
    priors   = list(
      beta.flat = TRUE,
      phi.unif  = phibounds,
      sigma.sq.ig = c(5, 0.5), # inverse -gamma (shape, scale):: scale identifies centre; shape higher = more centered .. assuming tau ~ sigma
      tau.sq.ig = c(5, 0.5),  # inverse gamma (shape, scale) :: invGamma( 3,1) -> modal peaking < 1, center near 1, long tailed
      nu.unif = nubounds
    )

    model = spLM( z ~ 1, coords=as.matrix(xy), starting=starting, tuning=tuning, priors=priors, cov.model="matern",
      n.samples=n.samples, verbose=TRUE )

    ##recover beta and spatial random effects
    mm <- spRecover(model, start=burn.in*n.samples )
    mm.pred <- spPredict(mm, pred.covars=pCovars, pred.coords=as.matrix(locsout), start=burn.in*n.samples )
    res = apply(mm.pred[["p.y.predictive.samples"]], 2, mean)

    u = apply(mm$p.theta.recover.samples, 2, mean)
    vrange = geoR::practicalRange("matern", phi=1/u["phi"], kappa=u["nu"]  )

    spb = list( model=model, recover=mm,
      range=vrange, varSpatial=u["sigma.sq"], varObs=u["tau.sq"],  phi=1/u["phi"], kappa=u["nu"] )  # output using geoR nomenclature

    if (plotdata) {
      x11()
      # to plot variogram
      x = seq( 0, vrange* 2, length.out=100 )
      acor = geoR::matern( x, phi=1/u["phi"], kappa=u["nu"] )
      acov = u["tau.sq"] +  u["sigma.sq"] * (1- acor)  ## geoR is 1/2 of gstat and RandomFields gamma's
      plot( acov ~ x , col="orange", type="l", lwd=2, ylim=c(0,max(acov)*1.1) )
      abline( h=u["tau.sq"] + u["sigma.sq"]  )
      abline( h=u["tau.sq"] )
      abline( h=0 )
      abline( v=0 )
      abline( v=vrange )

      round(summary(mm$p.theta.recover.samples)$quantiles,2)
      round(summary(mm$p.beta.recover.samples)$quantiles,2)
      mm.w.summary <- summary(mcmc(t(mm$p.w.recover.samples)))$quantiles[,c(3,1,5)]

      plot(z, mm.w.summary[,1], xlab="Observed w", ylab="Fitted w",
          xlim=range(w), ylim=range(mm.w.summary), main="Spatial random effects")
      arrows(z, mm.w.summary[,1], w, mm.w.summary[,2], length=0.02, angle=90)
      arrows(z, mm.w.summary[,1], w, mm.w.summary[,3], length=0.02, angle=90)
      lines(range(z), range(z))

      obs.surf <- mba.surf(cbind(xy, z), no.X=500, no.Y=500, extend=T)$xyz.est
      image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
      points(xy)
      contour(obs.surf, add=T)

      x11()
      require(lattice)
      levelplot( res ~ locsout[,1] + locsout[,2], add=T )
    }

    return( res )
  }

}
   
