
hivemod_variogram = function( xy, z, plotdata=FALSE, edge=c(1/3, 1), methods=c("fast"), maxdist=NA, nbreaks = 15, functionalform="matern", eps=1e-6 ) {

  #\\ estimate empirical variograms (actually correlation functions) and then model them using a number of different approaches .. mostly using Matern as basis
  #\\ returns empirical variogram and parameter estimates, and the models themselves
  #\\ expect xy = c(p/lon, p/lat), z= variable

  #\\ NOTE:: the default parameterization is as in spBayes and gstat which is:

  #\\ matern covariogram (||x||) = sigma^2 * (2^{nu-1} * Gamma(nu) )^{-1} * (phi*||x||)^{nu} * K_{nu}(phi*||x||)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is the range parameter  
  # -------------------------
  if ( 0 ) {
   # just for debugging / testing ... and example of access method:
   bioLibrary("bio.utilities", "bio.spacetime", "hivemod")
   require(sp)
   data(meuse)
    xy = meuse[, c("x", "y")]
    mz = log( meuse$zinc )
    mm = lm( mz ~ sqrt( meuse$dist ) )
    z = residuals( mm)

    plotdata=TRUE
    maxdist = NA
    edge=c(1/3, 1)
    nbreaks = 15

        # tests
    gr = hivemod_variogram( xy, z, methods="geoR" )
    gs = hivemod_variogram( xy, z, methods="gstat" )
    grf = hivemod_variogram( xy, z, methods="RandomFields" )
    gsp = hivemod_variogram( xy, z, methods="spBayes" )
    ginla = hivemod_variogram( xy, z, methods="inla" )

    # tests:
    out = gsp
    nd = nrow(out$spBayes$recover$p.theta.samples)
    rr = rep(NA, nd )
    for (i in 1:nd) rr[i] = geoR::practicalRange("matern", phi=1/out$spBayes$recover$p.theta.samples[i,3], kappa=out$spBayes$recover$p.theta.samples[i,4] )
    hist(rr)  # range estimate

    hist( out$spBayes$recover$p.theta.samples[,1] ) #"sigma.sq"
    hist( out$spBayes$recover$p.theta.samples[,2] ) # "tau.sq"
    hist( out$spBayes$recover$p.theta.samples[,3] ) # 1/phi
    hist( out$spBayes$recover$p.theta.samples[,4] ) # nu

    out = hivemod_variogram( xy, z )
    (out$geoR$range)
    out = hivemod_variogram( xy, z, nbreaks=30 )
    (out$geoR$range)

    out = hivemod_variogram( xy, log(z), nbreaks=30 )
    (out$geoR$range)
    out = hivemod_variogram( xy, log(z) )
    (out$geoR$range)
    require(mgcv)
    og = gam( log(z) ~ s( x) + s(y) + s(x,y), data=xy )
    zr = residuals(og)
    out = hivemod_variogram( xy, zr )  # remove spatial trend results in no variogram, as would be expected
    (out$geoR$range)
    og = gam( log(z) ~ s( elev ) , data=meuse )
    zr = residuals(og)
    out = hivemod_variogram( xy, zr )  # remove spatial trend results in no variogram, as would be expected
    (out$geoR$range)

    require(geoR)
    # plot( out$geoR$vgm )
    # lines( out$geoR$fit, lwd=2, col="slateblue" )
    xRange = c( 0, max(out$geoR$range*2.1 ) )
    yRange = c( 0, max(out$geoR$vgm$v )*1.05 )
    plot ( out$geoR$vgm$v ~ out$geoR$vgm$u, pch=20, xlim=xRange, ylim=yRange, ylab="Semivariance", xlab="Distance" )
      abline( h=0,  col="gray", lwd=2 )
      abline( h= (out$geoR$varSpatial + out$geoR$varObs), lty="dashed", col="slategray"  )
      abline( h=  out$geoR$varObs , lty="dashed", col="slategray")
      abline( v=out$geoR$range, lty="dotted", col="slateblue" )
      abline( v=0,  col="gray", lwd=2 )
      x = seq( 0, 2*out$geoR$range, length.out=100 )
      acor = geoR::matern( x, phi=out$geoR$phi, kappa=out$geoR$kappa  )
      acov = out$geoR$varObs +  out$geoR$varSpatial * (1- acor)
      lines( acov ~ x , col="blue", lwd=2 )
  } 



  nc_max = 5  # max number of iterations

  xy = as.data.frame(xy)
  names(xy) =  c("plon", "plat" ) # arbitrary
  xy_n = nrow(xy)

  out = list()
  out$varZ = var( z, na.rm=TRUE )  # this is the scaling factor for semivariance .. diving by sd, below reduces numerical floating point issues

  xr = range( xy$plon, na.rm=TRUE )
  yr = range( xy$plat, na.rm=TRUE )
  drange = min( diff( xr), diff( yr)  )

  # if max dist not given, make a sensible choice
  if ( is.na(maxdist)) {
    maxdist = drange * 0.1  # default
  } else if ( maxdist=="all") {
    maxdist = drange
  }
  out$drange = drange
  out$maxdist = maxdist

  # a small error term to prevent some errors in GRMF methods
  derr = out$maxdist / 10^6
  xy$plon = xy$plon + runif(xy_n, -derr, derr)
  xy$plat = xy$plat + runif(xy_n, -derr, derr)


  if ( "fast" %in% methods)  {
    # gives a fast stable empirical variogram

    require( RandomFields ) ## max likilihood
    RFoptions( allowdistanceZero=TRUE ) #, modus_operandi="easygoing" )

    rownames( xy) = 1:nrow(xy)  # seems to require rownames ...
    Yyy <- RFspatialPointsDataFrame( coords=xy, data=z, RFparams=list(vdim=1, n=1) )
    vario = RFempiricalvariogram( data=Yyy )
    
    # remove the (0,0) point -- force intercept
    todrop = which( !is.finite(vario@emp.vario )) # occasionally NaN's are created!
    todrop = unique( c(1, todrop) )
    vg = vario@emp.vario[-todrop]
    vx = vario@centers[-todrop]
    
    mvg = max(vg, na.rm=TRUE)
    mvx = max(vx, na.rm=TRUE)
    eps = 1e-6
    lower =c(eps,eps,eps, eps)
    upper =c(mvg, mvg, mvx, 2)
    #nonlinear est
    par = c(tau.sq=mvg*0.05, sigma.sq=mvg*0.95, phi=mvx/5, nu=0.5) 
    o = try( optim( par=par, vg=vg, vx=vx, method="L-BFGS-B", lower=lower, upper=upper,
      control=list(maxit=100),
      fn=function(par, vg, vx){ 
        vgm = par["tau.sq"] + par["sigma.sq"]*(1-fields::Matern(d=vx, range=par["phi"], smoothness=par["nu"]) )
        dy = sum( (vg - vgm)^2) # vario normal errors, no weights , etc.. just the line
      } ) 
    )
    
    if ( !inherits(o, "try-error")) { 
      if ( o$convergence==0 ) {
        par = o$par 
        out$fast = list( fit=o, vgm=vario, range=NA, nu=par["nu"], phi=par["phi"],
          varSpatial=par["sigma.sq"], varObs=par["tau.sq"] ) 
        rg = try(geoR::practicalRange("matern", phi=out$fast$phi, kappa=out$fast$nu ))
        if (! inherits(rg, "try-error") ) {
          out$fast$range = rg
        } else {
          out$fast$range = 0
        }
      } else {
        out = try( hivemod_variogram( xy=xy, z=z, methods="gstat") )
        if (!inherits(out, "try-error") )
        out$fast =out$gstat
        out$gstat =NULL
        return(out)
      }
    }
 
  
    if( 0) {
      scale = out$fast$phi * (sqrt(out$fast$nu*2) )  
      plot(vario, model=RMmatern( nu=out$fast$nu, var=out$fast$varSpatial, scale=scale) + RMnugget(var=out$fast$varObs) )

      RFoptions(seed=0, modus='easygoing' ) 
      rmod = ~ RMmatern( nu=NA, var=NA, scale=NA) + RMnugget(var=NA, scale=NA)

      rfit = RFfit(rmod, data=Yyy)

    }
  
    return(out)
  }




  # ------------------------
  if ("fields" %in% methods){
    require(fields)
    vg = vgram( xy, z, N=nbreaks, dmax=out$maxdist * 3 )
    smoothness =nu = 0.5 # 0.5 == exponential
    # theta = range paramter
    theta.grid = 10^seq( -6, 6, by=0.5) * out$maxdist
    lambda.grid = 10^seq( -9, 1, by=0.5) * out$maxdist
    
    res =NULL  
    fsp = try( MLESpatialProcess(xy, z, theta.grid=theta.grid, lambda.grid=lambda.grid,
      cov.function = "stationary.cov",  cov.args = list(Covariance = "Matern", smoothness = nu), 
      ngrid = 10, niter = 15, tol = 0.01, Distance = "rdist", nstep.cv = 50 ) )
    if (! inherits(fsp, "try-error") ) {
      res = fsp$pars 
    } else {
      fsp = MLE.Matern(xy, z, smoothness=nu, theta.grid =theta.grid )
      if( is.finite(sum(fsp$pars))) res = fsp$pars 
    }
    if (is.null(res)) return(NULL)

    vgm = Matern( d=vg$centers, range=res["theta"], smoothness=nu )    
    nugget = res["sigma"]^2 
    sill = res["rho"]  
    cvg = data.frame( cbind( x=vg$centers, cvgm= (nugget + sill * (1-vgm)) ))
    out$fields = list( fit=fsp, vgm=cvg, range=NA, nu=nu, phi=res["theta"] ,
        varSpatial=sill, varObs=nugget  )  # fields::"range" == range parameter == phi
    
    out$fields$range = geoR::practicalRange("matern", phi=out$fields$phi, kappa=out$fields$nu  )

    if( 0){
      x11()
      plot( cvg, type="b", ylim=range( c(0, cvg$cvgm) ) )
      abline( h=out$fields$varSpatial + out$fields$varObs)  
      abline( h=out$fields$varObs )
      abline( v=out$fields$range )
      
      x11()
      lambda.MLE<- fsp$pars["sigma"]^2/fsp$pars["rho"]  # ratio nugget / sill variance
      fsp2<- Krig( xy, z, Covariance="Matern", theta=fsp$pars["theta"], smoothness=nu, lambda= lambda.MLE)
      surface(fsp2)
      
      x11()
      fsp.p<- predictSurface(fsp2, lambda= lambda.MLE, nx=200, ny=200, )
      surface(fsp.p, type="I")
      
      x11()
      fsp.p2<- predictSurfaceSE(fsp2)
      surface(fsp.p2, type="C")
    }
  }
  

  # ------------------------

  if ("gstat" %in% methods){
    #\\ covariogram (||x||) = tau^2 * (2^{nu-1} * Gamma(nu) )^{-1} * (phi*||x||)^{nu} * K_{nu}(phi*||x||)
    #\\ gstat::kappa == spBayes::nu
    #\\ gstat::range == spBayes::phi {aka, "scale parameter"}

    require(gstat)
    require(sp)

    vrange = 0.5*out$maxdist # starting est of range
    distx = vrange * 0.9 ## back it up a bit to enter smoothly into the loop
    nc = 0
    while ( distx < vrange ) {
      nc = nc  + 1
      distx = distx * 1.25 # gradually increase distx until solution found
      vEm = try( variogram( z~1, locations=~plon+plat, data=xy, cutoff=distx, width=distx/nbreaks, cressie=TRUE ) ) # empirical variogram
      if (inherits(vEm, "try-error") ) return(NULL)
      vMod0 = vgm(psill=0.75, model="Mat", range=distx, nugget=0.25, kappa=1 ) # starting model parameters
      #vMod0 = vgm("Mat")
      vFitgs =  try( fit.variogram( vEm, vMod0, fit.kappa =TRUE, fit.sills=TRUE, fit.ranges=TRUE ) ) ## gstat's kappa is the Bessel function's "nu" smoothness parameter
      if (inherits(vFitgs, "try-error") )  return(NULL)
      vrange = max(1, geoR::practicalRange("matern", phi=vFitgs$range[2], kappa=vFitgs$kappa[2]  ) )
      if (nc > nc_max ) break()
    }
    if (inherits(vFitgs, "try-error") )  return(NULL)

    out$gstat = list( fit=vFitgs, vgm=vEm, range=NA, nu=vFitgs$kappa[2], phi=vFitgs$range[2],
        varSpatial=vFitgs$psill[2], varObs=vFitgs$psill[1]  )  # gstat::"range" == range parameter == phi
    
    out$gstat$range = geoR::practicalRange("matern", phi=out$gstat$phi, kappa=out$gstat$nu  )

    if (plotdata) {
      x11()
      plot(vEm, model=vFitgs, add=T)
      x11()
      plot( gamma ~ dist, data=out$gstat$vgm, xlim=c(0,out$maxdist), 
           ylim=c(0,max(out$gstat$vgm$gamma)*1.1), col="blue", pch=20 )
      abline( h=out$gstat$varSpatial + out$gstat$varObs ) 
      abline( h=out$gstat$varObs )
      abline( v=out$gstat$range )
      x = seq( 0, out$maxdist, length.out=100 )
      acor = geoR::matern( x, phi=out$gstat$phi, kappa=out$gstat$nu  )
      acov = out$gstat$varObs + out$gstat$varSpatial * (1- acor)
      lines( acov~x , col="red" )
    
      if (0) {
        # looks at the predictions
        gs <- gstat(id = "z", formula = z~1, locations=~plon+plat, data=xy, maxdist=distx, nmin=10, force=TRUE, model=vFitgs )
        # variogram of residuals
        data(meuse.grid)
        meuse.grid$plon = meuse.grid$x
        meuse.grid$plat = meuse.grid$y

        preds <- predict(gs, newdata=meuse.grid )
        spplot(preds)

      }
    }
    return(out)
  }


  # -------------------------

  if ("geoR" %in% methods) {
    # weighted least squares
    #  he Matern model (correlation function, rho) is defined as:
    #  rho(u;phi,kappa) =(2^(kappa-1) Gamma(kappa))^(-1) (u/phi)^kappa K_kappa(u/phi)
     # where phi and kappa are parameters and K_kappa(...) denotes the
     # modified Bessel function of the third kind of order kappa.  The
     # family is valid for phi > 0 and kappa > 0.
     #\\ default covariogram (||x||) = tau^2 * (2^{nu-1} * Gamma(nu) )^{-1} * (phi*||x||)^{nu} * K_{nu}(phi*||x||)
      #\\ geoR:: rho(h) = (1/(2^(kappa-1) * Gamma(kappa))) * ((h/phi)^kappa) * K_{kappa}(h/phi)

    require( geoR )
    vrange = 0.5 * out$maxdist
    distx = vrange * 0.9 ## back it up a bit to enter smoothly into the loop
    nc = 0
    while ( distx  < vrange ) {
      nc = nc + 1
      distx = distx * 1.25
      vEm = try( variog( coords=xy, data=z, uvec=nbreaks, max.dist=distx ) )
      if  (inherits(vEm, "try-error") )  return(NULL)
      vMod = try( variofit( vEm, nugget=0.5*out$varZ, kappa=0.5, cov.model="matern", ini.cov.pars=c(0.5*out$varZ, distx/4) ,
        fix.kappa=FALSE, fix.nugget=FALSE, max.dist=distx, weights="cressie" ) )
        # kappa is the smoothness parameter , also called "nu" by others incl. RF
      if  (inherits(vMod, "try-error") )  return(NULL)
      # maximum likelihood method does not work well with Matern
      ML = FALSE
      if (ML) {
        vMod = likfit( coords=xy, data=z, cov.model="matern", ini.cov.pars=vMod$cov.pars,
        fix.kappa=FALSE, fix.nugget=FALSE, kappa=vMod$kappa, nugget=vMod$nugget, lik.method = "REML" )
      }
     vrange = vMod$practicalRange
     if (nc > nc_max ) break()
    }
  
    out$geoR = list( fit=vMod, vgm=vEm, model=vMod, range=vMod$practicalRange,
              varSpatial= vMod$cov.pars[1], varObs=vMod$nugget, 
              nu=vMod$kappa,  phi=vMod$cov.pars[2] )

    if (plotdata) {
      # not rescaled ...
      x11()
      plot(vEm)
      lines(vMod)
      x11()
      plot( out$geoR$vgm )
      x11()
      plot( out$geoR$vgm$v ~ c(out$geoR$vgm$u), pch=20 , 
           xlim=c(0,out$maxdist), ylim=c(0, 1.25) )
      abline( h=out$geoR$varSpatial + out$geoR$varObs)  
      abline( h=out$geoR$varObs )
      abline( v=out$geoR$range )
      x = seq( 0, max(out$geoR$vgm$u*out$maxdist), length.out=100 )
      acor = geoR::matern( x, phi=out$geoR$phi, kappa=out$geoR$nu  )
      acov = out$geoR$varObs +  out$geoR$varSpatial * (1-acor)  ## geoR is 1/2 of gstat and RandomFields gamma's
      lines( acov ~ x , col="orange" )
    }

    return(out)

  }


  # -------------------------


  if ("RandomFields" %in% methods) {
    require( RandomFields ) ## max likilihood
    rownames( xy) = 1:nrow(xy)  # seems to require rownames ...
    rfdata <- RFspatialPointsDataFrame(
      coords = xy,
        data = z,
        RFparams=list(vdim=1, n=1)
    )
   # RandomFields:  Cov(h) = v * Orig(A*h/s) ; s=scale, h=dist, A=aniso, v=variance, Original model (data scale)

   # where nu > 0 and K_nu is the modified Bessel function of second kind and distance r >= 0 between two pointsd
   # The Matern covariance model is given by: C(h) = v * phi(A*h/s).
   #  Cov(r) = 2^{1- nu} Gamma(nu)^{-1} (sqrt{2nu} r)^nu K_nu(sqrt{2nu} r)
   # "phi" = sqrt{2nu}/?? NOt clear ...
   
   # RFoptions(
   #   allowdistanceZero=TRUE,
    #  modus_operandi="precise", #‘"careless"’,‘"sloppy"’, ‘"easygoing"’, ‘"normal"’, ‘"precise"’,        ‘"pedantic"’, ‘"neurotic"’
   #   bin_dist_factor=out$maxdist/2,
      #bins=nbreaks,
      #critical=TRUE, 
   #   approx_zero=0.05, #  Value below which a correlation is considered to be essentially zero.
   #   spConform=TRUE # FALSE is faster
    #)

    model = ~ RMmatern( nu=NA, var=NA, scale=NA) + RMnugget(var=NA, scale=NA)
    
    o = RFfit(model, data=rfdata, allowdistanceZero=TRUE )
    oo=summary(o)

    out$RandomFields = list ( fit=o, vgm=o[2], model=oo, range=NA,
              varSpatial=oo$param["value", "matern.var"],
              varObs=oo$param["value", "nugget.var"],
              phi=(oo$param["value", "matern.s"] )/(sqrt(oo$param["value", "matern.nu"]*2) ), 
              nu=oo$param["value", "matern.nu"], # RF::nu == geoR:: kappa (bessel smoothness param)
              error=NA )

    out$RandomFields$range = geoR::practicalRange("matern", phi=out$RandomFields$phi, kappa=out$RandomFields$nu  )

    if (plotdata) {
      x11()
      py = as.vector(out$RandomFields$vgm@emp.vario) 
      px = out$RandomFields$vgm@centers
      plot(  py ~ px, pch=20, ylim=c(0, 1.25) )
      abline( h=out$RandomFields$varSpatial + out$RandomFields$varObs  )
      abline( h=out$RandomFields$varObs )
      abline( v=out$RandomFields$range )

      x = seq( 0, max(px ), length.out=100 )
      acor = geoR::matern( x, phi=out$RandomFields$phi, kappa=out$RandomFields$nu  )
      acov = out$RandomFields$varObs  + out$RandomFields$varSpatial*(1- acor)
      lines( acov~x , col="red" )

      # compare with:
        data(meuse.grid)
        meuse.grid$plon = meuse.grid$x
        meuse.grid$plat = meuse.grid$y


      x11()
      plot(o)
      RFinterpolate(o, ...  )
    }
    return(out)

  }

  # -------------------------

  if ("spBayes" %in% methods) {
    # note spBayes::phi = 1/ gstat::phi  
    require(spBayes)
    library(MBA)
    require( geoR )
    geoR = hivemod_variogram( xy, z, methods="geoR" )
    rbounds = c( median( diff(  geoR$geoR$vgm$u) )/2, geoR$geoR$range *1.5 )
    phibounds = range( -log(0.05) / rbounds ) ## approximate
    nubounds = c(1e-3, geoR$geoR$nu * 1.5 )# Finley et al 2007 suggest limiting this to (0,2)
    # Finley, Banerjee Carlin suggest that kappa_geoR ( =nu_spBayes ) > 2 are indistinguishable .. identifiability problems cause slow solutions
    n.samples = 5000
    starting = list( phi=median(phibounds), sigma.sq=0.51, tau.sq=0.51, nu=1.1  ) # generic start
    #starting = list( phi=1/2, sigma.sq=res.geoR$geoR$varSpatial, tau.sq=res.geoR$geoR$varObs, nu=30  ) # generic start
    tuning   = list( phi=starting$phi/12, sigma.sq=starting$sigma.sq/12, tau.sq=starting$tau.sq/12, nu=starting$nu/12 ) # MH variance to get acceptance rante bet 30-40%
    priors   = list(
      beta.flat = TRUE,
      phi.unif  = phibounds,
      sigma.sq.ig = c(5, 0.5),  # inverse -gamma (shape, scale):: scale identifies centre; shape higher = more centered .. assuming tau ~ sigma
      tau.sq.ig = c(5, 0.5),  # inverse gamma (shape, scale) :: invGamma( 3,1) -> modal peaking < 1, center near 1, long tailed
      nu.unif = nubounds
    )

    model = spLM( z ~ 1, coords=as.matrix(xy), starting=starting, tuning=tuning, priors=priors, cov.model="matern",
      n.samples=n.samples, verbose=TRUE )

    burn.in <- 0.2*n.samples

    ##recover beta and spatial random effects
    m.1 <- spRecover(model, start=burn.in )

    u = apply(m.1$p.theta.recover.samples, 2, mean)
    u["phi"] = 1/u["phi"]

    vrange = geoR::practicalRange("matern", phi=u["phi"], kappa=u["nu"]  )

    out$spBayes = list( model=model, recover=m.1,
      range=vrange, varSpatial=u["sigma.sq"], varObs=u["tau.sq"], 
      phi=u["phi"], nu=u["nu"] )  # output using geoR nomenclature

    if (plotdata) {
      x11()
      x = seq( 0, vrange* 2, length.out=100 )
      acor = geoR::matern( x, phi=u["phi"], kappa=u["nu"] )
      acov = u["tau.sq"] +  u["sigma.sq"] * (1- acor)  ## geoR is 1/2 of gstat and RandomFields gamma's
      plot( acov ~ x , col="orange", type="l", lwd=2, ylim=c(0,max(acov)*1.1) )
      abline( h=u["tau.sq"] + u["sigma.sq"]  )
      abline( h=u["tau.sq"] )
      abline( h=0 )
      abline( v=0 )
      abline( v=vrange )

      if(0) {
        round(summary(m.1$p.theta.recover.samples)$quantiles,2)
        round(summary(m.1$p.beta.recover.samples)$quantiles,2)
        m.1.w.summary <- summary(mcmc(t(m.1$p.w.recover.samples)))$quantiles[,c(3,1,5)]

        plot(z, m.1.w.summary[,1], xlab="Observed w", ylab="Fitted w",
            xlim=range(z), ylim=range(m.1.w.summary), main="Spatial random effects")
        arrows(z, m.1.w.summary[,1], z, m.1.w.summary[,2], length=0.02, angle=90)
        arrows(z, m.1.w.summary[,1], z, m.1.w.summary[,3], length=0.02, angle=90)
        lines(range(z), range(z))

        par(mfrow=c(1,2))
        obs.surf <-   mba.surf(cbind(xy, z), no.X=100, no.Y=100, extend=T)$xyz.est
        image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
        points(xy)
        contour(obs.surf, add=T)
      }
    }

    return(out)

  }


  # -------------------------


  if ("inla" %in% methods){
    require(INLA)
    require(lattice)

    inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors

    locs0  = as.matrix( xy )
    xy$b0 = 1  # intercept for inla

    vRange = out$maxdist * 0.1

    M0.domain = inla.nonconvex.hull( locs0 )
    MESH = inla.mesh.2d (
      loc=locs0, # locations of data points
      boundary = M0.domain,
      max.edge = edge * vRange
    )

#    kappa0 = sqrt(8) / vRange
#    tau0 = 1/ ( sqrt(4*pi) * kappa0 * vPsill )
    alpha = 2 # -> alpha-1 == nu (inla fixes it at 1,2,or 3)
    SPDE = inla.spde2.matern( MESH, alpha=alpha )
    spatial.field <- inla.spde.make.index('spatial.field', n.spde=SPDE$n.spde )

    # projection matrix A to translate from mesh nodes to data nodes
    A = inla.spde.make.A( mesh=MESH, loc=locs0 )

    # data stack for occurence (PA)
    Z = inla.stack(
        tag="data",
        data=list( z=z ) ,
        A=list(A, 1 ),
        effects=list( spatial.field=spatial.field, xy )  # b0 is the intercept
    )

    RES <- inla(  z ~ 0 + b0+ f( spatial.field, model=SPDE ), family="gaussian",
        data=inla.stack.data(Z),
        # control.compute=list(dic=TRUE),
        control.results=list(return.marginals.random=TRUE ),
        # control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        # control.fixed = list(expand.factor.strategy='inla') ,
        control.predictor=list(A=inla.stack.A(Z), compute=TRUE, link=1 ) ,
        # control.inla = list( h=1e-4, tolerance=1e-10),
        # control.inla=list(strategy="laplace", npoints=21, stencil=7 ) ,
        verbose = FALSE
    )

    oo = inla.spde2.result(RES, "spatial.field", SPDE, do.transf=TRUE)

    inames = c( "mode", "mean", "sd", "quant0.025", "quant0.25", "quant0.5",  "quant0.75", "quant0.975", "low", "high" )

    # Range parameter .. ie, sqrt(8)/exp(oo$summary.log.kappa$mean)
    im = oo$marginals.range.nominal[[1]]
    iRange = c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im )) )

    # "Spatial variance/error ('partial sill variance')"
    im = oo$marginals.variance.nominal[[1]]
    iVar =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im )) )

    # kappa.inla  == 1/phi.geoR
    im = oo$marginals.kappa[[1]]
    iKappa =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im ) ) )

    # tau
    im = oo$marginals.tau[[1]]
    iTau =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im ) ) )

    ## Non-spatial ("observation") error ('nugget variance')
    iprec = grep ( "Precision.*observ.*", names(RES$marginals.hyperpar), ignore.case=TRUE )
    im = inla.tmarginal( function(x) {1/x}, RES$marginals.hyperpar[[ iprec ]] )
    iNugget =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im ) ) )

    inla.summary = as.matrix( rbind( iKappa, iTau, iRange, iVar, iNugget ) )
    rownames( inla.summary) = c( "kappa", "tau", "range", "spatial error", "observation error" )
    colnames( inla.summary) = inames

    out$inla = list(summary=inla.summary, 
      mesh=MESH, res=RES, range.inla90=inla.summary[["range","mean"]] ,
      varSpatial=inla.summary[["spatial error","mean"]] , 
      varObs=inla.summary[["observation error","mean"]] ,
      phi = 1/inla.summary[["kappa","mean"]] , nu=alpha-1, error=NA )

    # kappa{geoR} = lambda{INLA} == alpha-1 {INLA} and alpha=2 by default in INLA
    out$inla$range = geoR::practicalRange("matern", phi=out$inla$phi, kappa=out$inla$nu  )


    if (plotdata) {
      require( geoR )
      x = seq( 0,  out$inla$range * 1.5, length.out=100 )
      svar =  out$inla$varObs + out$inla$varSpatial * (1-geoR::matern( x, phi=out$inla$phi, kappa=out$inla$nu  ))
      plot( svar~x, type="l" )
      abline( h=out$inla$varObs + out$inla$varSpatial )
      abline( h=out$inla$varObs )
      abline( v=out$inla$range.inla90  )
      abline( v=out$inla$range, col="red"  )
    }
  
    return(out)

  }


  # -------------------------

  if ("bayesx" %in% methods){
    library("R2BayesX")
    # by default, bayesx fixes nu=1.5  , see: bayesx.term.options( bs="kr", method="REML" )
    # phi = max(distance) / const, such that Corr(distance=const) = 0.001; 
    # i.e. range at distance where covar ~0.999 .. but not sure how to recover the correct phi/range from this ...
    nu = 1.5 

    fm <- bayesx( z ~ sx(plon, plat, nu=nu,  bs="kr" ), family="gaussian", method="REML", data =xy )

    out$bayesx = list( fit=fitted(fm), range=NA, model=fm, 
        varSpatial = fm$smooth.hyp[,"Variance"] , 
        varObs = fm$fixed.effects[1,"Std. Error"] , 
        nu =nu, phi=1/fm$smooth.hyp[,"Smooth Par."] )

    out$bayesx$range = geoR::practicalRange("matern", phi=out$bayesx$phi, kappa=out$bayesx$nu  )
    out

    if(0){
      plot( fm, term = "sx(plon,plat)", image=TRUE, contour=TRUE )
      # summary(fm)
      # lattice::levelplot( out ~ plon+plat, data=k, aspect="iso" )
      x = seq( 0, out$maxdist, length.out=100 )
      acor = geoR::matern( x, phi=out$bayesx$phi, kappa=out$bayesx$nu  )
      acov = out$bayesx$varObs  + out$bayesx$varSpatial*(1- acor)
      plot( acov~x , col="red", type="b", xlim=c(0, out$maxdist * 1.5), 
        ylim=c(0,(out$bayesx$varSpatial + out$bayesx$varObs) *1.5) )
      abline( h=out$bayesx$varSpatial + out$bayesx$varObs  )
      abline( h=out$bayesx$varObs )
      abline( v=out$bayesx$range )
    }

    return(out)
  }


  # -------------------------

  if ("jags" %in% methods){
    require(rjags)
    require(jagsUI)
    # assume nu = 1 (due to identifiability issues)

    print( "Slow ... 7.5 min for meuse test data" )

    jagsmodel = paste0("
    model{
      for(i in 1:N){
        y[i] ~ dnorm(mu[i], prec)
        mu[i] = beta0 + errorSpatial[i]
        muSpatial[i] = 0
      }
      prec = 1.0/ (tausq + sigmasq )
      invCOVAR = inverse(COVAR)
      errorSpatial ~ dmnorm( muSpatial, invCOVAR)
      for(i in 1:N) {
        COVAR[i,i] = sigmasq
        for(j in 1:(i-1)) {
          COVAR[i,j] = sigmasq * exp(-( DIST[i,j]/phi))
          COVAR[j,i] = COVAR[i,j]
        } 
      }
      tausq = 1/tausq.inv
      tausq.inv ~ dgamma(0.1,0.1)
      sigmasq = 1/sigmasq.inv
      sigmasq.inv ~ dgamma(2,1)
      phi ~ dgamma(1,0.1)
      beta0 ~ dnorm(0,0.0001)
    } ")

  fn = tempfile()
  cat( jagsmodel, file=fn )
  distances = as.matrix(dist( xy, diag=TRUE, upper=TRUE))
  Data = list( N=length(z), DIST=distances, y=z )
  fit = jagsUI::jags(data=Data, 
       parameters.to.save=c("phi", "sigmasq", "tausq"),
       model.file=fn,
       n.iter=1000,
       n.chains=3,
       n.burnin=100,
       n.thin=5,
       parallel=TRUE,
       DIC=FALSE)
   
   if (0) {
     summary(fit)
     plot(fit)
     gelman.plot(fit$samples)
    # geweke.plot(fit$samples)
    #update(fit, n.iter=2000, n.thin=20 )
      acf( fit$sims.list$phi)
      acf( fit$sims.list$sigmasq)
      acf( fit$sims.list$tausq)

# JAGS output for model '/tmp/RtmpXoFpO7/file539863b2591e', generated by jagsUI.
# Estimates based on 3 chains of 2000 iterations,
# burn-in = 200 iterations and thin rate = 10,
# yielding 540 total samples from the joint posterior. 
# MCMC ran in parallel for 7.461 minutes at time 2016-08-10 18:05:54.

#          mean    sd  2.5%   50% 97.5% overlap0 f  Rhat n.eff
# phi     1.252 0.597 0.513 1.094 2.784    FALSE 1 1.014   238
# sigmasq 0.441 0.084 0.292 0.436 0.621    FALSE 1 1.008   356
# tausq   0.145 0.096 0.030 0.121 0.374    FALSE 1 1.013   279
# $phi
# [1] 222.3888

# $sigmasq
# [1] 0.08303909

# $tausq
# [1] 0.02724344

# $range
# [1] 889.2267

    }
    print (summary(fit))

    out$jags = list(
      fit = fit, model=jagsmodel,
      phi = 1 / fit$summary["phi", "mean"],
      sigmasq = fit$summary["sigmasq", "mean"] ,
      tausq = fit$summary["tausq", "mean"] 
    )
    out$jags$range = geoR::practicalRange("matern", phi=out$jags$phi, kappa=1  )
    return(out)
  }


  # -------------------------

  if ("TMB" %in% methods){
   
  }

  # -------------------------

  if ("LaplacesDemon" %in% methods){
    
    cat("TODO:: modify to RCPP some of the 
        1) matrix methods and \n 
        2) optimizers for more speed and even \n 
        3) LA.bio for streamlining and \n 
        4) dmvn sparse form? \n")

    require(LaplacesDemonCpp)
    
    Data = list(
      eps = 1e-6,
      N = length(z),  # required for LaplacesDemon
      DIST=as.matrix(dist( xy, diag=TRUE, upper=TRUE)), # distance matrix between knots
      y=z  
    )
    Data$mon.names = c( "LP", paste0("yhat[",1:Data$N,"]" ) )
    Data$parm.names = as.parm.names(list(tausq=0, sigmasq=0, phi=0, nu=0 ))
    Data$pos = list(
      tausq = grep("tausq", Data$parm.names),
      sigmasq = grep("sigmasq", Data$parm.names),
      phi = grep("phi", Data$parm.names),
      nu = grep("nu", Data$parm.names)
    )
    Data$PGF = function(Data) {
      #initial values .. get them near the center of mass
      tausq = rgamma (1, 1, 5) # 0 to 1.5 range
      sigmasq = rgamma (1, 1, 5)
      phi = rgamma (1, 1, 1)  # 0 to 500 range
      nu = runif(1, 0.5, 4)
      return( c( tausq, sigmasq, phi, nu ))
    }
    Data$PGF  = compiler::cmpfun(Data$PGF)
    Data$Model = function(parm, Data) {
      tausq = parm[Data$pos$tausq] = LaplacesDemonCpp::interval_random(parm[Data$pos$tausq], Data$eps, 1, 0.01 )
      sigmasq = parm[Data$pos$sigmasq]= LaplacesDemonCpp::interval_random(parm[Data$pos$sigmasq], Data$eps, 1, 0.01 )
      phi = parm[Data$pos$phi]= LaplacesDemonCpp::interval_random(parm[Data$pos$phi], Data$eps, Inf, 0.01 )
      nu = parm[Data$pos$nu] = LaplacesDemonCpp::interval_random(parm[Data$pos$nu], 0.1, 15.0, 0.01 )
      # corSpatial = exp(-Data$DIST/phi)^nu   ## spatial correlation .. simple exponential model
      # corSpatial = geoR::matern( Data$DIST, phi=1/phi, kappa=nu )   ## spatial correlation .. matern from geoR
      # wikipedia Matern parameterization: 
      # C_{\nu }(d) = \sigma ^{2}{\frac {2^{1-\nu }}{\Gamma (\nu )}} 
      #  {\Bigg (}{\sqrt {2\nu }}{\frac {d}{\rho }}{\Bigg )}^{\nu } K_{\nu }
      #  {\Bigg (}{\sqrt {2\nu }}{\frac {d}{\rho }}{\Bigg )}
      e <- sqrt(2*nu) * Data$DIST / phi
      corSpatial = {2^{1-nu}}/gamma(nu) * (e^nu) * besselK(x=e, nu=nu) 
      diag(corSpatial) = 1
      # corSpatial = zapsmall(corSpatial)
      if ( !is.positive.definite(corSpatial)) {
        cat("correlation matrix is not positive definite, adding a bit of noise ...\n")
        corSpatial = as.positive.definite(corSpatial) 
        # browser()
      }

      eSp = rmvn( 1, rep(0, Data$N), sigmasq*corSpatial )# psill
      eObs = rnorm( Data$N, 0, sqrt(tausq) ) # nugget error
      tausq.prior = dgamma(tausq, 1, 1, log=TRUE) # 0-1.55 range
      sigmasq.prior = dgamma(sigmasq, 1, 1, log=TRUE)
      phi.prior = dgamma(phi, 1, 1, log=TRUE)
      nu.prior = dnorm(nu, 10, 0.1, log=TRUE)
      yhat = eObs + eSp # local iid error + spatial error
      LL = sum(dnorm(Data$y, yhat, sqrt(sigmasq+tausq), log=TRUE)) ## Log Likelihood
      LP = sum(LL, sigmasq.prior, tausq.prior, phi.prior, nu.prior) ### Log-Posterior
      Modelout = list(LP=LP, Dev=-2*LL, Monitor=c(LP, yhat), yhat=yhat, parm=parm)
      return(Modelout)
    }
    Data$Model = compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed
    
    parm0=Data$PGF(Data)
  
    f = LaplaceApproximation(Data$Model, Data=Data, parm=parm0, Method="BFGS", Iterations=1000, CPUs=4, Stop.Tolerance=1.0E-9 )


    if (plotdata) {

      f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Fit.LA=f,
        Iterations=10000, Thinning=100, Status=1000, Covar=f$Covar, CPUs=8 )

      parm0 = as.initial.values(f)
      f0 = LaplacesDemon(Data$Model, Data=Data, Initial.Values=parm0, CPUs=8 )
      mu = f$Summary1[,1]
      f0 = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Fit.LA=f,
        Iterations=5000, Thinning=1, Status=1000, Algorithm="IM", Specs=list(mu=mu), 
        Covar=f$Covar, CPUs=8 )

      f0 = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Fit.LA=f,
        Iterations=10000, Thinning=100, Status=1000, Covar=f$Covar, CPUs=8 )

      Consort(f0)
      plot(f0, Data=Data)
      m = f0$Summary2[grep( "\\<yhat\\>", rownames( f0$Summary2 ) ),]
      m = f$Summary2[grep( "\\<yhat\\>", rownames( f$Summary2 ) ),]
      # m = f$Summary2[grep( "muSpatial", rownames( f$Summary2 ) ),]
      plot( Data$y ~ m[, "Mean"]  )


    }

    out$LaplacesDemon = list( fit=f, vgm=NA, model=Data$Model, range=NA,
      varSpatial=f$Summary2["sigmasq", "Mean"] , 
      varObs=f$Summary2["tausq", "Mean"], 
      nu=f$Summary2["nu", "Mean"],  
      phi = ( f$Summary2["phi", "Mean"]  / sqrt(2*f$Summary2["nu", "Mean"] ) ) 
    )   ## need to check parameterization...
 
    out$LaplacesDemon$range = geoR::practicalRange("matern", phi=out$LaplacesDemon$phi, kappa=out$LaplacesDemon$nu)
 
   # print( out$LaplacesDemon )


    if (plotdata) {
      x11()
      x = seq( 0,  out$LaplacesDemon$range * 1.25, length.out=100 )
      svar =  out$LaplacesDemon$varObs + out$LaplacesDemon$varSpatial * (1-geoR::matern( x, phi=out$LaplacesDemon$phi, kappa=out$LaplacesDemon$nu  ))
      plot( svar~x, type="l", ylim=c(0, max(svar)) )
      abline( h=out$LaplacesDemon$varObs + out$LaplacesDemon$varSpatial )
      abline( h=out$LaplacesDemon$varObs )
      abline( v=out$LaplacesDemon$range, col="red"  )
    }
  
  return(out)
  
  }

}


