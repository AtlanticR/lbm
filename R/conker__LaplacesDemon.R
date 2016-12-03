
conker__LaplacesDemon = function( p, x, pa ) {
   #\\ this is the core engine of conker .. localised space-time habiat modelling
 
    if ( exists("conker_local_model_distanceweighted", p) ) {
      if (p$conker_local_model_distanceweighted) {
        Hmodel = try( gam( p$conker_local_modelformula, data=x, weights=weights, optimizer=c("outer","optim")  ) )
      } else {
        Hmodel = try( gam( p$conker_local_modelformula, data=x, optimizer=c("outer","optim")  ) )
      }
    } else {
        Hmodel = try( gam( p$conker_local_modelformula, data=x ) )
    } 
    if ( "try-error" %in% class(Hmodel) ) return( NULL )


    out = try( predict( Hmodel, newdata=newdata, type="response", se.fit=T ) )

    if ( "try-error" %in% class( out ) ) return( NULL )

    newdata$mean = as.vector(out$fit)
    newdata$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

    if (exists( "conker_quantile_bounds", p)) {
      tq = quantile( x[,p$variables$Y], probs=p$conker_quantile_bounds, na.rm=TRUE  )
      bad = which( newdata$mean < tq[1] | newdata$mean > tq[2]  )
      if (length( bad) > 0) {
        newdata$mean[ bad] = NA
        newdata$sd[ bad] = NA
      }
    }

    ss = summary(Hmodel)
    conker_stats = list( sdTotal=sd(Y[], na.rm=T), rsquared=ss$r.sq, ndata=ss$n ) # must be same order as p$statsvars

    return( list( predictions=newdata, conker_stats=conker_stats ) )


    if(0) {

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
          varSpatial=f$Summary2["sigmasq", "Mean"] *out$varZ,
          varObs=f$Summary2["tausq", "Mean"]*out$varZ,
          nu=f$Summary2["nu", "Mean"],
          phi = out$maxdist * ( f$Summary2["phi", "Mean"]  / sqrt(2*f$Summary2["nu", "Mean"] ) )
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
