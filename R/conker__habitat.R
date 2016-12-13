
conker__habitat = function( p, x, pa ) {
   #\\ this is the core engine of conker .. localised space-time habiat modelling
 
  if (0) {
    if (!exists("nsims", p)) p$nsims = 5000
    if (!exists("habitat.threshold.quantile", p)) p$habitat.threshold.quantile = 0.05 # quantile at which to consider zero-valued abundance
  }

  if ( exists("conker_local_model_distanceweighted", p) ) {
    if (p$conker_local_model_distanceweighted) {
      Hmodel = try( gam( p$conker_local_modelformula, data=x, family=binomial, weights=weights, optimizer=c("outer","optim")  ) )
    } else {
      Hmodel = try( gam( p$conker_local_modelformula, data=x, family=binomial, optimizer=c("outer","optim")  ) )
    }
  } else {
      Hmodel = try( gam( p$conker_local_modelformula, data=x, family=binomial ) )
  } 
  if ( "try-error" %in% class(Hmodel) ) return( NULL )


  if ( exists("conker_local_model_distanceweighted", p) ) {
    if (p$conker_local_model_distanceweighted) {
      Amodel = try( gam( p$conker_local_modelformula, data=x, weights=weights, optimizer=c("outer","optim")  ) )
    } else {
      Amodel = try( gam( p$conker_local_modelformula, data=x, optimizer=c("outer","optim")  ) )
    }
  } else {
      Amodel = try( gam( p$conker_local_modelformula, data=x ) )
  } 
  if ( "try-error" %in% class(Amodel) ) return( NULL )

  x$P = try( predict( Hmodel, newdata=x, type="response", se.fit=FALSE ) ) 
  x$A = try( predict( Amodel, newdata=x, type="response", se.fit=FALSE ) ) 
  x$Yhat = x$P * x$A

  rsq = cor( x$Yhat, x[,p$variables$Y], use="pairwise.complete.obs" )^2
  if (rsq < p$conker_rsquared_threshold ) return(NULL)

  Hmodel.coef = mvtnorm::rmvnorm(p$nsims, coef(Hmodel), Hmodel$Vp, method="chol")
  rm( Hmodel); gc()
  Hsim = family(M)$linkinv( predict(Hmodel, newdata=pa, type="lpmatrix") %*% t(Hmodel.coef) )
  rm( Hmodel.coef); gc()
  oops = which( is.na(Hsim) )
  if (length(oops) > 0)  Hsim[oops ] = 0  # assume to be zero
  
  pa$logitmean = apply( Hsim, 1, mean, na.rm=T )
  pa$logitsd = apply( Hsim, 1, sd, na.rm=T )

  Amodel.coef = mvtnorm::rmvnorm(p$nsims, coef(Amodel), Amodel$Vp, method="chol")
  rm(Amodel); gc()
  Asim = family(M)$linkinv( predict(Amodel, newdata=pa, type="lpmatrix") %*% t(Amodel.coef) )
  rm (Amodel.coef); gc()
  oops = which( is.na(Asim) )
  if (length(oops) > 0)  Asim[oops ] = 0  # assume to be zero
  
  # Do not extrapolate: trim to XX% quantiles to be a little more conservative
  oopu =  which( Asim > p$qs[2] )
  if (length(oopu) > 0)  Asim[ oopu ] = p$qs[2]

  oopl =  which( Asim < p$qs[1]  )
  if (length(oopl) > 0)  Asim[ oopl ] = 0  # below detection limits

  Asim = Asim * Hsim  # Asim now becomes weighted by Pr of habitat

  pa$mean = as.vector( apply( Asim, 1, mean, na.rm=T ) )
  pa$sd =  as.vector( apply( Asim, 1, sd, na.rm=T ) )

  # iHabitat = which( pa$logitmean > p$habitat.threshold.quantile & (pa$logitmean - 2 * pa$logitsd) > 0 )

  ss = summary(hmod)
  conker_stats = list( sdTotal=sd(Y[], na.rm=T), rsquared=rsq, ndata=nrow(x) ) # must be same order as p$statsvars

  return( list( predictions=pa, conker_stats=conker_stats ) )  

}
