 
  lstfilter_summary_inla_spde2 = function( RES, SPDE, spatial.field.name="spatial.field" ) {
    #\\ extract summary stats from SPDE results of the random field 
    #\\ parameters are  on user scale, low/high is 95% HPD

    oo = inla.spde2.result(RES, spatial.field.name, SPDE, do.transf=TRUE)
   
    inames = c( "mode", "mean", "sd", "quant0.025", "quant0.25", "quant0.5",  "quant0.75", "quant0.975", "low", "high" )

    # Range parameter .. ie, sqrt(8)/exp(oo$summary.log.kappa$mean) 
    im = oo$marginals.range.nominal[[1]]
    iRange = c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im )) )

    # "Spatial variance/error ('partial sill variance')"
    im = oo$marginals.variance.nominal[[1]]
    iVar =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im )) )
    
    # kappa
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

    return( inla.summary )
  } 
  

