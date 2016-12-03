
  require(geostatsp)
  require(numDeriv)

  data(swissRain)
  str(swissRain)
  str(swissAltitude)
  names(swissRain)
  names(swissAltitude)


  # ML kriging with geostatsp

    # formula specifies fixed effects (right hand side (i.e., "CHE_alt") is vector in "data"; or layer/brick/stack in a Raster object passed in "covariates"; or name of one element if covariates is a list of Raster objects ... useful if rasters are of different resolutions/projections
    # locations=90 means raster with 90 cells in x direction, square cells; or as a Raster object
    # shape=1 # Matern shape parameter
    # boxcox=0.5 is square-root transform
    swissRes = lgm( data=swissRain, formula=rain~ CHE_alt, locations=90, covariates=swissAltitude,
                   shape=1, fixShape=TRUE, boxcox=0.5, fixBoxcox=TRUE, aniso=TRUE )
    
    swissRes = lgm( data=swissRain, formula=rain~ CHE_alt, locations=90, covariates=swissAltitude,
                   shape=0.2, fixShape=TRUE, boxcox=0.5, fixBoxcox=FALSE, aniso=TRUE, 
                   parscale=c(range=1000, shape=1, boxcox=1, nugget=0.1), reml=TRUE )



    # names(swissRes)
    # str(swissRes)  # "RasterStack"
    # parscale=c(range=10^4), param=c(range=10^5) 
 
    swissRes$summary
    
    names( swissRes$predict )

    plot(swissRes$predict[["predict"]], main="predicted rain") # "predict" is on the natural scale
    plot(swissRes$predict[["fixed"]], main="predicted rain") # "fixed" is the fixed effects
    plot(swissBorder, add=TRUE)
     
    # conditional probability that rainfall > 30 mm
    exc30 = excProb( swissRes, 30, nuggetInPrediction = TRUE )


    # simulating from the conditional distribution
    sim = geostatsp::RFsimulate( 
            model=swissRes$param, data=swissRes$resid["resid"], err.model=swissRes$param["nugget"],
            x=raster( extent(swissRain), nrow=100, ncol=100) )

    plot(sim)




    ---

install.packages("spdep")
install.packages("numDeriv")
install.packages("fields")
install.packages("rgdal")
install.packages("pixmap")
source("http://www.math.ntnu.no/inla/givemeINLA.R")
install.packages("geostatsinla",repos="http://r-forge.r-project.org")


    # Generalized Linear geostats models

    require(geostatsp)

    data(loaloa)

    plot( loaloa)
    plot( elevationLoa )
    plot( eviLoa )  # vegetation
    plot( ltLoa )  # landtype

    require(geostatsinla)


    eL = elevationLoa -750
    eLow  = reclassify( eL, c(0, Inf,0 ) )
    eHigh = reclassify( eL, c(-Inf, 0, 0 ) )

    rcl = rbind( c(9,8), # savanna to woody savanna
                c(5,2), c(11,2), # wetlands and mixed forest to forest
                c(12,14), c(13,14) ) # crop and urban to crop/natural mosaic
    (rcl)
    lTypes = reclassify( ltLoa, rcl )
    levels( lTypes ) = levels(ltLoa) 

    covars = list( elLow=eLow, elHigh=eHigh, land=lTypes, evi=eviLoa )

    names(loaloa)

    lfit =  glgm( loaloa, formula=y~factor(land) + evi + elHigh + elLow, 
                 family="binomial", Ntrials=loaloa$N, cells=150, covariates=covars, 
                 shape=1, buffer=50000, priorCI=list( sd=0.2,4), range=c(20000, 5e5) )
    # loaloa$y = number of positive samples
    # loaloa$N = number of samples taken in each village
    # prior 95% intervals ... glgm creates Gamma priors for precision 1/sigma^2 
    # and scaled range parameter phi/delta ; delta is cell size
    
    lft$paramaters$summary

    lfit$raster  # RasterStack with posterior means, SD, quantiles for the random effects U(s) and predicted values on link scale 
    plot( lfit$raster[["predict.invlogit"]] ) #posterior means lambda(s) 
    plot( lfit$raster[["random.mean"]]  ) # random effects means 

        

    




