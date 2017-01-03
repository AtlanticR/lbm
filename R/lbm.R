
lbm = function( p, DATA,  storage.backend="bigmemory.ram", continue=FALSE) {
  #\\ localized modelling of space and time data to predict/interpolate upon a grid OUT
  #\\ overwrite = FALSE restarts from a saved state
  #\\ speed ratings: bigmemory.ram (1), ff (2), bigmemory.filebacked (3)
  #\\ do a more relaxed second stage: turned off for now .. this would make sense to try but 
  #\\ it can be costly in terms of time in production runs .. use only for research purposes

  # TODO: splancs::kernel3d as a method ? .. for count data?
  # TODO: habitat methods
  # TODO: gaussian process // see lbm_interpolate_xy_simple
  #       .. the convoSPAT seems almost fast enough
  # TODO: LaplacesDemon method
  # TODO: look at bayesX a little more carefully.
  # TODO: MBA mba.surf method? ... seems very fast

  if(0) {

    p = bio.bathymetry::bathymetry.parameters( )
    p$lbm_local_modelengine = "krige"  # about 5 X faster than bayesx-mcmc method
    p$storage.backend="bigmemory.ram"
    p = bio.bathymetry::bathymetry.parameters( p=p, DS="lbm" )
    continue=FALSE
    DATA='bathymetry.db( p=p, DS="lbm.inputs" )'
  

    p = bio.substrate::substrate.parameters() # reset to defaults
    p$lbm_local_modelengine = "krige" 
    p$storage.backend="bigmemory.ram"  # filebacked metods are still too slow ..
    p = bio.bathymetry::bathymetry.parameters( p=p, DS="lbm" )
    continue=FALSE
    p$clusters = rep("localhost",  detectCores() )
    DATA = 'substrate.db( p=p, DS="substrate.lbm" )'
  
  
    p = bio.temperature::temperature.parameters( current.year=2016 )
    p$lbm_local_modelengine="gam"  # smoothest
    # p$lbm_local_modelengine="twostep" 
    # p$lbm_local_modelengine="uked" # (universal) kriging with external drift
    p$storage.backend="bigmemory.ram"
    p = bio.temperature::temperature.parameters( DS="lbm", p=p )
    continue=FALSE
    DATA='hydro.db( p=p, DS="temperature.lbm" )'

  }


  if (continue | exists( "ptr", p) ) {
  
    message( "Continuing from an interrupted start. This only works if your data objects still exist in RAM or on disk. " )
    message( "There is no guarantee that this is the case and so your mileage will vary ... ")
    if (exists("storage.backend", p) && p$storage.backend !="bigmemory.ram" ) p = lbm_db( p=p, DS="load.parameters" )  # ie. restart with saved parameters
    RLibrary( p$libs )
    lbm_db(p=p, DS="statistics.status.reset" )

  } else {


    p$time.start =  Sys.time()

    p$savedir = file.path(p$project.root, "lbm", p$spatial.domain )
    message( paste( "In case something should go wrong, intermediary outputs will be placed at:", p$savedir ) )
    if ( !file.exists(p$savedir)) dir.create( p$savedir, recursive=TRUE, showWarnings=FALSE )

    p$stloc = file.path( p$project.root, "lbm", p$spatial.domain, "tmp" )
    message( paste( "Temporary files are being created at:", p$stloc ) )
    if ( !file.exists(p$stloc)) dir.create( p$stloc, recursive=TRUE, showWarnings=FALSE )

    p$libs = unique( c( p$libs, "sp", "rgdal", "parallel", "RandomFields", "geoR" ) )

    if (!exists("storage.backend", p)) p$storage.backend = storage.backend
    if (any( grepl ("ff", p$storage.backend)))         p$libs = c( p$libs, "ff", "ffbase" )
    if (any( grepl ("bigmemory", p$storage.backend)))  p$libs = c( p$libs, "bigmemory" )

    if (p$lbm_local_modelengine=="bayesx")  p$libs = c( p$libs, "R2BayesX" )
    if (p$lbm_local_modelengine %in% c("gam", "mgcv", "habitat") )  p$libs = c( p$libs, "mgcv" )
    if (p$lbm_local_modelengine %in% c("LaplacesDemon") )  p$libs = c( p$libs, "LaplacesDemonCpp" )
    if (p$lbm_local_modelengine %in% c("inla") )  p$libs = c( p$libs, "INLA" )
    if (p$lbm_local_modelengine %in% c("fft", "gaussianprocess2Dt") )  p$libs = c( p$libs, "fields" )
    if (p$lbm_local_modelengine %in% c("gaussianprocess") )  p$libs = c( p$libs  )
    if (p$lbm_local_modelengine %in% c("spate") )  p$libs = c( p$libs, "spate" )
    if (p$lbm_local_modelengine %in% c("splancs") )  p$libs = c( p$libs, "splancs" )
    if (p$lbm_local_modelengine %in% c("twostep") )  p$libs = c( p$libs, "mgcv", "fields" )
    if (p$lbm_local_modelengine %in% c("krige") ) {
      if (p$lbm_krige_engine %in% c("default", "fields")) p$libs = c( p$libs, "fields" )
      if (p$lbm_krige_engine %in% c("gstat")) p$libs = c( p$libs, "gstat" )
    }  


    p$libs = unique( p$libs )
    RLibrary( p$libs )

    if (p$storage.backend=="bigmemory.ram") {
      if ( length( unique(p$clusters)) > 1 ) stop( "More than one unique cluster server was specified .. the RAM-based method only works within one server." )
    }

    p = lbm_parameters(p=p) # fill in parameters with defaults where possible
    p = lbm_db( p=p, DS="filenames" )
    p$ptr = list() # location for data pointers

    # set up the data and problem using data objects

    tmpfiles = unlist( p$cache)
	  for (tf in tmpfiles) if (file.exists( tf)) file.remove(tf)

    p$variables$all = NULL
    if (exists("lbm_local_modelformula", p))  {
      p$variables$local_all = all.vars( p$lbm_local_modelformula )
      p$variables$local_cov = intersect( p$variables$local_all, p$variables$COV ) 
      p$variables$all = unique( c( p$variables$all, p$variables$local_all ) )
    }
    if (exists("lbm_global_modelformula", p)) {
      p$variables$global_all = all.vars( p$lbm_global_modelformula )
      p$variables$global_cov = intersect( p$variables$global_all, p$variables$COV )      
      p$variables$all = unique( c( p$variables$all, p$variables$global_all ) )
    }

    # permit passing a function rather than data directly .. less RAM usage in parent call
    if (class(DATA)=="character") assign("DATA", eval(parse(text=DATA) ) )
    testvars = c(p$variables$Y, p$variables$COV, p$variables$TIME, p$variables$LOC)
    withdata = which(is.finite( (rowSums(DATA$input[, testvars] )) ) )
    if (length(withdata) < 1) stop( "Missing data or insufficient data")
    DATA$input = DATA$input[withdata, ]
    rm(withdata)

    # number of time slices
    if (!exists("nt", p)) {
      p$nt = 1  # default to 1 == no time
      if (exists( "ny", p)) p$nt = p$nt * p$ny  # annual time slices
      if (exists( "nw", p)) p$nt = p$nt * p$nw  # sub-annual time slices
    }

    # prediction times  for 2D methods, treat time as independent timeslices
    if ( !exists("ts", p)) p$ts = 1

    # require knowledge of size of stats output before create S, which varies with a given type of analysis
    othervars = c( )
    if (p$lbm_local_modelengine == "habitat") othervars = c( )
    if (exists("TIME", p$variables) )  othervars = c( "ar_timerange", "ar_1" )
    p$statsvars = unique( c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "range", "phi", "nu", othervars ) )

    message( "Initializing temporary storage of data and outputs (will take a bit longer on NFS clusters) ... ")
    message( "These are large files (4 to 6 X 5GB), esp. prediction grids (5 min .. faster if on fileserver), so be patient. ")
    lbm_db( p=p, DS="cleanup" )

    if (exists("lbm_global_modelengine", p)) {
      # to add global covariate model ??  .. simplistic this way but faster ~ kriging with external drift
      lbm_db( p=p, DS="global_model.redo", B=DATA$input )
    }

    # NOTE:: must not sink the following memory allocation into a deeper funcion as 
    # NOTE:: bigmemory RAM seems to lose the pointers if they are not made simultaneously 
    # init output data objects
    # statistics storage matrix ( aggregation window, coords ) .. no inputs required
    sbox = list( 
      plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$lbm_distance_statsgrid ),
      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$lbm_distance_statsgrid ) )

      # statistics coordinates
      Sloc = as.matrix( expand.grid( sbox$plons, sbox$plats ))
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Sloc = big.matrix(nrow=nrow(Sloc), ncol=ncol(Sloc), type="double"  )
          p$bm$Sloc[] = Sloc
          p$ptr$Sloc  = bigmemory::describe( p$bm$Sloc  )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Sloc  = p$cache$Sloc
          bigmemory::as.big.matrix( Sloc, type="double", backingfile=basename(p$bm$Sloc), descriptorfile=basename(p$cache$Sloc), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Sloc = ff( Sloc, dim=dim(Sloc), file=p$cache$Sloc, overwrite=TRUE )
        }
      rm( sbox )

      S = matrix( NaN, nrow=nrow(Sloc), ncol=length( p$statsvars ) ) # NA forces into logical
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$S = big.matrix(nrow=nrow(Sloc), ncol=length( p$statsvars ), type="double"  )
          p$bm$S[] = S
          p$ptr$S  = bigmemory::describe( p$bm$S )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$S  = p$cache$S
          bigmemory::as.big.matrix( S, type="double", backingfile=basename(p$bm$S), descriptorfile=basename(p$cache$S), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$S = ff( S, dim=dim(S), file=p$cache$S, overwrite=TRUE )
        }


      Sflag = matrix( 0L, nrow=nrow(Sloc), ncol=1 )  # 0L is the todo flag
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Sflag = big.matrix(nrow=nrow(Sloc), ncol=1, type="double" )
          p$bm$Sflag[] = 0L # TODO flag
          p$ptr$Sflag  = bigmemory::describe( p$bm$Sflag )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Sflag  = p$cache$Sflag
          bigmemory::as.big.matrix( Sflag, type="double", backingfile=basename(p$bm$Sflag), descriptorfile=basename(p$cache$Sflag), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Sflag = ff( Sflag, dim=dim(Sflag), file=p$cache$Sflag, overwrite=TRUE )
        }

      rm(S, Sflag, Sloc)

      # dependent variable
      Yraw = as.matrix(DATA$input[, p$variables$Y ])
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Yraw = big.matrix( nrow=nrow(Yraw), ncol=1, type="double"  )
          p$bm$Yraw[] = Yraw
          p$ptr$Yraw  = bigmemory::describe( p$bm$Yraw )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Yraw  = p$cache$Yraw
          bigmemory::as.big.matrix( Yraw, type="double", backingfile=basename(p$bm$Yraw), descriptorfile=basename(p$cache$Yraw), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Yraw = ff( Yraw, dim=dim(Yraw), file=p$cache$Yraw, overwrite=TRUE )
        }
      rm(Yraw)

      # limits based on quantiles to permit in predictions
      Yraw = lbm_attach( p$storage.backend, p$ptr$Yraw )
      p$qs0 = quantile( Yraw[], probs=p$lbm_quantile_bounds, na.rm=TRUE  )

     # default just copy Yraw ... but if covars are modelled then overwrite with residuals (below)
      Ydata = Yraw[]
      if (exists("lbm_global_modelengine", p)) {
        # to add global covariate model ??  .. simplistic this way but faster
        covmodel = lbm_db( p=p, DS="global_model")
        if (!is.null(covmodel)) {
          Ydata = predict(covmodel, type="response", se.fit=FALSE )
          Ydata = Yraw[] - Ydata
        }
        covmodel =NULL; gc()
      }

      # data to be worked upon .. either the raw data or covariate-residuals
      Ydata = as.matrix( Ydata )
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Y = big.matrix( nrow=nrow(Ydata), ncol=1, type="double"  )
          p$bm$Y[] = Ydata
          p$ptr$Y  = bigmemory::describe( p$bm$Y )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Y  = p$cache$Y
          bigmemory::as.big.matrix( Ydata, type="double", backingfile=basename(p$bm$Y), descriptorfile=basename(p$cache$Y), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Y = ff( Ydata, dim=dim(Ydata), file=p$cache$Y, overwrite=TRUE )
        }
      rm(Ydata)

      Y = lbm_attach( p$storage.backend, p$ptr$Y )
      p$qs = quantile( Y[], probs=p$lbm_quantile_bounds, na.rm=TRUE  )


      if (p$lbm_local_modelengine == "habitat") {
        logitY = lbm_db( p=p, DS="presence.absense" )
          if (p$storage.backend == "bigmemory.ram" ) {
            if (!exists("habitat.threshold.quantile", p)) p$habitat.threshold.quantile = 0.01
            p$bm$Ylogit = big.matrix( nrow=nrow(logitY), ncol=1, type="double" )
            p$bm$Ylogit[] = logitY
            p$ptr$Ylogit  = bigmemory::describe( p$bm$Ylogit )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ylogit  = p$cache$Ylogit
            bigmemory::as.big.matrix( logitY, type="double", backingfile=basename(p$bm$Ylogit), descriptorfile=basename(p$cache$Ylogit), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ylogit = ff( Ylogit, dim=dim(logitY), file=p$cache$Ylogit, overwrite=TRUE )
          }
        rm(logitY)
      }


     # data coordinates
      Yloc = as.matrix( DATA$input[, p$variables$LOCS ])
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Yloc = big.matrix( nrow=nrow(Yloc), ncol=ncol(Yloc), type="double" )
          p$bm$Yloc[] = Yloc
          p$ptr$Yloc = bigmemory::describe( p$bm$Yloc )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Yloc  = p$cache$Yloc
          bigmemory::as.big.matrix( Yloc, type="double", backingfile=basename(p$bm$Yloc), descriptorfile=basename(p$cache$Yloc), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Yloc = ff( Yloc, dim=dim(Yloc), file=p$cache$Yloc, overwrite=TRUE )
        }
      rm(Yloc)


      # independent variables/ covariate
      if (exists("COV", p$variables)) {
        Ycov = as.matrix(  DATA$input[ , p$variables$COV ] )
          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Ycov = big.matrix( nrow=nrow(Ycov), ncol=ncol(Ycov), type="double")
            p$bm$Ycov[] = Ycov
            p$ptr$Ycov  = bigmemory::describe( p$bm$Ycov )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ycov  = p$cache$Ycov
            bigmemory::as.big.matrix( Ycov, type="double", backingfile=basename(p$bm$Ycov), descriptorfile=basename(p$cache$Ycov), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ycov = ff( Ycov, dim=dim(Ycov), file=p$cache$Ycov, overwrite=TRUE )
          }
        rm(Ycov)
      }


      # data times
      if ( exists("TIME", p$variables) ) {
        Ytime = as.matrix(  DATA$input[, p$variables$TIME ] )
          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Ytime = big.matrix( nrow=nrow(Ytime), ncol=ncol(Ytime), type="double"  )
            p$bm$Ytime[] = Ytime
            p$ptr$Ytime  = bigmemory::describe( p$bm$Ytime )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ytime  = p$cache$Ytime
            bigmemory::as.big.matrix( Ytime, type="double", backingfile=basename(p$bm$Ytime), descriptorfile=basename(p$cache$Ytime), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ytime = ff( Ytime, dim=dim(Ytime), file=p$cache$Ytime, overwrite=TRUE )
          }
        Ytime =NULL; gc()
      }


      if (exists("COV", p$variables)) {
        # this needs to be done as Prediction covars need to be structured as lists
        if (!exists("Pcov", p$ptr) ) p$ptr$Pcov = list()
        for ( covname in p$variables$COV ) {
          Pcovdata = as.matrix( DATA$output$COV[[covname]] )
          attr( Pcovdata, "dimnames" ) = NULL
          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Pcov = list()
            p$bm$Pcov[[covname]] = big.matrix( nrow=nrow(Pcovdata), ncol=ncol(Pcovdata), type="double"  )
            p$bm$Pcov[[covname]][] = Pcovdata
            p$ptr$Pcov[[covname]]  = bigmemory::describe( p$bm$Pcov[[covname]] )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Pcov[[covname]]  = p$cache$Pcov[[covname]]
            bigmemory::as.big.matrix( Pcovdata, type="double", backingfile=basename(p$bm$Pcov[[covname]]), descriptorfile=basename(p$cache$Pcov[[covname]]), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Pcov[[covname]] = ff( Pcovdata, dim=dim(Pcovdata), file=p$cache$Pcov[[covname]], overwrite=TRUE )
          }
          Pcovdata = NULL; gc()
        }
      }


      # predictions and associated stats
      P = matrix( NaN, nrow=nrow(DATA$output$LOCS), ncol=p$nt )
        # predictions
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$P = big.matrix( nrow=nrow(P), ncol=ncol(P), type="double" )
          p$bm$P[] = P
          p$ptr$P  = bigmemory::describe( p$bm$P )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$P  = p$cache$P
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$P), descriptorfile=basename(p$cache$P), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$P = ff( P, dim=dim(P), file=p$cache$P, overwrite=TRUE )
        }

      # count of prediction estimates
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Pn = big.matrix( nrow=nrow(P), ncol=ncol(P), type="double" )
          p$bm$Pn[] = P
          p$ptr$Pn = bigmemory::describe( p$bm$Pn )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Pn  = p$cache$Pn
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Pn), descriptorfile=basename(p$cache$Pn), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Pn = ff( P, dim=dim(P), file=p$cache$Pn, overwrite=TRUE )
        }

      # sd of prediction estimates
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Psd = big.matrix( nrow=nrow(P), ncol=ncol(P), type="double" )
          p$bm$Psd[] = P
          p$ptr$Psd =bigmemory::describe( p$bm$Psd )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Psd  = p$cache$Psd
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Psd), descriptorfile=basename(p$cache$Psd), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Psd = ff( P, dim=dim(P), file=p$cache$Psd, overwrite=TRUE )
        }


        if (p$lbm_local_modelengine == "habitat") {
          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Plogit= big.matrix( nrow=nrow(P), ncol=ncol(P) , type="double" )
            p$bm$Plogit[] = P
            p$ptr$Plogit = bigmemory::describe(p$bm$Plogit )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Plogit  = p$cache$Plogit
            bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Plogit), descriptorfile=basename(p$cache$Plogit), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Plogit = ff( P, dim=dim(P), file=p$cache$Plogit, overwrite=TRUE )
          }

          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Plogitsd= big.matrix( nrow=nrow(P), ncol=ncol(P) , type="double" )
            p$bm$Plogitsd[] = P
            p$ptr$Plogitsd = bigmemory::describe(p$bm$Plogitsd )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Plogitsd  = p$cache$Plogitsd
            bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Plogitsd), descriptorfile=basename(p$cache$Plogitsd), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Plogitsd = ff( P, dim=dim(P), file=p$cache$Plogitsd, overwrite=TRUE )
          }
        }

      # prediction coordinates
      Ploc = as.matrix( DATA$output$LOCS )
      attr( Ploc, "dimnames" ) = NULL
         if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Ploc = big.matrix( nrow=nrow(Ploc), ncol=ncol(Ploc), type="double" )
            p$bm$Ploc[] = Ploc
            p$ptr$Ploc  = bigmemory::describe( p$bm$Ploc )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ploc  = p$cache$Ploc
            bigmemory::as.big.matrix( Ploc, type="double", backingfile=basename(p$bm$Ploc), descriptorfile=basename(p$cache$Ploc), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ploc = ff( Ploc, dim=dim(Ploc), file=p$cache$Ploc, overwrite=TRUE )
          }
      Ploc = DATA = NULL; gc()

      if (exists("lbm_global_modelengine", p) ) {
      # create prediction suface with covariate-based additive offsets

        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$P0= big.matrix( nrow=nrow(P), ncol=ncol(P) , type="double" )
          p$bm$P0[] = P
          p$ptr$P0 = bigmemory::describe(p$bm$P0 )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$P0  = p$cache$P0
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$P0), descriptorfile=basename(p$cache$P0), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$P0 = ff( P, dim=dim(P), file=p$cache$P0, overwrite=TRUE )
        }

        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$P0sd= big.matrix( nrow=nrow(P), ncol=ncol(P) , type="double" )
          p$bm$P0sd[] = P
          p$ptr$P0sd = bigmemory::describe(p$bm$P0sd )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$P0sd  = p$cache$P0sd
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$P0sd), descriptorfile=basename(p$cache$P0sd), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$P0sd = ff( P, dim=dim(P), file=p$cache$P0sd, overwrite=TRUE )
        }

        P=NULL; gc()

        # test to see if all covars are static as this can speed up the initial predictions
        message( "Predicting global effect of covariates at each prediction location ... ")
        message( "depending upon the size of the prediction grid and number of cpus (~1hr?)..")
        p$timec0 =  Sys.time()
        nc_cov =NULL
        for (i in p$variables$COV ) {
          pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[i]] )
          nc_cov = c( nc_cov,  ncol(pu) )
        }
        p$all.covars.static = ifelse( any(nc_cov > 1),  FALSE, TRUE )
        if (p$all.covars.static) {
          p = make.list( list( tindex=1:p$nt) , Y=p ) # takes about 28 GB per run .. adjust cluster number temporarily
          lbm_db( p=p, DS="global.prediction.surface" )
        } else {
          if (!exists("no.clusters.covars") ) p$no.clusters.covars = 4
          p$clusters0 = p$clusters
          p$clusters = p$clusters[p$no.clusters.covars]
          p = make.list( list( tindex=1:p$nt) , Y=p ) # takes about 28 GB per run .. adjust cluster number temporarily
          parallel.run( lbm_db, p=p, DS="global.prediction.surface" )
          p$clusters= p$clusters0
        }
        p$timec1 =  Sys.time()
        message( paste( "Time taken to predict covariate surface (mins):", round(difftime( p$timec1, p$timec0 , units="mins"), 3) ) )
    }

    P = NULL; gc() # yes, repeat in case covs are not modelled

  
    lbm_db( p=p, DS="statistics.Sflag" )

    Y = lbm_attach( p$storage.backend, p$ptr$Y )
    Yloc = lbm_attach( p$storage.backend, p$ptr$Yloc )

    Yi = 1:length(Y) # index with useable data
    bad = which( !is.finite( Y[]))
    if (length(bad)> 0 ) Yi[bad] = NA

    # data locations
    bad = which( !is.finite( rowSums(Yloc[])))
    if (length(bad)> 0 ) Yi[bad] = NA

  # data locations
    if (exists("COV", p$variables)) {
      Ycov = lbm_attach( p$storage.backend, p$ptr$Ycov )
      if (length(p$variables$COV)==1) {
        bad = which( !is.finite( Ycov[] ))
      } else {
        bad = which( !is.finite( rowSums(Ycov[])))
      }
      if (length(bad)> 0 ) Yi[bad] = NA
      Yi = na.omit(Yi)
    }

    # data locations
    if (exists("TIME", p$variables)) {
      Ytime = lbm_attach( p$storage.backend, p$ptr$Ytime )
      bad = which( !is.finite( Ytime[] ))
      if (length(bad)> 0 ) Yi[bad] = NA
      Yi = na.omit(Yi)
    }
    bad = NULL

    Yi = as.matrix(Yi)
      if (p$storage.backend == "bigmemory.ram" ) {
        p$bm$Yi = big.matrix( nrow=nrow(Yi), ncol=ncol(Yi), type="double" )
        p$bm$Yi[] = Yi
        p$ptr$Yi  = bigmemory::describe( p$bm$Yi )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Yi  = p$cache$Yi
        bigmemory::as.big.matrix( Yi, type="double", backingfile=basename(p$bm$Yi), descriptorfile=basename(p$cache$Yi), backingpath=p$stloc )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Yi = ff( Yi, dim=dim(Yi), file=p$cache$Yi, overwrite=TRUE )
      }
    rm(Yi)

    if ( !exists("lbm_distance_scale", p)) {
      Yloc = lbm_attach( p$storage.backend, p$ptr$Yloc )
      p$lbm_distance_scale = min( diff(range( Yloc[,1]) ), diff(range( Yloc[,2]) ) ) / 10
      message( paste( "Crude distance scale:", p$lbm_distance_scale ) )
    }
    if ( !exists("lbm_distance_min", p)) p$lbm_distance_min = mean( c(p$lbm_distance_prediction, p$lbm_distance_scale /20 ) )
    if ( !exists("lbm_distance_max", p)) p$lbm_distance_max = mean( c(p$lbm_distance_prediction*10, p$lbm_distance_scale * 2 ) )

    if ( !exists("sampling", p))  {
      # fractions of distance scale  to try in local block search
      p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )
    }

    lbm_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data objects
    message( "Finished. Moving onto analysis... ")
    gc()

    p <<- p  # push to parent in case a manual restart is possible

  } 
  


  # -------------------------------------
  # localized space-time modelling/interpolation/prediction
  message ("Current status is updated to file:" )
  message (p$lbm_current_status )

  p$timei0 =  Sys.time()
  o = lbm_db( p=p, DS="statistics.status" )
  # o = lbm_db(p=p, DS="statistics.status.reset" )
  p = make.list( list( locs=sample( o$todo )) , Y=p ) # random order helps use all cpus
  # lbm_interpolate (p=p )
  parallel.run( lbm_interpolate, p=p )
  p$timei1 =  Sys.time()
  message( paste( "Time taken for main interpolations (mins):", round( difftime( p$timei1, p$timei0, units="mins" ),3) ) )
  gc()


    if (0) {
      Ploc = lbm_attach( p$storage.backend, p$ptr$Ploc )
      P = lbm_attach( p$storage.backend, p$ptr$P )
      j = which( P[] > 5 & P[] < 1000 )

      lattice::levelplot( log(P[j,1])~Ploc[j,1]+Ploc[j,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
      for (i in 1:p$nt) {
        print( lattice::levelplot( P[j,i] ~ Ploc[j,1] + Ploc[j,2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      }
      Sloc = lbm_attach( p$storage.backend, p$ptr$Sloc )
      S = lbm_attach( p$storage.backend, p$ptr$S )
      lattice::levelplot(S[,1]~Sloc[,1]+Sloc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
    }



  # stage 2... revisit eac location in case it was due to a bad subsample 
  toredo = lbm_db( p=p, DS="flag.incomplete.predictions" )
  if ( !is.null(toredo) && length(toredo) > 0) { 
    Sflag = lbm_attach( p$storage.backend, p$ptr$Sflag )
    Sflag[toredo]=0L
    p = make.list( list( locs=sample( toredo )) , Y=p ) # random order helps use all cpus
    parallel.run( lbm_interpolate, p=p )
  }
  p$timei2 =  Sys.time()
  message( paste( "Time taken to stage 2 interpolations (mins):", 
                 round( difftime( p$timei2, p$timei1, units="mins" ),3) ) )


 # stage 3 .. redo, using solutions with relaxed spatial extent
  o = lbm_db(p=p, DS="statistics.status.reset" ) 
  if (length(o$todo) > 0) {
    p$lbm_distance_prediction = p$lbm_distance_prediction * 2
    p$lbm_distance_max = p$lbm_distance_max * 2
    p$lbm_distance_scale = p$lbm_distance_scale*2 # km ... approx guess of 95% AC range 
    p = make.list( list( locs=sample( o$todo )) , Y=p ) # random order helps use all cpus
    parallel.run( lbm_interpolate, p=p )
  }
  p$timei3 =  Sys.time()
  message( paste( "Time taken to stage 3 interpolations (mins):", 
                 round( difftime( p$timei3, p$timei2, units="mins" ),3) ) )

 
  # stage 4 .. last resort .. use a simple but failsafe method to interopalte the rest
  toredo = lbm_db( p=p, DS="flag.incomplete.predictions" )
  if ( !is.null(toredo) && length(toredo) > 0) { 
    Sflag = lbm_attach( p$storage.backend, p$ptr$Sflag )
    Sflag[toredo]=0L
    p = p$lbm_local_modelengine = "fft"  
    p$lbm_fft_filter = "spatial.process"
    p = bio.bathymetry::bathymetry.parameters( p=p, DS="lbm" )
    p = make.list( list( locs=sample( toredo )) , Y=p ) # random order helps use all cpus
    parallel.run( lbm_interpolate, p=p )
  }
  p$timei4 =  Sys.time()

  message( paste( "Time taken to stage 2 interpolations (mins):", 
                 round( difftime( p$timei4, p$timei3, units="mins" ),3) ) )

  # save solutions to disk (again .. overwrite)
  message( "Saving results to disk .. " )
  lbm_db( p=p, DS="lbm.prediction.redo" ) # save to disk for use outside lbm*
  lbm_db( p=p, DS="stats.to.prediction.grid.redo") # save to disk for use outside lbm*

  message ("Finished! \n")
  
  if ( p$storage.backend !="bigmemory.ram" ) {
    resp = readline( "To delete temporary files, type <YES>:  ")
    if (resp=="YES") {
      lbm_db( p=p, DS="cleanup" )
    } else {
      message( "Leaving temporary files alone in case you need to examine them or restart a process.")
      message( "You can delete them by running: lbm_db( p=p, DS='cleanup' ), once you are done.")
    }
  }

  p$time.end =  Sys.time()
  message( paste( "Time taken for full analysis (mins):", round( difftime( p$time.end, p$time.start, units="mins" ),3) ) )

  return( p )

}

