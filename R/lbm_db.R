
  lbm_db = function( ip=NULL, DS, p, B=NULL, yr=NULL, ret="mean", tasks=NULL  ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates
    #// B is the xyz or xytz data or the function to get the data to work upon

    # --------------------------
    if (!exists("savedir", p)) {
      p$savedir = file.path(p$project.root, "modelled", p$variables$Y, p$spatial.domain )
    }
    
    if (DS %in% "filenames" ) {
      # input data stored as a bigmemory file to permit operations with min memory usage
      # split into separate components to reduce filelocking conflicts
      
      p$cache =list()
      p$cache$Yraw =    file.path( p$stloc, "input.Yraw.cache" ) # raw data
      p$cache$Y =     file.path( p$stloc, "input.Y.cache" ) # residuals of covar model or raw data if none
      p$cache$Yloc =  file.path( p$stloc, "input.Yloc.cache" )
      p$cache$Yi =    file.path( p$stloc, "input.Yi.cache" ) # index of useable data
      
      p$cache$Ylogit = file.path( p$stloc, "Ylogit.cache" )

      p$cache$P =     file.path( p$stloc, "predictions.cache" )
      p$cache$Psd =   file.path( p$stloc, "predictions_sd.cache" )
      p$cache$Pn =    file.path( p$stloc, "predictions_n.cache" )
      
      if ( exists( "COV", p$variables)) {
        p$cache$Ycov =  file.path( p$stloc, "input.Ycov.cache"  )
        p$cache$Pcov =  list()
        for (cov in p$variables$COV) p$cache$Pcov[[cov]] = file.path( p$stloc, paste("predictions_cov", cov, "cache", sep=".") )
      }

      if (exists( "TIME", p$variables)){
        p$cache$Ytime = file.path( p$stloc, "input.Ytime.cache" )
      }

      p$cache$Ploc =  file.path( p$stloc, "predictions_loc.cache" )

      if (exists("lbm_global_modelengine", p) ) {
        p$cache$P0 = file.path( p$stloc, "P0.cache" )
        p$cache$P0sd = file.path( p$stloc, "P0sd.cache" )
      }

      p$cache$Plogit = file.path( p$stloc, "Plogit.cache" )
      p$cache$Plogitsd = file.path( p$stloc, "Plogitsd.cache" )

      p$cache$S =     file.path( p$stloc, "statistics.cache" )
      p$cache$Sloc =  file.path( p$stloc, "statistics_loc.cache" )
      p$cache$Sflag =     file.path( p$stloc, "statistics_flag.cache" )

      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$bm = p$cache
        for ( i in names(p$bm) ) {
          if ( i=="Pcov" ) {
            for (j in p$variables$COV) p$bm$Pcov[[j]] = gsub(".cache$", ".bigmemory", p$bm$Pcov[[j]] )
          } else {
            p$bm[[i]] = gsub(".cache$", ".bigmemory", p$bm[[i]] )
          }
        }
      }

      if (p$storage.backend == "bigmemory.ram" ) {
        p$bm=list() # initial storage of ram objects
      }      
      
      return(p)
    }

    # --------------------------

    if (DS=="save.parameters")  {
      fns = file.path( p$savedir, "p.rdata" )
      save( p, file=fns )
      message( "Saved parameters to file:")
      message( fns )
    }

    if (DS=="load.parameters")  {
      fns = file.path( p$savedir, "p.rdata" )
      if (file.exists( fns)) load( fns )
      return(p)
    }

    # --------------------------
    if (DS %in% "cleanup" ) {
      for (fn in unlist(p$cache) ) if (length(fn)>0) if (file.exists(fn)) file.remove(fn)
      for (fn in unlist(p$bm) ) if (length(fn)>0)  if (file.exists(fn)) file.remove(fn)
      return( "done" )
    }


    if (DS=="initialize") {

      p$savedir = file.path(p$project.root, "modelled", p$variables$Y, p$spatial.domain )
      message( paste( "In case something should go wrong, intermediary outputs will be placed at:", p$savedir ) )
      if ( !file.exists(p$savedir)) dir.create( p$savedir, recursive=TRUE, showWarnings=FALSE )
      
      p$lbm_current_status = file.path( p$savedir, "lbm_current_status" ) 

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

      # prediction times for space.annual methods, treat time as independent timeslices
      if ( !exists("prediction.ts", p)) p$prediction.ts = 1

      # require knowledge of size of stats output before create S, which varies with a given type of analysis
      othervars = c( )
      if (p$lbm_local_modelengine == "habitat") othervars = c( )
      if (exists("TIME", p$variables) )  othervars = c( "ar_timerange", "ar_1" )
      p$statsvars = unique( c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "range", "phi", "nu", othervars ) )

      message("---")
      message( "Initializing temporary storage of data and outputs files... ")
      message( "These are large files (4 to 6 X 5GB), it will take a minute ... ")
      lbm_db( p=p, DS="cleanup" )

      
      if ( "globalmodel" %in% tasks ) {
        if (exists("lbm_global_modelengine", p)) {
          # to add global covariate model ??  .. simplistic this way but faster ~ kriging with external drift
          lbm_db( p=p, DS="global_model.redo", B=DATA$input )
        }
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
          if (p$lbm_global_modelengine=="gam") require(mgcv)
     
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
          p$bm$Pcov = list()
         
          for ( covname in p$variables$COV ) {
            Pcovdata = as.matrix( DATA$output$COV[[covname]] )
            attr( Pcovdata, "dimnames" ) = NULL
            if (p$storage.backend == "bigmemory.ram" ) {
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
          message("---")
          message( "Predicting global effect of covariates at each prediction location ... ")
          message( "depending upon the size of the prediction grid and number of cpus (~1hr?).. ")

          p$timec_covariates_0 =  Sys.time()
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
          p$time_covariates = round(difftime( Sys.time(), p$timec_covariates_0 , units="hours"), 3)
          message( paste( "Time taken to predict covariate surface (hours):", p$time_covariates ) )
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
        message("---")
        message( paste( "Crude distance scale:", p$lbm_distance_scale, "" ) )
      }
      if ( !exists("lbm_distance_min", p)) p$lbm_distance_min = mean( c(p$lbm_distance_prediction, p$lbm_distance_scale /20 ) )
      if ( !exists("lbm_distance_max", p)) p$lbm_distance_max = mean( c(p$lbm_distance_prediction*10, p$lbm_distance_scale * 2 ) )

      if ( !exists("sampling", p))  {
        # fractions of distance scale  to try in local block search
        p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )
      }

      p <<- p  # push to parent in case a manual restart is needed
      
    }



    # -----------------

    if ( DS %in% c( "statistics.status", "statistics.status.reset") ) {
    
      Sflag = lbm_attach( p$storage.backend, p$ptr$Sflag )
      ioutside = which( Sflag[]==2L )
      itodo = which( Sflag[]==0L )       # 0 = TODO
      idone = which( Sflag[]==1L )       # 1 = completed
      iland = which( Sflag[]==3L )       # 3 = land
      iproblems = which( Sflag[] == 9L ) # 9 not completed due to a failed attempt
      out = list(problematic=iproblems, todo=itodo, completed=idone, outside=ioutside, land=iland,
                 n.total=length(Sflag), n.land=length(iland),
                 n.todo=length(itodo), n.problematic=length(iproblems), 
                 n.outside=length(which(is.finite(ioutside))),
                 n.complete=length(idone) )
      out$prop_incomp=out$n.todo / ( out$n.todo + out$n.complete)

      if ( DS=="statistics.status.reset" ) {
        # to reset all rejected locations 
        if (length(which(is.finite(out$problematic))) > 0) {
          Sflag[out$problematic] = 0L  # to reset all the problem flags to todo
          out$problematic=which( Sflag[] == 9L )
          out$n.problematic = 0
          out$todo = which( Sflag[]==0L )
          out$n.todo = length(which( Sflag[]==0L ))
        }
      }

      message( paste("Proportion to do:", round(out$prop_incomp,5), "\n" )) 
      return( out )
  
      if (0) {
        Yloc = lbm_attach( p$storage.backend, p$ptr$Yloc )
        Sloc = lbm_attach( p$storage.backend, p$ptr$Sloc )
      
        plot( Yloc[], pch=".", col="grey" ) # data locations
        bnds = try( lbm_db( p=p, DS="boundary" ) )
        if ( !is.null(bnds)) {
          lines( bnds$polygon[] , col="green", pch=2 )
          points( Sloc[which(bnds$inside.polygon==1),], pch=".", col="orange", cex=5 )
        }
        points( Sloc[which( Sflag[]== 0L),], pch=".", col="blue", cex=5 )
        points( Sloc[which( Sflag[]== 1L),], pch=".", col="purple", cex=5 )
        points( Sloc[which( Sflag[]== 2L),], pch=".", col="red", cex=5 )
        points( Sloc[which( Sflag[]== 3L),], pch=".", col="yellow", cex=5 )
        points( Sloc[which( Sflag[]== 9L),], pch=".", col="magenta", cex=5 )

      }
    }


    #-------------------

    if ( DS %in% c( "statistics.Sflag" ) ) {
      # create location specific flags for analysis, etc..
      if (exists( "boundary", p) && p$boundary) {
        p$timeb0 =  Sys.time()
        message( "Defining boundary polygon for data .. this reduces the number of points to analyse")
        message( "but takes a few minutes to set up ...")
        lbm_db( p=p, DS="boundary.redo" ) # ~ 5 min on nfs
      # last set of filters to reduce problem size
        Sflag = lbm_attach( p$storage.backend, p$ptr$Sflag )
        bnds = try( lbm_db( p=p, DS="boundary" ) )
        if (!is.null(bnds)) {
          if( !("try-error" %in% class(bnds) ) ) {
            outside = which( bnds$inside.polygon == 0 ) # outside boundary
            if (length(outside)>0) Sflag[outside] = 2L
        }}
        bnds = NULL
        p$timeb1 =  Sys.time()
        message( paste( "Time taken to estimate spatial bounds (mins):", round( difftime( p$timeb1, p$timeb0, units="mins" ),3) ) )
      }

      if ( exists("depth.filter", p) && p$depth.filter ) {
        # additionaldepth-based filter:
        # assuming that there is depth information in Pcov, match Sloc's and filter out locations that fall on land
        if ( "z" %in% p$variables$COV ){
          z = lbm_attach( p$storage.backend, p$ptr$Pcov[["z"]] )[]
          Pabove = which( z < p$depth.filter ) # negative = above land
          Pbelow = which( z >= p$depth.filter )

          Ploc = lbm_attach( p$storage.backend, p$ptr$Ploc )
          Sloc = lbm_attach( p$storage.backend, p$ptr$Sloc )

          pidA = array_map( "xy->1", Ploc[Pabove,], gridparams=p$gridparams ) 
          pidB = array_map( "xy->1", Ploc[Pbelow,], gridparams=p$gridparams )
          sid  = array_map( "xy->1", Sloc[], gridparams=p$gridparams )

          below = which( is.finite( match( sid, pidB ) )) 
          above = which( is.finite( match( sid, pidA ) ))
          
          Sflag = lbm_attach( p$storage.backend, p$ptr$Sflag )
          if (length(below) > 0 ) Sflag[below] = 0L
          if (length(above) > 0 ) Sflag[above] = 3L

          if (0) {
            Yloc = lbm_attach( p$storage.backend, p$ptr$Yloc )
            plot( Yloc[], pch=".", col="grey" ) # data locations
            bnds = try( lbm_db( p=p, DS="boundary" ) )
            if (!is.null(bnds)) {
              if ( !("try-error" %in% class(bnds) ) ) {
                points( Sloc[which(bnds$inside.polygon==1),], pch=".", col="orange" )
                lines( bnds$polygon[] , col="green", pch=2 )
              }
            }
            points( Sloc[which( Sflag[]== 0L),], pch=".", col="blue" )
          }
        }
      }

    }

    #---------------------
   
    if (DS== "flag.incomplete.predictions") {
      # statistics locations where estimations need to be redone 
      P = lbm_attach( p$storage.backend, p$ptr$P )
      noP = which( !is.finite( P[]) )
      uP = NULL
      if( length(noP)>0 ) {
        Sloc = lbm_attach( p$storage.backend, p$ptr$Sloc )
       
        Sloc_nplat = ceiling( diff( p$corners$plat) / p$lbm_distance_statsgrid)
        Sloc_nplon = ceiling( diff( p$corners$plon) / p$lbm_distance_statsgrid)

        Ploc = lbm_attach( p$storage.backend, p$ptr$Ploc )
        uS = array_map( "2->1", round( cbind(Sloc[,1]-p$origin[1], Sloc[,2]-p$origin[2])/p$lbm_distance_statsgrid)+1, c(Sloc_nplon, Sloc_nplat) )
        uP = array_map( "2->1", round( cbind(Ploc[noP,1]-p$origin[1], Ploc[noP,2]-p$origin[2])/p$lbm_distance_statsgrid)+1, c(Sloc_nplon, Sloc_nplat) ) 
        inrange = which( (uP >= min(uS)) & (uP <= max(uS)) )
        if (length( inrange) > 0) uP = uP[inrange] 
        uP = unique(uP)
      }
      return(uP) 
    }

    #---------------------
   
    if (DS %in% c( "boundary.redo", "boundary" ) )  {

      fn =  file.path(p$savedir, "boundary.rdata" )
      if (DS=="boundary") {
        boundary = NULL
        if( file.exists(fn)) load( fn)
        return( boundary )
      }

      # data:
      Y = lbm_attach(  p$storage.backend, p$ptr$Y )
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        Ycov = lbm_attach(  p$storage.backend, p$ptr$Ycov )
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov[]) )
        } else {
          bad = which( !is.finite( rowSums(Ycov[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
      }

      ii = na.omit(hasdata)
      Yloc = lbm_attach(  p$storage.backend, p$ptr$Yloc )
      yplon = round( ( Yloc[ii,1] - p$origin[1] )/p$pres) + 1
      yplat = round( ( Yloc[ii,2] - p$origin[2] )/p$pres) + 1
      uu = unique( array_map( "2->1", cbind(yplon, yplat), c(p$nplons, p$nplats) ) )
      vv = array_map( "1->2", uu, c(p$nplons, p$nplats) )
      
      ww = cbind( (vv[,1] - 1) * p$pres + p$origin[1], (vv[,2] - 1) * p$pres + p$origin[2] )

      if (!exists("lbm_nonconvexhull_alpha", p)) p$lbm_nonconvexhull_alpha=20
      boundary=list( polygon = non_convex_hull( ww, alpha=p$lbm_nonconvexhull_alpha, plot=FALSE ) )
      
      # statistical output locations
      Sloc = lbm_attach(  p$storage.backend, p$ptr$Sloc )
      boundary$inside.polygon = point.in.polygon( Sloc[,1], Sloc[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
      
      save( boundary, file=fn, compress=TRUE )
      plot( Yloc[ii,], pch=".", col="grey" ) # data locations
      points( Sloc[which(boundary$inside.polygon==1),], pch=".", col="orange" )
      lines( boundary$polygon[] , col="green", pch=2 )
      message( "Check the map of data and boundaries. ")
      message( "If not suitable, set another value for p$lbm_nonconvexhull_alpha value (radius; distance) ")
      message( "and re-run lbm() " )
      return( fn )
    }

    # -----
  
    if (DS %in% c("global_model", "global_model.redo") ) {

      fn.global_model = file.path( p$savedir, paste( "global_model", p$lbm_global_modelengine, "rdata", sep=".") )

      if (DS =="global_model") {
        global_model = NULL
        if (file.exists( fn.global_model ))  load(fn.global_model)
        return(global_model)
      }  

      good = which( is.finite (rowSums(B[ , c(p$variables$Y,p$variables$COV) ])) )
      if (length(good)>0) B= B[good,]

      # as a first pass, model the time-independent factors as a user-defined model
      if (p$lbm_global_modelengine=="gam") {
        require(mgcv)
        
        global_model = try( 
          gam( formula=p$lbm_global_modelformula, data=B, optimizer=c("outer","bfgs"), family=p$lbm_global_family ) ) 
      }

      if ( "try-error" %in% class(global_model) ) stop( "The covariate model was problematic" )
      print( summary( global_model ) )
      save( global_model, file= fn.global_model, compress=TRUE )

      return (fn.global_model)
    }

    # -----

    if (DS %in% c("global.prediction.surface") ) {
      if (exists( "libs", p)) RLibrary( p$libs )
      if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns
      global_model = lbm_db( p=p, DS="global_model") 
      if (is.null(global_model)) stop("Covariate model not found.")

      P0 = lbm_attach( p$storage.backend, p$ptr$P0 ) 
      P0sd = lbm_attach( p$storage.backend, p$ptr$P0sd ) 
      
      for ( iip in ip ) {
        it = p$runs$tindex[iip]
        pa = NULL # construct prediction surface
        for (i in p$variables$COV ) {
          pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[i]] )
          ncpu = ncol(pu)
          if ( ncpu== 1 ) {
            pa = cbind( pa, pu[] ) # ie. a static variable (space)
          } else if( ncpu == p$ny )  {
            iy = round( (it-1) / p$nw ) + 1
            pa = cbind( pa, pu[,iy] ) # ie., annual data (space.annual)
          } else if ( ncpu == p$nt) {
            pa = cbind( pa, pu[,it] ) # ie. same time dimension as predictive data (space.annual.seasonal)
          }
        }
        pa = as.data.frame( pa )
        names(pa) = p$variables$COV
        
        if ( any( p$variables$LOCS %in%  all.vars( p$lbm_global_modelformula ) ) ) {
          Ploc = lbm_attach( p$storage.backend, p$ptr$Ploc )
          pa = cbind(pa, Ploc[])
          names(pa) = c( p$variables$COV, p$variables$LOCS )
        }

        if (p$lbm_global_modelengine=="gam") {
          Pbaseline = try( predict( global_model, newdata=pa, type="response", se.fit=TRUE ) ) 
          pa = NULL
          gc()
          if (!inherits(Pbaseline, "try-error")) {
            P0[,it] = Pbaseline$fit
            P0sd[,it] = Pbaseline$se.fit
          } 
          Pbaseline = NULL; gc()
        } else  {
          stop ("This global model method requires a bit more work .. ")
        }
        
        if (p$all.covars.static) {
          # if this is true then this is a single cpu run and all predictions for each time slice is the same
          # could probably catch this and keep storage small but that would make the update math a little more complex
          # this keeps it simple with a quick copy
          if (p$nt  > 1 ) {
            for (j in ip[2:p$nruns]){
              P0[,j] = P0[,1]
              P0sd[,j] = P0sd[,1]
            } 
          }
          global_model =NULL
          return(NULL) 
        }
      } # end each timeslice
      global_model =NULL
      message( "Done ... moving onto the rest of the analysis...")
    }


    # -----
    
    if (DS %in% c("lbm.prediction.redo", "lbm.prediction") )  {

      if (DS=="lbm.prediction") {
        if (! exists("TIME", p$variables)) {
          fn = file.path( p$savedir, paste("lbm.prediction",  ret, "rdata", sep="." ) )
        } else {
          fn = file.path( p$savedir, paste("lbm.prediction",  ret, yr, "rdata", sep="." ) ) 
        }
        if (file.exists(fn) ) load(fn) 
        if (ret=="mean") return (P)
        if (ret=="sd") return( V)
      }

      PP = lbm_attach( p$storage.backend, p$ptr$P )
      PPsd = lbm_attach( p$storage.backend, p$ptr$Psd )
      if (exists("lbm_global_modelengine", p)) {
        P0 = lbm_attach( p$storage.backend, p$ptr$P0 )
        P0sd = lbm_attach( p$storage.backend, p$ptr$P0sd )
      }

      if ( exists("TIME", p$variables)) {
        # outputs are on yearly breakdown
        for ( r in 1:p$ny ) {
          y = p$yrs[r]
          fn1 = file.path( p$savedir, paste("lbm.prediction", "mean", y, "rdata", sep="." ) )
          fn2 = file.path( p$savedir, paste("lbm.prediction", "sd",   y, "rdata", sep="." ) )
          if (exists("nw", p)) {
            col.ranges = (r-1) * p$nw + (1:p$nw) 
            P = PP  [,col.ranges]
            V = PPsd[,col.ranges] # simpleadditive independent errors assumed
          } else {
            P = PP[,r]
            V = PPsd[,r]
          }
          if (exists("lbm_global_modelengine", p) ) {
            ## maybe add via simulation ? ... 
            uu = which(!is.finite(P[]))
            if (length(uu)>0) P[uu] = 0 # permit covariate-base predictions to pass through .. 
            P = P[] + P0[,r] 
            vv = which(!is.finite(V[]))
            if (length(vv)>0) V[vv] = 0 # permit covariate-base predictions to pass through ..
            V = sqrt( V[]^2 + P0sd[,r]^2) # simple additive independent errors assumed
          }
          save( P, file=fn1, compress=T )
          save( V, file=fn2, compress=T )
          print ( paste("Year:", y)  )
        } 
      } else {
          fn1 = file.path( p$savedir, paste("lbm.prediction", "mean", "rdata", sep="." ) )
          fn2 = file.path( p$savedir, paste("lbm.prediction", "sd",   "rdata", sep="." ) )
          P = PP[]
          V = PPsd[]
          if (exists("lbm_global_modelengine", p) ) {
            uu = which(!is.finite(P[]))
            if (length(uu)>0) P[uu] = 0 # permit covariate-base predictions to pass through ..
            P = P[] + P0[] 
            vv = which(!is.finite(V[]))
            if (length(vv)>0) V[vv] = 0 # permit covariate-base predictions to pass through ..
            V = sqrt( V[]^2 + P0sd[]^2) # simple additive independent errors assumed
          }
          save( P, file=fn1, compress=T )
          save( V, file=fn2, compress=T )
      }
    }
    
    if(0) {
      i = 100
      Z = smooth.2d( Y=P[,i], x=Ploc[], ncol=p$nplats, nrow=p$nplons, cov.function=stationary.cov, Covariance="Matern", range=p$lbm_lowpass_phi, nu=p$lbm_lowpass_nu )
      lattice::levelplot( P[,i] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
    }
    # 

    # ----------------

    if (DS %in% c("stats.to.prediction.grid.redo", "stats.to.prediction.grid") ) {
      
      # TODO:: parallelize this

      fn = file.path( p$savedir, paste( "lbm.statistics", "rdata", sep=".") )
      if (DS=="stats.to.prediction.grid") {
        stats = NULL
        if (file.exists(fn)) load(fn)
        return(stats)
      }
    
      Ploc = lbm_attach( p$storage.backend, p$ptr$Ploc )
      S = lbm_attach( p$storage.backend, p$ptr$S )
      Sloc = lbm_attach( p$storage.backend, p$ptr$Sloc )

      Sloc_nplat = ceiling( diff( p$corners$plat) / p$lbm_distance_statsgrid)
      Sloc_nplon = ceiling( diff( p$corners$plon) / p$lbm_distance_statsgrid)

      stats = matrix( NaN, ncol=length( p$statsvars ), nrow=nrow( Ploc) )  # output data .. ff does not handle NA's .. using NaN for now
      colnames(stats)=p$statsvars

      for ( i in 1:length( p$statsvars ) ) {
        print(i)
        # linear interpolation
        u = as.image( S[,i], x=Sloc[,], na.rm=TRUE, nx=Sloc_nplon, ny=Sloc_nplat )
        stats[,i] = as.vector( fields::interp.surface( u, loc=Ploc[] ) ) # linear interpolation
      }

      # lattice::levelplot( stats[,1] ~ Ploc[,1]+Ploc[,2])
 
      boundary = try( lbm_db( p=p, DS="boundary" ) )
      if (!is.null(boundary)) {
        if( !("try-error" %in% class(boundary) ) ) {
        inside.polygon = point.in.polygon( Ploc[,1], Ploc[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
        o = which( inside.polygon == 0 ) # outside boundary
        if (length(o) > 0) stats[o,] = NA         
      }}

    
      if (0){
        i = 1
        ii = which (is.finite(stats[,i]))
        lattice::levelplot( stats[ii,i] ~ Ploc[ii,1]+Ploc[ii,2])
      }
 
      if ( exists("depth.filter", p) && p$depth.filter ) {
        # stats is now with the same indices as Pcov, Ploc, etc..
        if ( "z" %in% p$variables$COV ){
          depths = lbm_attach( p$storage.backend, p$ptr$Pcov[["z"]] )
          land = which( depths[] < p$depth.filter )
          if (length(land)>0) stats[land,] = NA 
          rm(land); gc()
        }
      }

      save( stats, file=fn, compress=TRUE )

    }

    #-------------

    if (DS=="presence.absense") {

      Y = lbm_attach( p$storage.backend, p$ptr$Y )
      z = which( Y == 0) # assumed to be real zeros
      i = which( Y >  0)  # positive values
      
      # determine quantiles .. 
      Yq = rep( 0, length(Y) )
      Yq[z] = 1
    
      pr = ecdf(Y[i])( Y[i] )
      ix = which( pr ==1 )
      if ( !( length(ix) %in% c(0, length(x)) ))  pr[ix] = max( pr[-ix] )
      Yq[i] = pr

      s01 = which( Yq < p$habitat.threshold.quantile )  # buffer zone
      s0 = unique( c(s01, sz ) )
      s1 = which( Yq >= p$habitat.threshold.quantile )

      # determine presence-absence
      Ybin =  rep( NA, length(Y) )
      Ybin[s1] = 1  
      Ybin[s0] = 0
      Ybin[z] = 0

      return(Ybin)
    }

  }


