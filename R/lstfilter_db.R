
  lstfilter_db = function( ip=NULL, DS, p, B=NULL, yr=NULL, ret="mean"  ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates
    #// B is the xyz or xytz data or the function to get the data to work upon

    # --------------------------

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

      if (exists("lstfilter_global_modelengine", p) ) {
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
      save(p, file=file.path( p$savedir, "p.rdata") )
      message( "Saved parameters:")
      message( file.path( p$savedir, "p.rdata") )
    }

    if (DS=="load.parameters")  {
      load( file.path( p$savedir, "p.rdata") )
      return(p)
    }

    # --------------------------
    if (DS %in% "cleanup" ) {
      for (fn in unlist(p$cache) ) if (length(fn)>0) if (file.exists(fn)) file.remove(fn)
      for (fn in unlist(p$bm) ) if (length(fn)>0)  if (file.exists(fn)) file.remove(fn)
      return( "done" )
    }

    # -----------------
    
    if ( DS=="statistics.status" ) {
      # find locations for statistic computation and trim area based on availability of data
      # stats:
      bnds = try( lstfilter_db( p=p, DS="boundary" ) )
      ioutside = NA
      if (!is.null(bnds)) {
        if( !("try-error" %in% class(bnds) ) ) {
          ioutside = which( bnds$inside.polygon == 0 ) # outside boundary
      }}
      Sflag = lstfilter_attach( p$storage.backend, p$ptr$Sflag )
      itodo = setdiff( which( is.nan( Sflag[] )), ioutside)       # incomplete
      idone = setdiff( which( is.finite (Sflag[] )  ), ioutside)      # completed
      iskipped = which( is.infinite( Sflag[] )  ) # skipped due to problems or out of bounds
      iproblems = setdiff( iskipped, ioutside)    # not completed due to a failed attempt
      out = list(problematic=iproblems, skipped=iskipped, todo=itodo, completed=idone, outside=ioutside,
                 n.total=length(Sflag) , n.skipped=length(iskipped),
                 n.todo=length(itodo), n.problematic=length(iproblems), 
                 n.outside=length(which(is.finite(ioutside))),
                 n.complete=length(idone) )
      out$prop_incomp=out$n.todo / ( out$n.todo + out$n.complete)
      message( paste("Proportion to do:", round(out$prop_incomp,5), "\n" )) 
      return( out )
    }
    
  # -------------------

    if ( DS=="statistics.reset.problem.locations" ) {
        # to reset all rejected locations 
        Sflag = lstfilter_attach( p$storage.backend, p$ptr$Sflag )
        o = lstfilter_db( p=p, DS="statistics.status" )
        if (length(which(is.finite(o$skipped))) > 0) Sflag[o$skipped] = NaN  # to reset all the flags
        if (length(which(is.finite(o$outside))) > 0) Sflag[o$outside] = Inf  # flag area outside of data boundary to skip
      }

    #-------------------

    if ( DS %in% c( "statistics.Sflag" ) ) {
    
      if (exists( "boundary", p) && p$boundary) {
        p$timeb0 =  Sys.time()
        message( "Defining boundary polygon for data .. this reduces the number of points to analyse")
        message( "but takes a few minutes to set up ...")
        lstfilter_db( p=p, DS="boundary.redo" ) # ~ 5 min on nfs
      # last set of filters to reduce problem size
        Sflag = lstfilter_attach( p$storage.backend, p$ptr$Sflag )
        bnds = try( lstfilter_db( p=p, DS="boundary" ) )
        if (!is.null(bnds)) {
          if( !("try-error" %in% class(bnds) ) ) {
            to.ignore = which( bnds$inside.polygon == 0 ) # outside boundary
            if (length(to.ignore)>0) Sflag[to.ignore] = Inf
        }}
        bnds = NULL
        p$timeb1 =  Sys.time()
        message( paste( "Time taken to estimate spatial bounds (mins):", round( difftime( p$timeb1, p$timeb0, units="mins" ),3) ) )
      }

      if ( !is.null(p$depth.filter) ) {
        # assuming that there is depth information in Pcov, match Sloc's and filter out locations that fall on land
        if ( "z" %in% p$variables$COV ){
          Pland = which( lstfilter_attach( p$storage.backend, p$ptr$Pcov[["z"]] )[] < p$depth.filter )
          
          Ploc = lstfilter_attach( p$storage.backend, p$ptr$Ploc )
          Sloc = lstfilter_attach( p$storage.backend, p$ptr$Sloc )
          
          Swater = match( array_map( "2->1", cbind(Ploc[-Pland,1]-p$plons[1], Ploc[-Pland,2]-p$plats[1])/p$pres+1, c(p$nplons, p$nplats) ), 
                           array_map( "2->1", cbind(Sloc[,1]-p$plons[1], Sloc[,2]-p$plats[1])/p$pres+1, c(p$nplons, p$nplats) ) )
          Swater = unique( Swater )
          Swater = Swater[ is.finite(Swater)]

          Sland = setdiff( 1:nrow(Sloc), Swater )
          Sflag = lstfilter_attach( p$storage.backend, p$ptr$Sflag )
          if (length(Sland)>0) Sflag[Sland] = Inf
          if (length(Sland)>0) Sflag[Swater] = NaN
           
          rm(land, S_index); gc()
        }
      }

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
      Y = lstfilter_attach(  p$storage.backend, p$ptr$Y )
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        Ycov = lstfilter_attach(  p$storage.backend, p$ptr$Ycov )
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov[]) )
        } else {
          bad = which( !is.finite( rowSums(Ycov[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
      }

      ii = na.omit(hasdata)
      Yloc = lstfilter_attach(  p$storage.backend, p$ptr$Yloc )
      yplon = (grid.internal( Yloc[ii,1], p$plons ) - p$plons[1] )/p$pres + 1
      yplat = (grid.internal( Yloc[ii,2], p$plats ) - p$plats[1] )/p$pres + 1
      uu = unique( array_map( "2->1", cbind(yplon, yplat), c(p$nplons, p$nplats) ) )
      vv = array_map( "1->2", uu, c(p$nplons, p$nplats) )
      
      ww = cbind( (vv[,1] - 1) * p$pres + p$plons[1], (vv[,2] - 1) * p$pres + p$plats[1] )

      if (!exists("lstfilter_nonconvexhull_alpha", p)) p$lstfilter_nonconvexhull_alpha=20
      boundary=list( polygon = non_convex_hull( ww, alpha=p$lstfilter_nonconvexhull_alpha, plot=FALSE ) )
      
      # statistical output locations
      Sloc = lstfilter_attach(  p$storage.backend, p$ptr$Sloc )
      boundary$inside.polygon = point.in.polygon( Sloc[,1], Sloc[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
      
      save( boundary, file=fn, compress=TRUE )
      plot( Yloc[ii,], pch=".", col="grey" ) # data locations
      points( Sloc[which(boundary$inside.polygon==1),], pch=".", col="orange" )
      lines( boundary$polygon[] , col="green", pch=2 )
      message( "Check the map of data and boundaries. ")
      message( "If not suitable, set another value for p$lstfilter_nonconvexhull_alpha value (radius; distance) ")
      message( "and re-run lstfilter() " )
      return( fn )
    }

    # -----
  
    if (DS %in% c("global_model", "global_model.redo") ) {
      
      fn.global_model = file.path( p$savedir, paste( "global_model", p$lstfilter_global_modelengine, "rdata", sep=".") )

      if (DS =="global_model") {
        global_model = NULL
        if (file.exists( fn.global_model ))  load(fn.global_model)
        return(global_model)
      }  

      good = which( is.finite (rowSums(B[ , c(p$variables$Y,p$variables$COV) ])) )
      if (length(good)>0) B= B[good,]

      # as a first pass, model the time-independent factors as a user-defined model
      if (p$lstfilter_global_modelengine=="gam") {
        global_model = try( 
          gam( formula=p$lstfilter_global_modelformula, data=B, optimizer=c("outer","bfgs"), family=p$lstfilter_global_family ) ) 
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
      global_model = lstfilter_db( p=p, DS="global_model") 
      if (is.null(global_model)) stop("Covariate model not found.")

      P0 = lstfilter_attach( p$storage.backend, p$ptr$P0 ) 
      P0sd = lstfilter_attach( p$storage.backend, p$ptr$P0sd ) 
      
      for ( iip in ip ) {
        it = p$runs$tindex[iip]
        pa = NULL # construct prediction surface
        for (i in p$variables$COV ) {
          pu = lstfilter_attach( p$storage.backend, p$ptr$Pcov[[i]] )
          ncpu = ncol(pu)
          if ( ncpu== 1 ) {
            pa = cbind( pa, pu[] ) # ie. a static variable
          } else if( ncpu == p$ny )  {
            iy = trunc( (it-1) / p$nw ) + 1
            pa = cbind( pa, pu[,iy] ) # ie., annual data 
          } else if ( ncpu == p$nt) {
            pa = cbind( pa, pu[,it] ) # ie. same time dimension as predictive data
          }
        }
        pa = as.data.frame( pa )
        names(pa) = p$variables$COV
        if (p$lstfilter_global_modelengine=="gam") {
          Pbaseline = try( predict( global_model, newdata=pa, type="response", se.fit=TRUE ) ) 
          pa = NULL
          gc()
          if (!inherits(Pbaseline, "try-error")) {
            P0[,it] = Pbaseline$fit
            P0sd[,it] = Pbaseline$se.fit
          } 
          Pbaseline = NULL; gc()
        } else if (p$lstfilter_global_modelengine=="abc") {
          # ... other methods ..
        }
        
        if (p$all.covars.static) {
          # if this is true then this is a single cpu run and all predictions for each time slice is the same
          # could probably catch this and keep storage small but that would make the update math a little more complex
          # this keeps it simple with a quick copy
          for (j in ip[2:p$nruns]){
            P0[,j] = P0[,1]
            P0sd[,j] = P0sd[,1]
          } 
          global_model =NULL
          return(NULL) 
        }
      } # end each timeslice
      global_model =NULL
      message( "Done ... moving onto the rest of the analysis...")
    }


    # -----
    
    if (DS %in% c("lstfilter.prediction.redo", "lstfilter.prediction") )  {

      if (DS=="lstfilter.prediction")  {
        if (! exists("TIME", p$variables)) {
          fn = file.path( p$savedir, paste("lstfilter.prediction", ret, "rdata", sep="." ) )
        } else {
          fn = file.path( p$savedir, paste("lstfilter.prediction", ret, yr, "rdata", sep="." ) ) 
        }
        if (file.exists(fn) ) load(fn) 
        if (ret=="mean") return (P)
        if (ret=="sd") return( V)
      }

      PP = lstfilter_attach( p$storage.backend, p$ptr$P )
      PPsd = lstfilter_attach( p$storage.backend, p$ptr$Psd )
      if (exists("lstfilter_global_modelengine", p)) {
        P0 = lstfilter_attach( p$storage.backend, p$ptr$P0 )
        P0sd = lstfilter_attach( p$storage.backend, p$ptr$P0sd )
      }

      if ( exists("TIME", p$variables)) {
        # outputs are on yearly breakdown
        for ( r in 1:p$ny ) {
          y = p$yrs[r]
          fn1 = file.path( p$savedir, paste("lstfilter.prediction.mean",  y, "rdata", sep="." ) )
          fn2 = file.path( p$savedir, paste("lstfilter.prediction.sd",  y, "rdata", sep="." ) )
          if (exists("nw", p)) {
            col.ranges = (r-1) * p$nw + (1:p$nw) 
            P = PP  [,col.ranges]
            V = PPsd[,col.ranges] # simpleadditive independent errors assumed
          } else {
            P = PP[,r]
            V = PPsd[,r]
          }
          if (exists("lstfilter_global_modelengine", p) ) {
            ## maybe add via simulation ? ... 
            P = P + P0[,r] 
            V = sqrt( V^2 + P0sd[,r]^2) # simple additive independent errors assumed
          }
          save( P, file=fn1, compress=T )
          save( V, file=fn2, compress=T )
          print ( paste("Year:", y)  )
        } 
      } else {
          fn1 = file.path( p$savedir, paste("lstfilter.prediction.mean",  "rdata", sep="." ) )
          fn2 = file.path( p$savedir, paste("lstfilter.prediction.sd", "rdata", sep="." ) )
          P = PP
          V = PPsd
          if (exists("lstfilter_global_modelengine", p) ) {
            P = P + P0 
            V = sqrt( V^2 + P0sd^2) # simple additive independent errors assumed
          }
          save( P, file=fn1, compress=T )
          save( V, file=fn2, compress=T )
      }
    }
  
    # ----------------

    if (DS %in% c("stats.to.prediction.grid.redo", "stats.to.prediction.grid") ) {

      fn = file.path( p$savedir, paste( "lstfilter.statistics", "rdata", sep=".") )
      if (DS=="stats.to.prediction.grid") {
        stats = NULL
        if (file.exists(fn)) load(fn)
        return(stats)
      }
    
      Ploc = lstfilter_attach( p$storage.backend, p$ptr$Ploc )
      S = lstfilter_attach( p$storage.backend, p$ptr$S )
      Sloc = lstfilter_attach( p$storage.backend, p$ptr$Sloc )
    
      # locations of the new (output) coord system
      locsout = expand.grid( p$plons, p$plats ) # final output grid
      attr( locsout , "out.attrs") = NULL
      names( locsout ) = p$variables$LOCS

      stats = matrix( NaN, ncol=length( p$statsvars ), nrow=nrow( locsout) )  # output data
      colnames(stats)=p$statsvars

      # map of row, col indices of input data in the new (output) coordinate system
      l2M = cbind( ( Sloc[,1]-p$plons[1])/p$pres + 1, (Sloc[,2]-p$plats[1])/p$pres + 1)
     
      # matrix representation of the output surface
      M = matrix( NA, nrow=p$nplons, ncol=p$nplats) 
      
      for ( i in 1:length( p$statsvars ) ) {
        M = M[] * NA  # init
        M[l2M] = S[,i] # fill with data in correct locations
        Z = smooth.2d( Y=S[,i], x=Sloc[], ncol=p$nplats, nrow=p$nplons, cov.function=stationary.cov, Covariance="Matern", range=p$lstfilter_phi, nu=p$lstfilter_nu ) 
        stats[,i] = Z$z
      }
     # lattice::levelplot( stats[,1] ~ locsout[,1]+locsout[,2])
 
      boundary = try( lstfilter_db( p=p, DS="boundary" ) )
      if (!is.null(boundary)) {
        if( !("try-error" %in% class(boundary) ) ) {
        inside.polygon = point.in.polygon( locsout[,1], locsout[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
        o = which( inside.polygon == 0 ) # outside boundary
        if (length(o) > 0) stats[o,] = NA         
      }}

      # subset to match to Ploc
      good = match( array_map( "2->1", cbind(Ploc[,1]-p$plons[1], Ploc[,2]-p$plats[1])/p$pres+1, c(p$nplons, p$nplats) ), 
                    array_map( "2->1", cbind(locsout[,1]-p$plons[1], locsout[,2]-p$plats[1])/p$pres+1, c(p$nplons, p$nplats) ) )

      # Ploc_id = paste( Ploc[,1], Ploc[,2], sep="~" )
      # locsout_id = paste( locsout$plon, locsout$plat, sep="~" )
      # good = match( Ploc_id, locsout_id )

      bad = which( !is.finite(good))
      if (length(bad) > 0 ) good = good[-bad]
      stats = stats[ good, ]
      colnames(stats) = p$statsvars
      rm (good); gc()
     # lattice::levelplot( stats[,1] ~ Ploc[,1]+Ploc[,2])
 
      if ( !is.null(p$depth.filter) ) {
        # stats is now with the same indices as Pcov, Ploc, etc..
        if ( "z" %in% p$variables$COV ){
          depths = lstfilter_attach( p$storage.backend, p$ptr$Pcov[["z"]] )
          land = which( depths[] < p$depth.filter )
          if (length(land)>0) stats[land,] = NA 
          rm(land); gc()
        }
      }

      save( stats, file=fn, compress=TRUE )

    }

    #-------------

    if (DS=="presence.absense") {

      Y = lstfilter_attach( p$storage.backend, p$ptr$Y )
      z = which( Y == 0) # assumed to be real zeros
      i = which( Y >  0)  # positive values
      
      # determine quantiles
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


