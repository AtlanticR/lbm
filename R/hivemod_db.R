
  hivemod_db = function( ip=NULL, DS, p, B=NULL, yr=NULL, ret="mean"  ) {
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

      if (exists("hivemod_global_modelengine", p) ) {
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

    if ( DS %in% c( "statistics.status", "statistics.status.reset") ) {
    
      Sflag = hivemod_attach( p$storage.backend, p$ptr$Sflag )
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

      if ( DS=="statistics.reset.problem.locations" ) {
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
        Yloc = hivemod_attach( p$storage.backend, p$ptr$Yloc )
        Sloc = hivemod_attach( p$storage.backend, p$ptr$Sloc )
      
        plot( Yloc[], pch=".", col="grey" ) # data locations
        bnds = try( hivemod_db( p=p, DS="boundary" ) )
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
        hivemod_db( p=p, DS="boundary.redo" ) # ~ 5 min on nfs
      # last set of filters to reduce problem size
        Sflag = hivemod_attach( p$storage.backend, p$ptr$Sflag )
        bnds = try( hivemod_db( p=p, DS="boundary" ) )
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
          z = hivemod_attach( p$storage.backend, p$ptr$Pcov[["z"]] )[]
          Pabove = which( z < p$depth.filter ) # negative = above land
          Pbelow = which( z >= p$depth.filter )

          Ploc = hivemod_attach( p$storage.backend, p$ptr$Ploc )
          Sloc = hivemod_attach( p$storage.backend, p$ptr$Sloc )

          pidA = array_map( "2->1", trunc( cbind(Ploc[Pabove,1]-p$plons[1], Ploc[Pabove,2]-p$plats[1])/p$pres)+1, c(p$nplons, p$nplats) ) 
          pidB = array_map( "2->1", trunc( cbind(Ploc[Pbelow,1]-p$plons[1], Ploc[Pbelow,2]-p$plats[1])/p$pres)+1, c(p$nplons, p$nplats) )

          sid  = array_map( "2->1", trunc( cbind(Sloc[,1]-p$plons[1], Sloc[,2]-p$plats[1])/p$pres)+1, c(p$nplons, p$nplats) ) 
          
          below = which( is.finite( match( sid, pidB ) )) 
          above = which( is.finite( match( sid, pidA ) ))
          
          Sflag = hivemod_attach( p$storage.backend, p$ptr$Sflag )
          if (length(below) > 0 ) Sflag[below] = 0L
          if (length(above) > 0 ) Sflag[above] = 3L

          if (0) {
            Yloc = hivemod_attach( p$storage.backend, p$ptr$Yloc )
            plot( Yloc[], pch=".", col="grey" ) # data locations
            bnds = try( hivemod_db( p=p, DS="boundary" ) )
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
   
    if (DS %in% c( "boundary.redo", "boundary" ) )  {

      fn =  file.path(p$savedir, "boundary.rdata" )
      if (DS=="boundary") {
        boundary = NULL
        if( file.exists(fn)) load( fn)
        return( boundary )
      }

      # data:
      Y = hivemod_attach(  p$storage.backend, p$ptr$Y )
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        Ycov = hivemod_attach(  p$storage.backend, p$ptr$Ycov )
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov[]) )
        } else {
          bad = which( !is.finite( rowSums(Ycov[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
      }

      ii = na.omit(hasdata)
      Yloc = hivemod_attach(  p$storage.backend, p$ptr$Yloc )
      yplon = trunc( ( Yloc[ii,1] - p$plons[1] )/p$pres) + 1
      yplat = trunc( ( Yloc[ii,2] - p$plats[1] )/p$pres) + 1
      uu = unique( array_map( "2->1", cbind(yplon, yplat), c(p$nplons, p$nplats) ) )
      vv = array_map( "1->2", uu, c(p$nplons, p$nplats) )
      
      ww = cbind( (vv[,1] - 1) * p$pres + p$plons[1], (vv[,2] - 1) * p$pres + p$plats[1] )

      if (!exists("hivemod_nonconvexhull_alpha", p)) p$hivemod_nonconvexhull_alpha=20
      boundary=list( polygon = non_convex_hull( ww, alpha=p$hivemod_nonconvexhull_alpha, plot=FALSE ) )
      
      # statistical output locations
      Sloc = hivemod_attach(  p$storage.backend, p$ptr$Sloc )
      boundary$inside.polygon = point.in.polygon( Sloc[,1], Sloc[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
      
      save( boundary, file=fn, compress=TRUE )
      plot( Yloc[ii,], pch=".", col="grey" ) # data locations
      points( Sloc[which(boundary$inside.polygon==1),], pch=".", col="orange" )
      lines( boundary$polygon[] , col="green", pch=2 )
      message( "Check the map of data and boundaries. ")
      message( "If not suitable, set another value for p$hivemod_nonconvexhull_alpha value (radius; distance) ")
      message( "and re-run hivemod() " )
      return( fn )
    }

    # -----
  
    if (DS %in% c("global_model", "global_model.redo") ) {
      
      fn.global_model = file.path( p$savedir, paste( "global_model", p$hivemod_global_modelengine, "rdata", sep=".") )

      if (DS =="global_model") {
        global_model = NULL
        if (file.exists( fn.global_model ))  load(fn.global_model)
        return(global_model)
      }  

      good = which( is.finite (rowSums(B[ , c(p$variables$Y,p$variables$COV) ])) )
      if (length(good)>0) B= B[good,]

      # as a first pass, model the time-independent factors as a user-defined model
      if (p$hivemod_global_modelengine=="gam") {
        global_model = try( 
          gam( formula=p$hivemod_global_modelformula, data=B, optimizer=c("outer","bfgs"), family=p$hivemod_global_family ) ) 
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
      global_model = hivemod_db( p=p, DS="global_model") 
      if (is.null(global_model)) stop("Covariate model not found.")

      P0 = hivemod_attach( p$storage.backend, p$ptr$P0 ) 
      P0sd = hivemod_attach( p$storage.backend, p$ptr$P0sd ) 
      
      for ( iip in ip ) {
        it = p$runs$tindex[iip]
        pa = NULL # construct prediction surface
        for (i in p$variables$COV ) {
          pu = hivemod_attach( p$storage.backend, p$ptr$Pcov[[i]] )
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
        if (p$hivemod_global_modelengine=="gam") {
          Pbaseline = try( predict( global_model, newdata=pa, type="response", se.fit=TRUE ) ) 
          pa = NULL
          gc()
          if (!inherits(Pbaseline, "try-error")) {
            P0[,it] = Pbaseline$fit
            P0sd[,it] = Pbaseline$se.fit
          } 
          Pbaseline = NULL; gc()
        } else if (p$hivemod_global_modelengine=="abc") {
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
    
    if (DS %in% c("hivemod.prediction.redo", "hivemod.prediction") )  {

      if (DS=="hivemod.prediction")  {
        if (! exists("TIME", p$variables)) {
          fn = file.path( p$savedir, paste("hivemod.prediction", ret, "rdata", sep="." ) )
        } else {
          fn = file.path( p$savedir, paste("hivemod.prediction", ret, yr, "rdata", sep="." ) ) 
        }
        if (file.exists(fn) ) load(fn) 
        if (ret=="mean") return (P)
        if (ret=="sd") return( V)
      }

      PP = hivemod_attach( p$storage.backend, p$ptr$P )
      PPsd = hivemod_attach( p$storage.backend, p$ptr$Psd )
      if (exists("hivemod_global_modelengine", p)) {
        P0 = hivemod_attach( p$storage.backend, p$ptr$P0 )
        P0sd = hivemod_attach( p$storage.backend, p$ptr$P0sd )
      }

      if ( exists("TIME", p$variables)) {
        # outputs are on yearly breakdown
        for ( r in 1:p$ny ) {
          y = p$yrs[r]
          fn1 = file.path( p$savedir, paste("hivemod.prediction.mean",  y, "rdata", sep="." ) )
          fn2 = file.path( p$savedir, paste("hivemod.prediction.sd",  y, "rdata", sep="." ) )
          if (exists("nw", p)) {
            col.ranges = (r-1) * p$nw + (1:p$nw) 
            P = PP  [,col.ranges]
            V = PPsd[,col.ranges] # simpleadditive independent errors assumed
          } else {
            P = PP[,r]
            V = PPsd[,r]
          }
          if (exists("hivemod_global_modelengine", p) ) {
            ## maybe add via simulation ? ... 
            P = P + P0[,r] 
            V = sqrt( V^2 + P0sd[,r]^2) # simple additive independent errors assumed
          }
          save( P, file=fn1, compress=T )
          save( V, file=fn2, compress=T )
          print ( paste("Year:", y)  )
        } 
      } else {
          fn1 = file.path( p$savedir, paste("hivemod.prediction.mean",  "rdata", sep="." ) )
          fn2 = file.path( p$savedir, paste("hivemod.prediction.sd", "rdata", sep="." ) )
          P = PP
          V = PPsd
          if (exists("hivemod_global_modelengine", p) ) {
            P = P + P0 
            V = sqrt( V^2 + P0sd^2) # simple additive independent errors assumed
          }
          save( P, file=fn1, compress=T )
          save( V, file=fn2, compress=T )
      }
    }
  
    # ----------------

    if (DS %in% c("stats.to.prediction.grid.redo", "stats.to.prediction.grid") ) {

      fn = file.path( p$savedir, paste( "hivemod.statistics", "rdata", sep=".") )
      if (DS=="stats.to.prediction.grid") {
        stats = NULL
        if (file.exists(fn)) load(fn)
        return(stats)
      }
    
      Ploc = hivemod_attach( p$storage.backend, p$ptr$Ploc )
      S = hivemod_attach( p$storage.backend, p$ptr$S )
      Sloc = hivemod_attach( p$storage.backend, p$ptr$Sloc )
    
      # locations of the new (output) coord system
      locsout = expand.grid( p$plons, p$plats ) # final output grid
      attr( locsout , "out.attrs") = NULL
      names( locsout ) = p$variables$LOCS

      stats = matrix( NaN, ncol=length( p$statsvars ), nrow=nrow( locsout) )  # output data .. ff does not handle NA's .. using NaN for now
      colnames(stats)=p$statsvars

      # map of row, col indices of input data in the new (output) coordinate system
      l2M = cbind( ( Sloc[,1]-p$plons[1])/p$pres + 1, (Sloc[,2]-p$plats[1])/p$pres + 1)
     
      # matrix representation of the output surface
      M = matrix( NA, nrow=p$nplons, ncol=p$nplats) 
      
      for ( i in 1:length( p$statsvars ) ) {
        M = M[] * NA  # init
        M[l2M] = S[,i] # fill with data in correct locations
        
## TODO :: replace with fft method

        Z = smooth.2d( Y=S[,i], x=Sloc[], ncol=p$nplats, nrow=p$nplons, cov.function=stationary.cov, Covariance="Matern", range=p$hivemod_lowpass_phi, nu=p$hivemod_lowpass_nu ) 
        stats[,i] = Z$z
      }
     # lattice::levelplot( stats[,1] ~ locsout[,1]+locsout[,2])
 
      boundary = try( hivemod_db( p=p, DS="boundary" ) )
      if (!is.null(boundary)) {
        if( !("try-error" %in% class(boundary) ) ) {
        inside.polygon = point.in.polygon( locsout[,1], locsout[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
        o = which( inside.polygon == 0 ) # outside boundary
        if (length(o) > 0) stats[o,] = NA         
      }}

      # subset to match to Ploc
      good = match( array_map( "2->1", trunc( cbind(Ploc[,1]-p$plons[1], Ploc[,2]-p$plats[1])/p$pres)+1, c(p$nplons, p$nplats) ), 
                    array_map( "2->1", trunc( cbind(locsout[,1]-p$plons[1], locsout[,2]-p$plats[1])/p$pres)+1, c(p$nplons, p$nplats) ) )

      # Ploc_id = paste( Ploc[,1], Ploc[,2], sep="~" )
      # locsout_id = paste( locsout$plon, locsout$plat, sep="~" )
      # good = match( Ploc_id, locsout_id )

      bad = which( !is.finite(good))
      if (length(bad) > 0 ) good = good[-bad]
      stats = stats[ good, ]
      colnames(stats) = p$statsvars
      rm (good); gc()
     # lattice::levelplot( stats[,1] ~ Ploc[,1]+Ploc[,2])
 
      if ( exists("depth.filter", p) && p$depth.filter ) {
        # stats is now with the same indices as Pcov, Ploc, etc..
        if ( "z" %in% p$variables$COV ){
          depths = hivemod_attach( p$storage.backend, p$ptr$Pcov[["z"]] )
          land = which( depths[] < p$depth.filter )
          if (length(land)>0) stats[land,] = NA 
          rm(land); gc()
        }
      }

      save( stats, file=fn, compress=TRUE )

    }

    #-------------

    if (DS=="presence.absense") {

      Y = hivemod_attach( p$storage.backend, p$ptr$Y )
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


