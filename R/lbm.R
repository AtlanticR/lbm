

lbm = function( p, DATA,  storage.backend="bigmemory.ram", tasks=c("initiate", "stage1", "save") ) {

  #\\ localized modelling of space and time data to predict/interpolate upon a grid OUT
  #\\ overwrite = FALSE restarts from a saved state
  #\\ speed ratings: bigmemory.ram (1), ff (2), bigmemory.filebacked (3)

  # TODO: splancs::kernel3d as a method ? .. for count data?
  # TODO: habitat methods
  # TODO: gaussian process // see lbm_interpolate_xy_simple
  #       .. the convoSPAT seems almost fast enough
  # TODO: LaplacesDemon method
  # TODO: look at bayesX a little more carefully.
  # TODO: MBA mba.surf method? ... seems very fast


  if ( "continue" %in% tasks) {
    message( "Continuing from an interrupted start" ) 
    p = lbm_db( p=p, DS="load.parameters" )  # ie. restart with saved parameters
    RLibrary( p$libs )
    lbm_db(p=p, DS="statistics.status.reset" )
  } 
    
  if ( "initiate" %in% tasks) {
    lbm_db( p=p, DS="initialize", tasks=tasks )
    lbm_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data objects
    message( "Finished. Moving onto analysis... ")
    p <<- p  # push to parent in case a manual restart is needed
    gc()
  } 
  

  # -------------------------------------
  # localized space-time modelling/interpolation/prediction
  message("Monitor the status of modelling by looking at the output of the following file (e.g., 'watch -n 60 cat {directory}/lbm_current_status'" )
  message (p$lbm_current_status )
  

  if ("stage0" %in% tasks) {
    currentstatus = lbm_db( p=p, DS="statistics.status" )
    print( c( unlist( currentstatus[ c("n.total", "n.land", "n.todo", "n.problematic", "n.outside", "n.complete" ) ] ) ) )
    p = make.list( list( locs=sample( currentstatus$todo )) , Y=p ) # random order helps use all cpus
    lbm_interpolate (p=p )
    browser()
  }

  if ( "stage1" %in% tasks) {  
    p$timei1 =  Sys.time()
    # this is the basic run
    currentstatus = lbm_db( p=p, DS="statistics.status" )
    # currentstatus = lbm_db(p=p, DS="statistics.status.reset" )
    p = make.list( list( locs=sample( currentstatus$todo )) , Y=p ) # random order helps use all cpus
    # lbm_interpolate (p=p )
    parallel.run( lbm_interpolate, p=p )
    p$time_stage1 = round( difftime( Sys.time(), p$timei1, units="hours" ), 3 )
    message("---")
    message( paste( "Time taken for main stage 1, interpolations (hours):", p$time_stage1, "" ) )
    currentstatus = lbm_db( p=p, DS="statistics.status" )
    print( c( unlist( currentstatus[ c("n.total", "n.land", "n.todo", "n.problematic", "n.outside", "n.complete" ) ] ) ) )
    gc()
  }


  if ( "debug_pred_map" %in% tasks) {  
      Ploc = lbm_attach( p$storage.backend, p$ptr$Ploc )
      P = lbm_attach( p$storage.backend, p$ptr$P )
      j = which( P[] > 5 & P[] < 1000 )

      lattice::levelplot( (P[,1])~Ploc[,1]+Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
  }

  if (0) {
      lattice::levelplot( log(P[j,1])~Ploc[j,1]+Ploc[j,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
      for (i in 1:p$nt) {
        print( lattice::levelplot( P[j,i] ~ Ploc[j,1] + Ploc[j,2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      }
  }
  
  if ( "debug_stats_map" %in% tasks) {  
      Sloc = lbm_attach( p$storage.backend, p$ptr$Sloc )
      S = lbm_attach( p$storage.backend, p$ptr$S )
      lattice::levelplot(S[,1]~Sloc[,1]+Sloc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
  }


  if ( "stage2" %in% tasks) {
    p$timei2 =  Sys.time()
    message("---")
    message( "Starting stage 2: more permisssive/relaxed distance settings (spatial extent) " )
    currentstatus = lbm_db(p=p, DS="statistics.status.reset" ) 
    if (length(currentstatus$todo) > 0) {
      p$lbm_distance_prediction = p$lbm_distance_prediction *  p$lbm_multiplier_stage2
      p$lbm_distance_max = p$lbm_distance_max *  p$lbm_multiplier_stage2
      p$lbm_distance_scale = p$lbm_distance_scale* p$lbm_multiplier_stage2 # km ... approx guess of 95% AC range 
      p = make.list( list( locs=sample( currentstatus$todo )) , Y=p ) # random order helps use all cpus
      parallel.run( lbm_interpolate, p=p )
    }
    p$time_stage2 = round( difftime( Sys.time(), p$timei2, units="hours" ), 3)
    message("---")
    message( paste( "Time taken to stage 2 interpolations (hours):", p$time_stage2, "" ) )
    currentstatus = lbm_db( p=p, DS="statistics.status" )
    print( c( unlist( currentstatus[ c("n.total", "n.land", "n.todo", "n.problematic", "n.outside", "n.complete" ) ] ) ) )
    gc()
  }


  if ( "stage3" %in% tasks) {
    p$timei3 =  Sys.time()
    message("---")
    message( "Starting stage 3: simple TPS-based failsafe method to interpolate all the remaining locations " )
    toredo = lbm_db( p=p, DS="flag.incomplete.predictions" )
    if ( !is.null(toredo) && length(toredo) > 0) { 
      Sflag = lbm_attach( p$storage.backend, p$ptr$Sflag )
      Sflag[toredo]=0L
      p$lbm_local_modelengine = "tps"  
      p = bio.bathymetry::bathymetry.parameters( p=p, DS="lbm" )
      p = make.list( list( locs=sample( toredo )) , Y=p ) # random order helps use all cpus
      parallel.run( lbm_interpolate, p=p )
    }
    p$time_stage3 = round( difftime( Sys.time(), p$timei3, units="hours" ), 3)
    message( paste( "Time taken to stage 3 interpolations (hours):", p$time_stage3, "" ) )
    currentstatus = lbm_db( p=p, DS="statistics.status" )
    print( c( unlist( currentstatus[ c("n.total", "n.land", "n.todo", "n.problematic", "n.outside", "n.complete" ) ] ) ) )
    gc()
  }


  # save again, in case some timings/etc needed in a restart
  p <<- p  # push to parent in case a manual restart is possible
  lbm_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data objects


  if ("save" %in% tasks) {
    # save solutions to disk (again .. overwrite)
    message("---")
    message( "Saving predictions to disk .. " )
    lbm_db( p=p, DS="lbm.prediction.redo" ) # save to disk for use outside lbm*
 
    message( "Saving statistics to disk .. " )
    lbm_db( p=p, DS="stats.to.prediction.grid.redo") # save to disk for use outside lbm*

    message ("Finished! ")
  }

  if ( p$storage.backend !="bigmemory.ram" ) {
    resp = readline( "To delete temporary files, type <YES>:  ")
    if (resp=="YES") {
      lbm_db( p=p, DS="cleanup" )
    } else {
      message("---")
      message( "Leaving temporary files alone in case you need to examine them or restart a process. ")
      message( "You can delete them by running: lbm_db( p=p, DS='cleanup' ), once you are done. ")
    }
  }

  p$time_total = round( difftime( Sys.time(), p$time.start, units="hours" ),3)
  message("---")
  message( paste( "Time taken for full analysis (hours):", p$time_total, "" ) )

  return( p )
}

