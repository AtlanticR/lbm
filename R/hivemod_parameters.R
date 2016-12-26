

hivemod_parameters = function( p=NULL  ) {
  # some generic defaults
  if (is.null(p)) stop( "Parameter list is not structured properly" )
  
  if (!exists("hivemod_current_status", p))  p$hivemod_current_status = "/tmp/hivemod_current_status"

  if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )  # default if not given
  if( !exists( "storage.backend", p))  p$storage.backend="bigmemory.ram"
  if( !exists( "hivemod_variogram_method", p)) p$hivemod_variogram_method="fast"   # note GP methods are slow when there is too much data
  if (!exists( "hivemod_local_family", p)) p$hivemod_local_family = gaussian()
  if (!exists( "hivemod_global_family", p)) p$hivemod_global_family = gaussian()
  if (!exists( "hivemod_noise", p)) p$hivemod_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
  if (!exists( "hivemod_quantile_bounds", p)) p$hivemod_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
  if (!exists( "eps", p)) p$eps = 1e-6 # floating point precision 
  if (!exists( "boundary", p)) p$boundary = FALSE
  if (!exists( "depth.filter", p)) p$depth.filter = FALSE # depth is given as log(depth) so, choose andy stats locations with elevation > 1 m as being on land
  if (!exists( "hivemod_kernelmethods_use_all_data", p)) p$hivemod_kernelmethods_use_all_data =TRUE ## speed and RAM usage improvement is minimal (if any) when off, leave on or remove option and fix as on
  
  # used by "fields" GRMF functions
  if ( p$hivemod_local_modelengine %in% c("gaussianprocess2Dt", "gaussianprocess" )) {
    if (!exists("phi.grid", p) ) p$phi.grid = 10^seq( -6, 6, by=0.5) * p$hivemod_distance_scale # maxdist is aprox magnitude of the phi parameter
    if (!exists("lambda.grid", p) ) p$lambda.grid = 10^seq( -9, 3, by=0.5) # ratio of tau sq to sigma sq
  }

  if ( p$hivemod_local_modelengine %in% c("gam" )) { 
    # p$hivemod_gam_optimizer=c("outer","optim")
    # p$hivemod_gam_optimizer=c("outer","bfgs)
    if (!exists("hivemod_gam_optimizer", p)) p$hivemod_gam_optimizer="perf"
  }

  return(p)
}


