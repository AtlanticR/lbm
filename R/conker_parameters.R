

conker_parameters = function( p=NULL  ) {
  # some generic defaults
  if (is.null(p)) stop( "Parameter list is not structured properly" )

  if (!exists("clusters", p)) {
    message("'p$clusters' was not defined, using localhost with all local compute nodes...")
    p$clusters = rep("localhost", detectCores() )  # default if not given
  }

  if( !exists( "conker_variogram_method", p)) p$conker_variogram_method="fast"   # note GP methods are slow when there is too much data

  if (!exists( "conker_local_family", p)) p$conker_local_family = gaussian()
  if (!exists( "conker_global_family", p)) p$conker_global_family = gaussian()
  
  if (!exists( "conker_noise", p)) p$conker_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
  if (!exists( "conker_quantile_bounds", p)) p$conker_quantile_bounds = c(0.001, 0.999) # remove these extremes in interpolations

  return(p)
}


