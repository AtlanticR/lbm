

lstfilter_parameters = function( p=NULL  ) {
  # some generic defaults
  if (is.null(p)) stop( "Parameter list is not structured properly" )

  if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )  # default if not given
  if( !exists( "storage.backend", p))  p$storage.backend="bigmemory.ram"
  if( !exists( "lstfilter_variogram_method", p)) p$lstfilter_variogram_method="fast"   # note GP methods are slow when there is too much data
  if (!exists( "lstfilter_local_family", p)) p$lstfilter_local_family = gaussian()
  if (!exists( "lstfilter_global_family", p)) p$lstfilter_global_family = gaussian()
  if (!exists( "lstfilter_noise", p)) p$lstfilter_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
  if (!exists( "lstfilter_quantile_bounds", p)) p$lstfilter_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
  if (!exists( "eps", p)) p$eps = 1e-6 # floating point precision 
  if (!exists( "boundary", p)) p$boundary = FALSE
  if (!exists( "depth.filter", p)) p$depth.filter = FALSE # depth is given as log(depth) so, choose andy stats locations with elevation > 1 m as being on land
  if (!exists( "lstfilter_kernelmethods_use_all_data", p)) p$lstfilter_kernelmethods_use_all_data =TRUE ## speed and RAM usage improvement is minimal (if any) when off, leave on or remove option and fix as on

  return(p)
}


