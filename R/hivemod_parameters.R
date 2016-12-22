

hivemod_parameters = function( p=NULL  ) {
  # some generic defaults
  if (is.null(p)) stop( "Parameter list is not structured properly" )

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
  
  if (!exists("hivemod_fft_missingvalue", p) ) p$hivemod_fft_missingvalue = function(x) median(x, na.rm=TRUE) # function to determine value to use to fill missing data

  return(p)
}


