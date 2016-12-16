
hivemod_mesh_inla = function(locs, lengthscale=NULL, max.edge=NULL, bnd.offset=NULL, cutoff=NULL, convex=NULL, resolution=NULL, plotmesh=FALSE ) {
  #\\ create a mesh object for inla-based analysis
  #\\    note: "-" parameter values means to take as a proportion in inla

  require(INLA)

  if (is.null(lengthscale) ) lengthscale = max( diff(range( locs[,1])), diff(range( locs[,2]) )) / 10  # in absence of range estimate take 1/10 of domain size

  if (is.null(max.edge) ) max.edge = c( 0.25, 1)*lengthscale #   # max size of a triange (in, out) proportion of dist.max
  if (is.null(bnd.offset) ) bnd.offset = -c(  0.025,   0.05 )  # how much to extend inside and outside of boundary: proportion of dist.max
  if (is.null(cutoff) ) cutoff = -c(  0.005, 0.05 ) # min distance allowed between points: proportion of dist.max
  if (is.null(convex) ) convex= c( -0.04, -0.08 )  # convex hull radius
  if (is.null(resolution) ) resolution = 125  # discretization resolution for boundary delineation

  MESH = NULL

  MESH =try(
        inla.mesh.2d (
            loc=locs,
            max.edge = max.edge,
            offset = bnd.offset,
            cutoff = cutoff,
            boundary = list(
              inla.nonconvex.hull(locs, convex=convex[1], resolution=resolution ) ,
              inla.nonconvex.hull(locs, convex=convex[2], resolution=resolution ) )
        ), silent=TRUE
      )

  if ( "try-error" %in% class(MESH)) {

  MESH =try(
        inla.mesh.2d (
            loc=locs,
            max.edge = max.edge * 1.25,
            offset = bnd.offset,
            cutoff = cutoff
        ), silent=TRUE
      )

  }

  if ( "try-error" %in% class(MESH)) {

  MESH =try(
        inla.mesh.2d (
            loc=locs,
            max.edge = max.edge * 1.5 ,  # try a lower resolution
            offset = bnd.offset ,
            cutoff = cutoff
        ), silent=TRUE
      )

  }


  if ( "try-error" %in% class(MESH)) MESH = NULL

  if (!is.null( MESH)  & (plotmesh) ) plot( MESH)

  return( MESH )

}


