
rspde.inla = function (coords, kappa, variance = 1, alpha = 2, n = 1, verbose = TRUE, mesh, return.attributes = TRUE) {
  
  ## SPDE helper function to create random fields for use with INLA's SPDEs
  ## eg: rspde( matrix(runif(1e6), ncol=2), kappa=2 )
  
  t0 <- Sys.time()
  theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
  if (verbose) cat("theta =", theta, "\n")
  if (missing(mesh)) {
    mesh.pars <- c(0.5, 1, 0.1, 0.5, 1)/kappa
    if (verbose) cat("mesh.pars =", mesh.pars, "\n")
    attributes <- list(mesh = inla.mesh.2d(, coords[chull(coords), ], max.edge = mesh.pars[1:2], cutoff = mesh.pars[3], offset = mesh.pars[4:5]))
    if (verbose) cat("n.mesh =", attributes$mesh$n, "\n")
  } else {
    attributes <- list(mesh = mesh)
  }

  attributes$spde <- inla.spde2.matern(attributes$mesh, alpha = alpha)
  attributes$Q <- inla.spde2.precision(attributes$spde, theta = theta)
  attributes$Q <- inla.spde2.precision(attributes$spde, theta = theta)
  attributes$A <- inla.mesh.project(mesh = attributes$mesh,
  loc = coords)$A
  if (n == 1) result <- drop(attributes$A %*% inla.qsample(Q = attributes$Q, constr = attributes$spde$f$extraconstr))
  t1 <- Sys.time()
  result <- inla.qsample(n, attributes$Q, constr = attributes$spde$f$extraconstr)
  if (nrow(result) < nrow(attributes$A)) {
    result <- rbind(result, matrix(NA, nrow(attributes$A) - nrow(result), ncol(result)))
    dimnames(result)[[1]] <- paste("x", 1:nrow(result), sep = "")
    for (j in 1:ncol(result)) result[, j] <- drop(attributes$A %*% result[1:ncol(attributes$A), j])
  } else {
    for (j in 1:ncol(result)) result[1:nrow(attributes$A), j] <- drop(attributes$A %*% result[, j])
    result <- result[1:nrow(attributes$A), ]
  }
  t2 <- Sys.time()
  attributes$cpu <- c(prep = t1 - t0, sample = t2 - t1, total = t2 - t0)
  if (return.attributes) attributes(result) <- c(attributes(result), attributes)
  return ( drop(result) )
}





