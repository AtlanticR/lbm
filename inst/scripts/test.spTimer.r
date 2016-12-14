
  # spTimer test

if (0) {
  devtools::install_github("kaskr/adcomp/TMB")
  install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/testing")
  install.packages( "RandomFields" )
  install.packages( "RANN" )
}

set.seed(1)


# sim dimensions
nx = 200
ny = 100
nt = 10
nsample = 200
si = sample( 1:(nx*ny), nsample )

locs <- list( x= seq(1, nx), y=seq(1, ny) )
locsGrid <- expand.grid( x=locs$x, y=locs$y )

nugget = 0.2 # tau^2
psill = 0.8 # sigma^2
nu = 0.5  # Bessel param
phi = 0.75 # range parameter

# ------------------
require(gstat)
modG <- gstat::gstat(formula = z~1, locations = ~x+y, data=locsGrid, dummy=TRUE, beta = 0,
    model = gstat::vgm(psill=psill, nugget=nugget, model="Mat", range=1/phi, kappa=nu), nmax = 20)
rfG <- predict(modG, newdata=locsGrid, nsim = 1)
sp::gridded(rfG) = ~x+y
sp::spplot(rfG)

# re-estimate params
ss = rfG[ si, ]
locs_xy = ss[, c("x","y") ]
z = ss[,"sim1"]
vg = lstfilter_variogram( locs_xy, z, methods="gstat" )


# -------------------
require(geoR)
rfGeoR = geoR::grf( nsim=1, grid="reg", nx=nx, ny=ny, xlims=1:nx, ylims=1:ny,
  cov.model="matern", cov.pars=c(sigmasq=psill, phi=phi), kappa=nu )
plot(rfGeoR)
plot(variog(rfGeoR, max.dist=100))
lines.variomodel(rfGeoR)
image(rfGeoR)
z = rfG[si,3]

vg = lstfilter_variogram( locs_xy, z, methods="geoR" )


# -------------------
require(RandomFields)

