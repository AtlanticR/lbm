

if (0) {
  devtools::install_github("kaskr/adcomp/TMB")
  install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/testing")
  install.packages( "RandomFields" )
  install.packages( "RANN" )
}


# choose data:
if (using.meuse) {
  library(sp)
  data(meuse)
  data(meuse.grid)
  Data = meuse
  Data.pred = meuse.grid
  ## var names : y=zinc, x = dist
 nugget = 0.2 # tau^2
  psill = 0.8 # sigma^2
  nu = 0.5  # Bessel param
  phi = 100 # range parameter
  # range = geoR::practicalRange("matern", phi=100, kappa=0.5)  # 300 km

  # range = distance_matern( phi=100, nu=0.5, cor=0.95)
}

if (using.random) {
  set.seed(1)
  # sim dimensions
  nx = 200
  ny = 100
  nt = 10
  nsample = 200
  si = sample( 1:(nx*ny), nsample )
  Data <- list( x= seq(1, nx), y=seq(1, ny) )
  Data.pred <- expand.grid( x=Data$x, y=Data$y )
  nugget = 0.2 # tau^2
  psill = 0.8 # sigma^2
  nu = 0.5  # Bessel param
  phi = 100 # range parameter
  # range = geoR::practicalRange("matern", phi=100, kappa=0.5)  # 300 km

  # range = distance_matern( phi=100, nu=0.5, cor=0.95)
}

  #----------------
  # using gstat
  require(gstat)
  modG <- gstat::gstat(formula = z~1, locations = ~x+y, data=Data.pred, dummy=TRUE, beta = 0,
      model = gstat::vgm(psill=psill, nugget=nugget, model="Mat", range=1/phi, kappa=nu), nmax = 20)
  rfG <- predict(modG, newdata=Data.pred, nsim = 1)
  names( rfG ) = c("x", "y", "z" )

  # -------------------
  # using geoR
  require(geoR)
  rfGeoR = geoR::grf( 100, cov.model="matern", cov.pars=c(sigmasq=psill, phi=phi), kappa=nu )

  plot(rfGeoR)
  plot(variog(rfGeoR, max.dist=100))
  lines.variomodel(rfGeoR)
  image(rfGeoR)

  # re-estimate params
  ss = data.frame( rfG[ si, ] )
  xy = ss[,c("x","y")]
  z =  ss[,"z"]

  vgs = lbm_variogram( xy, z, methods="gstat" )
  vgr = lbm_variogram( xy, z, methods="geoR" )
  vsp = lbm_variogram( xy, z, methods="spBayes" )
  vrf = lbm_variogram( xy, z, methods="RandomFields" )
  vin = lbm_variogram( xy, z, methods="inla" )

}


# ---------------------------------------------------
# basic Universal Kriging with gstat

require(gstat)

coordinates(Data) = ~x+y
gridded(Data.pred) = ~x+y

# local universal kriging
g = gstat(id = "log_zinc", formula = log(zinc)~sqrt(dist), data=Data )
vgrm = fit.variogram(variogram(g), vgm(1, "Exp", 300, 1)) # variogram of residuals

# prediction from local neighbourhoods within radius of 170 m or at least 10 points
g = gstat(id="log_zinc", formula=log(zinc)~sqrt(dist), data=Data, maxdist=170, nmin=10, force=TRUE, model=vgrm)
gpred <- predict(g, Data.pred )
spplot(gpred)




ex.bayes <- krige.bayes(rfGeoR, loc=ex.grid,
  model = model.control(cov.m="matern", kappa=nu ),
  prior = prior.control(phi.discrete=seq(0, 0.7, l=51), phi.prior="reciprocal")
)
     #
     # Prior and posterior for the parameter phi
     plot(ex.bayes, type="h", tausq.rel = FALSE, col=c("red", "blue"))
     #
     # Plot histograms with samples from the posterior
     par(mfrow=c(3,1))
     hist(ex.bayes)
     par(mfrow=c(1,1))

     # Plotting empirical variograms and some Bayesian estimates:
     # Empirical variogram
     plot(variog(rfGeoR, max.dist = 1), ylim=c(0, 15))
     # Since rfGeoR is a simulated data we can plot the line with the "true" model
     lines.variomodel(rfGeoR, lwd=2)
     # adding lines with summaries of the posterior of the binned variogram
     lines(ex.bayes, summ = mean, lwd=1, lty=2)
     lines(ex.bayes, summ = median, lwd=2, lty=2)
     # adding line with summary of the posterior of the parameters
     lines(ex.bayes, summary = "mode", post = "parameters")
     # Plotting again the empirical variogram
     plot(variog(rfGeoR, max.dist=1), ylim=c(0, 15))
     # and adding lines with median and quantiles estimates
     my.summary <- function(x){quantile(x, prob = c(0.05, 0.5, 0.95))}
     lines(ex.bayes, summ = my.summary, ty="l", lty=c(2,1,2), col=1)

     # Plotting some prediction results
     op <- par(no.readonly = TRUE)
     par(mfrow=c(2,2), mar=c(4,4,2.5,0.5), mgp = c(2,1,0))
     image(ex.bayes, main="predicted values")
     image(ex.bayes, val="variance", main="prediction variance")
     image(ex.bayes, val= "simulation", number.col=1,
           main="a simulation from the \npredictive distribution")
     image(ex.bayes, val= "simulation", number.col=2,
           main="another simulation from \nthe predictive distribution")
     #
     par(op)
     ## End(Not run)



# -------------------
require(RandomFields)

