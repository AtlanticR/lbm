

library("R2BayesX")

require(sp)
require( rgeos) 

data(meuse)
meuse$ID = as.character( 1:nrow(meuse) )
rownames(meuse) = meuse$ID 
head(meuse)

coordinates(meuse) = c("x", "y") 

# kriging
fm1 <- bayesx( log(zinc) ~ sx(x,y, bs="kr" ), family = "gaussian", method = "REML", data = as.data.frame(meuse) )


# MRF / GK using areal/polygonal basis:
drange = range( c( diff(range( meuse$x )), diff(range(meuse$y)) ) )
spbuffer =  floor( min(drange)/ 25 )

# define boundary of points if no boundary -- could also use convex hull ...
meuse.boundary = gBuffer( gUnaryUnion( gBuffer( meuse, width=spbuffer, byid=TRUE) ), width=spbuffer)
plot(meuse.boundary)

meuse.boundary2 = bio.utilities::concave.hull( coordinates(meuse), ub=500 )
plot(meuse.boundary2)

# triangulate and tessilate
vd = deldir::deldir( meuse$x, meuse$y, z=meuse$ID )
w = deldir::tile.list(vd)

polys = vector(mode='list', length=length(w)) 
  for (i in seq(along=polys)) {
      pcrds = cbind(w[[i]]$x, w[[i]]$y)
      pcrds = rbind(pcrds, pcrds[1,])
      polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i) )
  }

SP = SpatialPolygons(polys)

row.names(meuse.boundary) = ""
SP = gIntersection(  meuse.boundary, SP, byid=TRUE ) # crop
row.names(SP) =  gsub( "[[:space:]]", "", row.names(SP) )

SP = SP[order(as.numeric(row.names(SP)))]
meuse = meuse[ order(as.numeric(meuse$ID)), ]
meuse.tess = SpatialPolygonsDataFrame(SP, data=as.data.frame(meuse) )  

# make sure it is OK
plot( coordinates(meuse)[,2] ~ coordinates(meuse.tess)[,2] )
plot( coordinates(meuse)[,1] ~ coordinates(meuse.tess)[,1] )
plot( meuse$zinc ~ meuse.tess$zinc )


mypalette = colorRampPalette(c("darkblue","blue3", "green", "yellow", "orange","red3", "darkred"), space = "Lab")(100)
# mypalette = rev( heat.colors( 150 ) )
#mypalette <- brewer.pal(8, "YlOrRd")

meuse$logzinc = log( meuse$zinc)
meuse.tess$logzinc = log( meuse.tess$zinc )

spplot( meuse.tess, "logzinc", col.regions=mypalette )
x11()
spplot( meuse, "logzinc", col.regions=mypalette )

bxmap = sp2bnd( meuse.tess )

fm2 <- bayesx( log(zinc) ~  sqrt(dist), family = "gaussian", method = "MCMC", data = as.data.frame(meuse) )

fm2 <- bayesx( log(zinc) ~  sx(sqrt(dist)) + sx(ID, bs="gk", map=bxmap), family = "gaussian", method = "MCMC", data = as.data.frame(meuse) )


summary(fm2)

## compare effects for elevation and inclination
plot( fm2, term = "sx(x,y)", image=TRUE )

## extract fitted effects
f <- fitted( fm2, term = "sx(ID)")
     
## now use plot3d
plot3d(f)
plot( fm2, map=bxmap)


##Take a look
require(MBA)

par(mfrow=c(1,2))
coords = coordinates(meuse)  
surf <- mba.surf(cbind(coords,f$Mean),no.X=100, no.Y=100, extend=FALSE)$xyz.est
image(surf, main="Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)



## load forest health data and tree location map
data("ForestHealth")
data("BeechBnd")

## estimate model without spatial effect
fm1 <- bayesx(defoliation ~  stand + fertilized + 
  humus + moisture + alkali + ph + soil + 
  sx(age) + sx(inclination) + sx(canopy) + 
  sx(year) + sx(elevation), family = "cumlogit", 
  method = "MCMC", data = ForestHealth,  iter = 1200, burnin = 200)
summary(fm1)
plot(fm1, term = c("sx(age)", "sx(inclination)", 
  "sx(canopy)", "sx(year)", "sx(elevation)"))

## now include spatial effect
## warning: long runtime
fm2 <- bayesx(defoliation ~  stand + fertilized + 
  humus + moisture + alkali + ph + soil + 
  sx(age) + sx(inclination) + sx(canopy) + sx(year) + 
  sx(elevation) + sx(id, bs = "gk", map = BeechBnd, full = TRUE),
  family = "cumlogit", method = "REML", data = ForestHealth)
summary(fm2)
plot(fm2, term = c("sx(age)", "sx(inclination)", 
  "sx(canopy)", "sx(year)", "sx(elevation)", "sx(id)"),
  map = BeechBnd, pos = "topleft")

## compare effects for elevation and inclination
plot(c(fm1, fm2), term = c("sx(elevation)", "sx(inclination)"))



## Not run: 
## more examples
set.seed(111)
n <- 500

## regressors
dat <- data.frame(x = runif(n, -3, 3), z = runif(n, -3, 3),
  w = runif(n, 0, 6), fac = factor(rep(1:10, n/10)))

## response
dat$y <- with(dat, 1.5 + sin(x) + cos(z) * sin(w) +
  c(2.67, 5, 6, 3, 4, 2, 6, 7, 9, 7.5)[fac] + rnorm(n, sd = 0.6))

## estimate models with
## bayesx MCMC and REML
## and compare with
## mgcv gam()
b1 <- bayesx(y ~ sx(x) + sx(z, w, bs = "te") + fac,
  data = dat, method = "MCMC")
b2 <- bayesx(y ~ sx(x) + sx(z, w, bs = "te") + fac,
  data = dat, method = "REML")
b3 <- gam(y ~ s(x, bs = "ps") + te(z, w, bs = "ps") + fac, 
  data = dat)

## summary statistics
summary(b1)
summary(b2)
summary(b3)


## plot the effects
op <- par(no.readonly = TRUE)
par(mfrow = c(3, 2))
plot(b1, term = "sx(x)")
plot(b1, term = "sx(z,w)")
plot(b2, term = "sx(x)")
plot(b2, term = "sx(z,w)")
plot(b3, select = 1)
vis.gam(b3, c("z","w"), theta = 40, phi = 40)
par(op)

## combine models b1 and b2
b <- c(b1, b2)

## summary
summary(b)

## only plot effect 2 of both models
plot(b, term = "sx(z,w)") 

## with residuals
plot(b, term = "sx(z,w)", residuals = TRUE) 

## same model with kriging
b <- bayesx(y ~ sx(x) + sx(z, w, bs = "kr") , method = "MCMC", data = dat)
plot(b)

