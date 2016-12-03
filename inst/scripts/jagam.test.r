
http://www.petrkeil.com/?p=2385

I will use example dataset bei from spatstat package. The data are positions of 3605 individual trees of Beilschmiedia pendula (Lauraceae) in a 1000 by 500 metre rectangular sampling region in the tropical rainforest of Barro Colorado Island. The data are stored in a point process pattern ppp object.


-- jags/jagam test  

  library(mgcv)   # fits GAMs
  library(spatstat) # the source of the example data
  library(raster) # for operations with rasters
  library(R2jags) # interface between R and JAGS

  par(mai=c(0.5,0.3,0.3,0))
  plot(bei, cex=0.1, main=NULL)

 # cropping the data so that they have exactly 500 x 1000 cells
  ext <- extent(0, 1000, 0, 500)         # spatial extent of the raster
  empty <- raster(ext, nrow=25, ncol=50) # empty raster
  
  # aggregating the point data into the raster
  xy <- data.frame(x = bei$x, y = bei$y)
  rst <- rasterize(xy, empty, fun = "count")
  
  # replacing the NA values by 0
  rst[is.na(rst)] <- 0
  
  # extracting the cell values and their coordinates to a data.frame
  coord <- xyFromCell(rst,1:ncell(rst))
  count <- extract(rst, 1:ncell(rst))
  all.data <- data.frame(coord, count=count)

  plot(rst, axes=FALSE)
  points(xy, cex=0.1)

  # the gam model with s() indicating that I fit splines
  space.only <- gam(count~s(x, y), data=all.data, family = "poisson")
  # extraction of the predictions
  preds.mgcv <- as.vector(predict(space.only, type = "response"))

  # putting the predictions into a raster
  rst.mgcv <- rst
  rst.mgcv[] <- preds.mgcv
    plot(rst.mgcv, axes=FALSE)
  points(xy, cex=0.1)



jags.ready <- jagam(count~s(x, y), 
  data=all.data, 
  family="poisson", 
  file="jagam.bug")

readLines("jagam.bug")

model.fit <- jags(data=jags.ready$jags.data, 
  model.file="jagam.bug",
  parameters.to.save=c("mu"),
  n.chains=3,
  n.iter=1000,
  n.burnin=500)

preds.mu <- as.vector(model.fit$BUGSoutput$mean$mu)

 rst.jags <- rst
  rst.jags[] <- preds.mu

 plot(rst.jags, axes=FALSE)
  points(xy, cex=0.1) # adding the positions of individual trees
