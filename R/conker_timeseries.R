
lstfilter_timeseries  = function( x, method="spec.pgram", quant=0.95, taper=0.05, kernel= kernel("modified.daniell", c(1,1)) , freq=1, plotdata=FALSE ) {
  #\\ estimate simple time series autocorrelation (serial)
  # x = sunspot.year
  u = NULL
  if (method=="spec.pgram") {
    # with spec.pgram, default is to taper 0.1 and remove linear trend
    u = spec.pgram ( x, detrend=TRUE, plot=FALSE, na.action=na.omit, taper=taper  )
  }
 
  if (method=="spec.ar") {
    # with spec.ar -- parametric AR fit/smooth via AIC then compute FFT on modeled results .. z must be ts 
    if (! is.ts(x) ) stop( "Data must be a ts for this to work" )
    u = spec.ar ( x , plot=FALSE, na.action=na.omit )
  }

  if (method=="fft") {
    # direct with FFT .. no detrending
    z = spec.taper(scale(x, TRUE, FALSE), p=taper )
    u = list()
    nx = length(z)
    nx2 = floor(nx/2)  # make even
    I = Mod(fft(z))^2 /nx 
    I[1L] = 0
    u$spec = (4/nx)*I[1:nx2]    # "scaled periodogram"
    u$freq = (0:(nx2-1))*freq/nx
  }

  if (method=="cpgram") {
    # direct copy from stats::cpgram .. just for reference
    if (! is.ts(x) ) stop( "Must be ts .. this is just to show method" )
    z = spec.taper(scale(x, TRUE, FALSE), p=taper )
    y <- Mod(fft(z))^2/length(z)
    y[1L] <- 0
    n <- length(z)
    z <- (0:(n/2)) * frequency(x)/n  
    if (length(z)%%2 == 0) {
        n <- length(z) - 1
        y <- y[1L:n]
        z <- z[1L:n]
    } else {
      y <- y[seq_along(z)]
    }
    u=list()
    u$spec = y
    u$freq = z
  }


  if (method == "inla" ) {
    # incomplete .. interpolate as a simple AR in INLA and then FFT ...
    require(INLA)
    nn = abs( diff( dat[,"x"]  ) )
    dd = median( nn[nn>0], na.rm=TRUE )
    dat$xiid = dat$x = jitter( dat$x, amount=dd / 20 ) # add noise as inla seems unhappy with duplicates in x?
    rsq = 0
    nw = length( which(is.finite( dat$y)))
    nw0 = nw + 1

    #    preds_ydata = list()
    #       preds_ydata[[ p$variables$Y ]] = NA ## ie. to predict
    #      PREDS = inla.stack( tag="preds", data=preds_ydata, A=preds_A, effects=preds_eff, remove.unused=FALSE )
    #DATA = inla.stack(DATA, PREDS )
    #    preds_stack_index = inla.stack.index( DATA, "preds")$data  # indices of predictions in stacked data
          
    r = inla( y ~ 0 + f( x, model="ar1" ), data = dat )
    dat$predictions =  r$summary.random$x[["mean"]]
    ar.pred =  r$summary.hyperpar["Rho for x", "mean" ]
    mm = glm( predictions~y, data=dat )
   
    # xmean = RES$summary.fitted.values[ stack_index, "mean"]
   
  }

  u$powerPr = cumsum( u$spec ) / sum( u$spec )
  u$quantileFreq = u$freq[ min( which( u$powerPr >= quant ) ) ]
  u$quantilePeriod = 1 / u$quantileFreq 

  if (plotdata) { 
    plot ( u$powerPr ~ u$freq, type="l" )
    abline( v=u$quantileFreq )
    abline( h=quant )
    legend( "bottomright", legend= paste( "Period =", round( u$quantilePeriod, digits=3 ) )) 
  }

  if ( 0 & is.ts( x ) ) {
    xm <- frequency(x)/2
    mp <- length(x) - 1
    crit <- 1.358/(sqrt(mp) + 0.12 + 0.11/sqrt(mp))
    oldpty <- par(pty = "s")
    on.exit(par(oldpty))
     ci.col = "blue"
    plot(z, cumsum(y)/sum(y), type = "s", xlim = c(0, xm), ylim = c(0, 
        1), xaxs = "i", yaxs = "i", xlab = "frequency", ylab = "")
    lines(c(0, xm * (1 - crit)), c(crit, 1), col = ci.col, lty = 2)
    lines(c(xm * crit, xm), c(0, 1 - crit), col = ci.col, lty = 2)
  }
 
  return (u) 

}
