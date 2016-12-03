
  # ----------------------------------------------------------------------------
  # randomly choose lon/lat combinations in a 10 x 10 minute grid
  #
  #   n is the number of longitude/latitude points to return
  #   gridsize is the size of the refernce grid box (ie., 10 minute squares)
  #   res is the resolution of the output ... default in seconds
  #   seed is the random number seed if manual control is desired
  #
  #   e.g. usage : coordinates = select.position.random.lonlat.grid(n=10)
  # ----------------------------------------------------------------------------

  select.position.random.lonlat.grid = function ( n=1, gridsize=10, seed=1 ) {
    if (seed != 1) .Random.seed = seed  # do we need to reset the random number generator?
    nx = gridsize * 60   # lowest resolution is in seconds
    labels = seq(1:nx)      # these contain the actual list of lons or lats

    # randomly fill a vector of length nx with random uniform numbers
      lonlat = NULL
      for (i in 1:n ) {
        lon = which.max( runif(nx) )
        lat = which.max( runif(nx) )
        lonlat = rbind(lonlat, cbind( labels[lon], labels[lat] ))
      }

      decimal = round(lonlat / 60, 3)
      minutes = trunc(decimal)
      seconds = round((decimal - minutes)*60,0)
      out = data.frame ( minutes=minutes[,1],seconds=seconds[,1],
                         minutes=minutes[,2], seconds=seconds[,2],
                         decimal=decimal )
    return (out)
  }


