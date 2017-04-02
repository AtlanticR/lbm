
lbm_predictionarea = function(p, sloc) {

  pa_w_n = length(p$pa_w)

  # determine prediction locations and time slices
  iwplon = round( (sloc[1]-p$origin[1])/p$pres + 1 + p$pa_w )
  iwplat = round( (sloc[2]-p$origin[2])/p$pres + 1 + p$pa_w )

  pa = NULL
  pa = data.frame( iplon = rep.int(iwplon, pa_w_n) , 
                   iplat = rep.int(iwplat, rep.int(pa_w_n, pa_w_n)) )

  bad = which( (pa$iplon < 1 & pa$iplon > p$nplons) | (pa$iplat < 1 & pa$iplat > p$nplats) )
  if (length(bad) > 0 ) pa = pa[-bad,]
  if (nrow(pa)< 5) return(NULL)
 
  pa$i = match( array_map( "2->1", pa[, c("iplon", "iplat")], gridparams=p$gridparams ), ploc_ids )
      
  bad = which( !is.finite(pa$i))
  if (length(bad) > 0 ) pa = pa[-bad,]
  pa_n = nrow(pa)
  if ( pa_n < 5) return(NULL)

  pa$plon = Ploc[ pa$i, 1]
  pa$plat = Ploc[ pa$i, 2]

  # prediction covariates i.e., independent variables/ covariates
  pvars = c("plon", "plat", "i")
  if (p$nloccov > 0) {
    # .. not necessary except when covars are modelled locally
    for (ci in 1:p$nloccov) {
      vn = p$variables$local_cov[ci]
      pu = NULL
      pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
      nts = ncol(pu)
      if ( nts== 1 ) {
        pvars = c( pvars, vn )
        pa[,vn] = pu[pa$i]  # ie. a static variable
      }
    }
  }
  pa = pa[, pvars]

  if ( exists("TIME", p$variables) ) {
    pa = cbind( pa[ rep.int(1:pa_n, p$nt), ], 
                    rep.int(p$prediction.ts, rep(pa_n, p$nt )) )
    names(pa) = c( pvars, p$variables$TIME )

    pa = cbind( pa, lbm_timecovars ( vars=p$variables$local_all, ti=pa[,p$variables$TIME]  ) )

    if (p$nloccov > 0) {
      # add time-varying covars .. not necessary except when covars are modelled locally
      for (ci in 1:p$nloccov) {
        vn = p$variables$local_cov[ci]
        pu = NULL
        pu = lbm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
        nts = ncol(pu)
        if ( nts == p$ny )  {
          pa$iy = pa$yr - p$yrs[1] + 1 #yr index
          pa[,vn] = pu[ cbind(pa$i, pa$iy) ]  
          message("Need to check that data order is correct")
        } else if ( nts == p$nt ) {
          pa$it = p$nw*(pa$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
          pa[,vn] = pu[ cbind(pa$i, pa$it) ]  
          message("Need to check that data order is correct")
        } else if (nts==1) { } #nothing to do .. already processed above }
      }
    }
  }

  return(pa)
}



