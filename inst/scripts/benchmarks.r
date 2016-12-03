
require(sp)
require(microbenchmark)

data(meuse)
coordsK = data.frame( plon=meuse$y, plat=meuse$x ) 
dKK = as.matrix(dist( coordsK, diag=TRUE, upper=TRUE)) # distance matrix between knots

 # testing speed of MVN functions:
  sigma = c(1, 1)
  kappa = 1 
  phi =1
  rhoKs = exp(-phi * dKK)^kappa   ## spatial correlation
  covKs = sigma[2]*sigma[2] * rhoKs

  nKs = nrow(dKK)
  muKs = rnorm(nKs)
  muKs2 = rnorm(nKs)
  muKK = muKs %*% t(muKs2)

  microbenchmark( 
    LaplacesDemonCpp::rnorm(10e6, muKs, muKs2), 
    stats::rnorm(10e6, muKs, muKs2), 
    times=10 )
                                        expr      min        lq      mean    median        uq      max neval cld
 LaplacesDemonCpp::rnorm(1e+07, muKs, muKs2)  86.3477  86.44021  94.46323  86.70997  88.10098 161.3076    10  a 
            stats::rnorm(1e+07, muKs, muKs2) 401.6346 401.80529 413.12531 403.29420 407.89803 492.5047    10   b

  
 microbenchmark( 
    LaplacesDemonCpp::dnorm_rcpp(muKs2,  muKs, log=TRUE), 
    stats::dnorm(muKs2,  muKs, log=TRUE), 
    times=5000 )

                                                 expr    min      lq     mean median     uq
 LaplacesDemonCpp::dnorm_rcpp(muKs2, muKs, log = TRUE) 19.931 21.7695 23.96264 22.416 23.330
                 stats::dnorm(muKs2, muKs, log = TRUE) 14.931 16.9110 18.80110 17.711 18.571
 


 microbenchmark( 
    {mvnfast::dmvn(muKK, rep(0, nKs), covKs, log=TRUE)}, 
    {mvnfast::dmvn(muKK, rep(0, nKs), covKs, log=TRUE,  ncores=8 )},
    {LaplacesDemonCpp::dmvn(muKK, rep(0, nKs), covKs, log=TRUE)},
    {mvtnorm::dmvnorm(muKK, rep(0, nKs), covKs, log=TRUE)},
    {FastGP::rcpp_log_dmvnorm( covKs, rep(0, nKs), muKK, istoep=FALSE )},
    times=100 )

                                                                       expr      min       lq
                {     mvnfast::dmvn(muKK, rep(0, nKs), covKs, log = TRUE) } 3.309083 4.380223
    {     mvnfast::dmvn(muKK, rep(0, nKs), covKs, log = TRUE, ncores = 8) } 4.594609 4.861577
       {     LaplacesDemonCpp::dmvn(muKK, rep(0, nKs), covKs, log = TRUE) } 3.659659 4.712023
             {     mvtnorm::dmvnorm(muKK, rep(0, nKs), covKs, log = TRUE) } 1.125608 1.527980
 {     FastGP::rcpp_log_dmvnorm(covKs, rep(0, nKs), muKK, istoep = FALSE) } 3.369025 4.592944
     mean   median       uq       max neval  cld
 4.251953 4.447984 4.477536  7.765641  1000  b  
 5.513903 4.907827 6.181261 17.140637  1000    d
 4.890866 4.881906 4.988855 81.152338  1000   c 
 1.720808 1.595274 1.647026 86.716230  1000 a   
 4.697876 4.983066 5.001937 11.633854  1000   c 


 nsamp = 1
 microbenchmark( 
    {mvnfast::rmvn(nsamp, rep(0, nKs),  sigma[2]*sigma[2]*exp(-phi*dKK)^kappa )},
    {LaplacesDemon::rmvn(nsamp, rep(0,nKs), sigma[2]*sigma[2]*exp(-phi*dKK)^kappa)},
    {LaplacesDemonCpp::rmvn(nsamp, rep(0,nKs), sigma[2]*sigma[2]*exp(-phi*dKK)^kappa)},
    {mvtnorm::rmvnorm(nsamp, rep(0,nKs), sigma[2]*sigma[2]*exp(-phi*dKK)^kappa, method="chol")},
    {FastGP::rcpp_rmvnorm(nsamp, sigma[2]*sigma[2]*exp(-phi*dKK)^kappa, rep(0,nKs))},
    {FastGP::rcpp_rmvnorm_stable(nsamp, sigma[2]*sigma[2]*exp(-phi*dKK)^kappa, rep(0,nKs))},
    times=1000 )


                                                                                                               expr      min       lq
                     {     mvnfast::rmvn(nsamp, rep(0, nKs), sigma[2] * sigma[2] * exp(-phi *          dKK)^kappa) }  818.457  838.933
               {     LaplacesDemon::rmvn(nsamp, rep(0, nKs), sigma[2] * sigma[2] *          exp(-phi * dKK)^kappa) } 4865.806 4887.672
            {     LaplacesDemonCpp::rmvn(nsamp, rep(0, nKs), sigma[2] * sigma[2] *          exp(-phi * dKK)^kappa) } 2236.022 2257.005
 {     mvtnorm::rmvnorm(nsamp, rep(0, nKs), sigma[2] * sigma[2] *          exp(-phi * dKK)^kappa, method = "chol") } 1437.710 1464.744
              {     FastGP::rcpp_rmvnorm(nsamp, sigma[2] * sigma[2] * exp(-phi *          dKK)^kappa, rep(0, nKs)) } 1026.571 1041.773
       {     FastGP::rcpp_rmvnorm_stable(nsamp, sigma[2] * sigma[2] *          exp(-phi * dKK)^kappa, rep(0, nKs)) } 1689.792 1709.479
      mean   median        uq       max neval    cld
  861.9695  848.285  858.7925  1694.562  1000 a     
 5074.1735 4911.876 4935.2620 81460.801  1000      f
 2311.4176 2270.596 2291.4995  3252.717  1000     e 
 1548.2172 1491.574 1519.7615  2880.406  1000   c   
 1141.9567 1048.901 1058.7975 79364.708  1000  b    
 1762.8241 1720.117 1737.7110  2941.062  1000    d  
> 

Note:: LaplacesDemonCpp::rmvn is now a copy of mvnfast::rmvn
