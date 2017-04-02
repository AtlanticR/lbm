
lbm_timecovars = function( vars, ti) {
    # where time exists and there are seasonal components, 
    # additional variables are created/needed here: cos.w, sin.w, etc.. 
    # for harmonic analysis: to add an offset to a trig function (b) must add cos to a sin function
    # y ~ a + c*sin(x+b)
    # y ~ a + c*sin(b)*cos(x) + c*cos(b)*sin(x)  
    #   .. as C*sin(x+b) = C*( cos(b) * sin(x) + sin(b) * cos(x) )
    # y ~ b0 + b1*x1 + b2*x2
    # where: 
    #   a = b0
    #   c^2 = b1^2 + b2^2 = c^2*(sin^2(b) + cos^2(b))
    #   c = sqrt(b1^2 + b2^2)
    #   b1/b2 = tan(b)  
    #   b = arctan(b1/b2)
    # more than 3 harmonics would not be advisable .. but you would add them here..
  out = data.frame( ti=ti )
  if ("yr" %in% vars)     out$yr = trunc( ti )
  if ("dyear" %in% vars)  out$dyear = ti - out$yr  # fractional year
  if ("cos.w" %in% vars)  out$cos.w  = cos( ti )
  if ("sin.w" %in% vars)  out$sin.w  = sin( ti )
  if ("cos.w2" %in% vars) out$cos.w2 = cos( 2*ti )
  if ("sin.w2" %in% vars) out$sin.w2 = sin( 2*ti )
  if ("cos.w3" %in% vars) out$cos.w3 = cos( 3*ti )
  if ("sin.w3" %in% vars) out$sin.w3 = sin( 3*ti )
  out$ti = NULL
  return(out)
}
