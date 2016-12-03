# the following is based upon work posted by
# James Thorsen for his course:
# FSH 507 "Spatio-temporal models for ecologists" in
# Spring quarter 2016 at University of Washington
# https://github.com/James-Thorson/2016_Spatio-temporal_models

if (0) {
  devtools::install_github("kaskr/adcomp/TMB")
  install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/testing")
  install.packages( "RandomFields" )
  install.packages( "RANN" )
}


# sim dimensions
nx = 200
ny = 100
nt = 10
nsample = 100

locs <- list( x= seq(1, nx), y=seq(1, ny) )
locsGrid <- expand.grid( x=locs$x, y=locs$y )

# unconditional simulation
require(gstat)
modG <- gstat::gstat(formula = z~1, locations = ~x+y, data=locsGrid, dummy=TRUE, beta = 0,
    model = vgm(1,"Exp",15), nmax = 20)
rfG <- predict(modG, newdata=locsGrid, nsim = 1)
sp::gridded(rfG) = ~x+y
sp::spplot(rfG)



#Simulate a Gaussian random field
library(fields)
obj = fields::matern.image.cov( grid=locs, theta=15, smoothness=0.5, setup=TRUE)
rfF = fields::sim.rf( obj)
# Now simulate another ...
# image.plot( locs$x, locs$y, rfF, axes=FALSE  )
# quilt.plot( locsGrid$x, locsGrid$y, rfF , nrow=nx, ncol=ny, add=F) # alternate plot
image(rfF)

# Simulate data
SimList =list()
SimList$loc_xy = locsGrid[ sample( 1:nrow(locsGrid), nsample ), ]

# Make triangulated mesh
mesh = INLA::inla.mesh.create( SimList$loc_xy )
spde = INLA::inla.spde2.matern(mesh)

# Area for each location
loc_extrapolation = expand.grid( "x"=seq(0,1,length=1e3), "y"=seq(0,1,length=1e3) )
NN_extrapolation = RANN::nn2( data=SimList$loc_xy, query=loc_extrapolation, k=1 )
a_s = table(factor(NN_extrapolation$nn.idx,levels=1:nrow(SimList$loc_xy))) / nrow(loc_extrapolation)

# Make inputs
Data = list("n_s"=nsample, "n_t"=nt,
            "a_s"=a_s, "c_i"=SimList$DF$c_i, "s_i"=SimList$DF$s_i-1, "t_i"=SimList$DF$t_i-1,
            "M0"=spde$param.inla$M0, "M1"=spde$param.inla$M1, "M2"=spde$param.inla$M2)
Params = list("beta0"=0, "ln_tau_O"=log(1), "ln_tau_E"=log(1), "ln_kappa"=1,
              "omega_s"=rep(0,mesh$n), "epsilon_st"=matrix(0,nrow=mesh$n,ncol=Data$n_t))
Random = c("omega_s", "epsilon_st")

# Compile
# Libraries and functions
library( TMB )
TMB::compile( fn )
dyn.load( dynlib(Version) )
Obj = MakeADFun( data=Data, parameters=Params, random=Random )
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list(trace=1, eval.max=1e4, iter.max=1e4))
Opt[["final_diagnostics"]] = data.frame( "Name"=names(Obj$par), "final_gradient"=Obj$gr(Opt$par))
Report = Obj$report()
unlist( Report[c('Range','SigmaO','SigmaE')] )
SD = sdreport( Obj, bias.correct=TRUE)


tmb.model = paste0( "

#include <TMB.hpp>
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_s );
  DATA_INTEGER( n_t );

  DATA_VECTOR( a_s );  // Area associated with location s
  DATA_VECTOR( c_i );  // counts for observation i
  DATA_FACTOR( s_i );  // Random effect index for observation i
  DATA_FACTOR( t_i );  // Random effect index for observation i

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_tau_O );
  PARAMETER( ln_tau_E );
  PARAMETER( ln_kappa );

  // Random effects
  PARAMETER_VECTOR( omega_s );
  PARAMETER_ARRAY( epsilon_st );

  // Objective funcction
  using namespace density;
  int n_i = c_i.size();
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_O) * exp(2*ln_kappa));
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau_E) * exp(2*ln_kappa));

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  jnll_comp(1) += SCALE( GMRF(Q), 1/exp(ln_tau_O) )( omega_s );
  for( int t=0; t<n_t; t++){
    jnll_comp(2) += SCALE( GMRF(Q), 1/exp(ln_tau_E) )( epsilon_st.col(t) );
  }

  // True density and abundance
  array<Type> log_d_st( n_s, n_t );
  for( int t=0; t<n_t; t++){
  for( int s=0; s<n_s; s++){
    log_d_st(s,t) = beta0 + omega_s(s) + epsilon_st(s,t);
  }}

  // Probability of data conditional on random effects
  for( int i=0; i<n_i; i++){
    if( !isNA(c_i(i)) ) jnll_comp(0) -= dpois( c_i(i), exp(log_d_st(s_i(i),t_i(i))), true );
  }

  // Reporting
  Type jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( Range );
  REPORT( SigmaE );
  REPORT( SigmaO );
  REPORT( log_d_st );

  return jnll;
}
")


fn = tempfile(pattern="tmb_", tmpdir=getwd(), fileext=".cpp" )
cat( tmb.model, file=fn)




