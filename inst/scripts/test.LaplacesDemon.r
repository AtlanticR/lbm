

# ---------------------------------------------------
# basic Universal Kriging
library(sp)
library(gstat)

data(meuse)
coordinates(meuse) = ~x+y
# local universal kriging
gmeuse = gstat(id = "log_zinc", formula = log(zinc)~dist, data = meuse)
vmeuse.res = fit.variogram(variogram(gmeuse), vgm(1, "Exp", 300, 1)) # variogram of residuals
# prediction from local neighbourhoods within radius of 170 m or at least 10 points
gmeuse = gstat(id = "log_zinc", formula = log(zinc)~sqrt(dist),
data = meuse, maxdist=170, nmin=10, force=TRUE, model=vmeuse.res)

data(meuse.grid)
gridded(meuse.grid) = ~x+y
predmeuse = predict(gmeuse, meuse.grid)
spplot(predmeuse)



# ----------------------------

library(spBayes)

data(meuse)
coordinates(meuse) = ~x+y

##Collect samples
y = log(meuse$zinc)
x = meuse$dist
fit <- glm( y~x-1, family="gaussian")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))


coords = coordinates(meuse)

m.1 <- spLM( y~1, coords=coords, knots=c(6,6),
             starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n, size=weights, prob=p)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(coords, label=paste("(",y,")",sep=""))

surf <- mba.surf(cbind(coords,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(coords, label=paste("(",weights,")",sep=""))

# ---------------------------------------------------
# Universal Kriging with prediction via LaplacesDemon
bioLibrary( "bio.spacetime", "conker")

require(LaplacesDemonCpp)

Data = conker_LaplacesDemon_setup( DS="spatial.test", Data ) # spatial + intercept


# maximum likelihood solution
f.ml = optim( par=Data$PGF(Data), fn=Data$Model.ML, Data=Data, control=list(maxit=5000, trace=1), method="BFGS"  )
names(f.ml$par ) = Data$parm.names

# penalized maximum likelihood .. better but still a little unstable depending on algorithm
f.pml = optim( par=Data$PGF(Data), fn=Data$Model.PML, Data=Data,  control=list(maxit=5000, trace=1), method="BFGS" , hessian=TRUE )
names(f.pml$par ) = Data$parm.names
#print(sqrt( diag( solve(f.pml$hessian) )) ) # assymptotic standard errors


f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=Data$PGF(Data), Iterations=100, Status=10, Thinning=1 )
f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=100, Status=10, Thinning=1 )

f <- LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Covar=NULL, Iterations=5000, Status=100, Thinning=25,
     Algorithm="CHARM", Specs=list(alpha.star=0.44))

f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=100, Status=1, Thinning=1, Algorithm="NUTS", Covar=f$Covar, Specs=list(A=10, delta=0.6, epsilon=NULL, Lmax=Inf))  # A=burnin, delta=target acceptance rate


f = LaplaceApproximation(Data$Model, Data=Data, parm=as.initial.values(f), Iterations=1000, Method="SPG", CPUs=8 ) # fast inital solution
f = LaplacesDemon.hpc(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=10, Thinning=10 )




f = LaplaceApproximation(Data$Model, Data=Data, parm=as.initial.values(f), Method="Roptim", method="BFGS", Stop.Tolerance=1e-9, Iterations = 50  )


f = LaplaceApproximation(Data$Model, Data=Data, parm=as.initial.values(f), Iterations=10, Method="LBFGS" ) # refine it
f = VariationalBayes(Data$Model, Data=Data, parm=as.initial.values(f), Iterations=10, Samples=20, CPUs=5, Covar=f$Covar ) # refine it again


f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=10 )
Initial.Values <- as.initial.values(f)
f <- LaplacesDemon(Data$Model, Data=Data, Initial.Values,
     Covar=f$Covar, Iterations=2000, Status=142, Thinning=20,
     Algorithm="AMWG", Specs=list(Periodicity=21))


f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=100, Status=1, Thinning=1, Covar=f$Covar )
f = VariationalBayes(Data$Model, Data=Data, parm=as.initial.values(f), Iterations=100, Samples=10, CPUs=5, Covar=f$Covar )
f = IterativeQuadrature(Data$Model, Data=Data, parm=as.initial.values(f), Iterations=10, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL), Covar=f$Covar )


inew = grep( "yP", rownames( f$Summary2 ) )
m = f$Summary2[inew,]
preds = coords.new 
preds$prediction = f$Summary2[inew, 1]




---
# space-time separable
data(demontexas)
Y = as.matrix(demontexas[1:20,c(18:30)])
X = cbind(1,as.matrix(demontexas[1:20,c(1,4)])) #Static predictors
plat = demontexas[1:20,2]
plon = demontexas[1:20,3]
S = nrow(Y) #Number of sites, or points in space
T = ncol(Y) #Number of time-periods
K = ncol(X) #Number of columns in design matrix X including the intercept
D.S = as.matrix(dist(cbind(plon,plat), diag=TRUE, upper=TRUE))
D.T = as.matrix(dist(cbind(c(1:T),c(1:T)), diag=TRUE, upper=TRUE))
mon.names = "LP"
parm.names = as.parm.names(list(muKs=rep(0,S), theta=rep(0,T),
beta=rep(0,K), phi=rep(0,2), sigma=rep(0,3)))
pos.muKs = grep("muKs", parm.names)
pos.theta = grep("theta", parm.names)
pos.beta = grep("beta", parm.names)
pos.phi = grep("phi", parm.names)
pos.sigma = grep("sigma", parm.names)

PGF = function(Data) {
  beta = rnorm(Data$nKs) + c(mean(Data$Y),rep(0,Data$nKs-1))
  phi = runif(2,1,5)
  sigma = runif(3)
  kappa = 1
  rho.t = 1
  covKs = sigma[2]^2 * exp(-phi[1] * Data$D.S)^kappa
  covKt = sigma[3]^2 * exp(-phi[2] * Data$D.T)^rho.t
  muKs = as.vector(mvtnorm::rmvnorm(1, rep(0,Data$S), covKs, method="chol" ))
  theta = as.vector(mvtnorm::rmvnorm(1, rep(0,Data$T), covKt, method="chol" ))
  return(c(muKs, theta, beta, phi, sigma))
}

Data = list(D.S=D.S, D.T=D.T, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y,
  mon.names=mon.names,
  parm.names=parm.names, pos.muKs=pos.muKs, pos.theta=pos.theta,
  pos.beta=pos.beta, pos.phi=pos.phi, pos.sigma=pos.sigma)

Model = function(parm, Data) {
  ### Hyperparameters
  muKs.mu = rep(0,Data$S)
  theta.mu = rep(0,Data$T)
  ### Parameters
  beta = parm[Data$pos.beta]
  muKs = parm[Data$pos.muKs]
  theta = parm[Data$pos.theta]
  kappa = 1; rho.t = 1

  parm[Data$pos.sigma] = sigma = interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.phi] = phi = interval(parm[Data$pos.phi], 1, 5)
  covKs = sigma[2]^2 * exp(-phi[1] * Data$D.S)^kappa
  covKt = sigma[3]^2 * exp(-phi[2] * Data$D.T)^rho.t
  ### Log-Prior
  beta.prior = sum(dnormv(beta, 0, 1000, log=TRUE))
  muKs.prior = mvtnorm::dmvnorm(muKs, muKs.mu, covKs, log=TRUE)
  theta.prior = mvtnorm::dmvnorm(theta, theta.mu, covKt, log=TRUE)
  sigma.prior = sum(dhalfcauchy(sigma, 25, log=TRUE))
  phi.prior = sum(dunif(phi, 1, 5, log=TRUE))
  ### Log-Likelihood
  Theta = matrix(theta, Data$S, Data$T, byrow=TRUE)
  mu = as.vector(tcrossprod(Data$X, t(beta))) + muKs + Theta
  LL = sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP = LL + beta.prior + muKs.prior + theta.prior + sigma.prior +
  phi.prior
  Modelout = list(LP=LP, Dev=-2*LL, Monitor=LP,
  yhat=rnorm(prod(dim(mu))) * sigma[1] + mu, parm=parm)
  return(Modelout)
}

Initial.Values = c(rep(0,S), rep(0,T), rep(0,2), rep(1,2), rep(1,3))


f = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )
Initial.Values = as.initial.values(f)
f = LaplacesDemon(Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=1 )
f = VariationalBayes(Model, Data=Data, parm=as.initial.values(f), Iterations=10000, Samples=1000, CPUs=5 )
f = IterativeQuadrature(Model, Data=Data, parm=as.initial.values(f), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )



---

# space time non-separable

data(demontexas)
Y = as.matrix(demontexas[1:10,c(18:30)])
X = cbind(1,as.matrix(demontexas[1:10,c(1,4)])) #Static predictors
plat = demontexas[1:10,2]
plon = demontexas[1:10,3]
S = nrow(Y) #Number of sites, or points in space
T = ncol(Y) #Number of time-periods
K = ncol(X) #Number of columns in design matrix X including the intercept
D.S = as.matrix(dist(cbind(rep(plon,T),rep(plat,T)), diag=TRUE,
upper=TRUE))
D.T = as.matrix(dist(cbind(rep(1:T,each=S),rep(1:T,each=S)), diag=TRUE,
upper=TRUE))
mon.names = "LP"
parm.names = as.parm.names(list(Xi=matrix(0,S,T), beta=rep(0,K),
phi=rep(0,2), sigma=rep(0,2), psi=0))
pos.Xi = grep("Xi", parm.names)
pos.beta = grep("beta", parm.names)
pos.phi = grep("phi", parm.names)
pos.sigma = grep("sigma", parm.names)
pos.psi = grep("psi", parm.names)
PGF = function(Data) {
beta = rnorm(Data$nKs) + c(mean(Data$Y),rep(0,Data$nKs-1) ) 
phi = runif(2,1,5)
sigma = runif(2)
psi = runif(1)
kappa = 1
rho.t = 1
covKs = sigma[2]*sigma[2] * exp(-(Data$D.S / phi[1])^kappa - (Data$D.T
/ phi[2])^rho.t - psi*(Data$D.S / phi[1])^kappa * (Data$D.T / phi[2])^rho.t)
Xi = as.vector(rmvn(1, rep(0,Data$S*Data$T), covKs))
return(c(Xi, beta, phi, sigma, psi))
}
Data = list(D.S=D.S, D.T=D.T, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y,
mon.names=mon.names,
parm.names=parm.names, pos.Xi=pos.Xi, pos.beta=pos.beta,
pos.phi=pos.phi, pos.sigma=pos.sigma, pos.psi=pos.psi)

Model = function(parm, Data) {
  ### Hyperparameters
  Xi.mu = rep(0,Data$S*Data$T)
  ### Parameters
  beta = parm[Data$pos.beta]
  Xi = parm[Data$pos.Xi]
  kappa = 1; rho.t = 1
  sigma = interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] = sigma
  parm[Data$pos.phi] = phi = interval(parm[Data$pos.phi], 1, 5)
  parm[Data$pos.psi] = psi = interval(parm[Data$pos.psi], 1e-100, Inf)
  covKs = sigma[2]*sigma[2] * exp(-(Data$D.S / phi[1])^kappa -
  (Data$D.T / phi[2])^rho.t -
  psi*(Data$D.S / phi[1])^kappa * (Data$D.T / phi[2])^rho.t)
  ### Log-Prior
  beta.prior = sum(dnormv(beta, 0, 1000, log=TRUE))
  Xi.prior = mvtnorm::dmvnorm(Xi, Xi.mu, covKs, log=TRUE)
  sigma.prior = sum(dhalfcauchy(sigma, 25, log=TRUE))
  phi.prior = sum(dunif(phi, 1, 5, log=TRUE))
  psi.prior = dhalfcauchy(psi, 25, log=TRUE)
  ### Log-Likelihood
  Xi = matrix(Xi, Data$S, Data$T)
  mu = as.vector(tcrossprod(Data$X, t(beta))) + Xi
  LL = sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP = LL + beta.prior + Xi.prior + sigma.prior + phi.prior + psi.prior
  Modelout = list(LP=LP, Dev=-2*LL, Monitor=LP,
  yhat=rnorm(prod(dim(mu)))* sigma[1] +mu, parm=parm)
  return(Modelout)
}

Initial.Values = c(rep(0,S*T), c(mean(Y),rep(0,K-1)), rep(1,2), rep(1,2), 1)


f = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )
f = LaplacesDemon(Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=1 )
f = VariationalBayes(Model, Data=Data, parm=as.initial.values(f), Iterations=10000, Samples=1000, CPUs=5 )
f = IterativeQuadrature(Model, Data=Data, parm=as.initial.values(f), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )




----


# spacetime dynamic
data(demontexas)
Y = as.matrix(demontexas[1:20,c(18:30)])
X = cbind(1,as.matrix(demontexas[1:20,c(1,4)])) #Static predictors
plat = demontexas[1:20,2]
plon = demontexas[1:20,3]
D = as.matrix(dist(cbind(plon,plat), diag=TRUE, upper=TRUE))
S = nrow(Y) #Number of sites, or points in space
T = ncol(Y) #Number of time-periods
K = ncol(X) #Number of columns in design matrix X including the intercept
mon.names = "LP"
parm.names = as.parm.names(list(muKs=rep(0,S), beta=matrix(0,K,T),
phi=rep(0,T), kappa=rep(0,T), rho.t=rep(0,T), sigma=rep(0,4),
tau=rep(0,K)))
pos.muKs = grep("muKs", parm.names)
pos.beta = grep("beta", parm.names)
pos.phi = grep("phi", parm.names)
pos.kappa = grep("kappa", parm.names)
pos.rho.t = grep("rho.t", parm.names)
pos.sigma = grep("sigma", parm.names)

pos.tau = grep("tau", parm.names)
PGF = function(Data) {
  beta = rnorm(Data$nKs*Data$T) + rbind(mean(Data$Y),
  matrix(0, Data$nKs-1, Data$T)) 
  phi = rhalfnorm(Data$T, 1)
  kappa = rhalfnorm(Data$T, 1)
  rho.t = rhalfnorm(Data$T, 1)
  covKs = rho.t[1]*rho.t[1]*exp(-phi[1]*Data$D)^kappa[1]
  muKs = as.vector(mvtnorm::rmvnorm(1, rep(0,Data$S), covKs, method="chol" ))
  sigma = runif(4)
  tau = runif(Data$nKs)
  return(c(muKs, beta, phi, kappa, rho.t, sigma, tau))
}

Data = list(D=D, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y, mon.names=mon.names, parm.names=parm.names,
  pos.muKs=pos.muKs, pos.beta=pos.beta, pos.phi=pos.phi,
  pos.kappa=pos.kappa, pos.rho.t=pos.rho.t, pos.sigma=pos.sigma,
  pos.tau=pos.tau)

Model = function(parm, Data) {
  ### Parameters
  beta = matrix(parm[Data$pos.beta], Data$nKs, Data$T)
  muKs = parm[Data$pos.muKs]
  parm[Data$pos.phi] = phi = interval(parm[Data$pos.phi], 1e-100, Inf)
  kappa = interval(parm[Data$pos.kappa], 1e-100, Inf)
  parm[Data$pos.kappa] = kappa
  rho.t = interval(parm[Data$pos.rho.t], 1e-100, Inf)
  parm[Data$pos.rho.t] = rho.t
  sigma = interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] = sigma
  parm[Data$pos.tau] = tau = interval(parm[Data$pos.tau], 1e-100, Inf)
  covKs = array(0, dim=c(Data$S, Data$S, Data$T))
  for (t in 1:Data$T) {
  covKs[ , ,t] = rho.t[t]^2 * exp(-phi[t] * Data$D)^kappa[t]}
  ### Log-Prior
  beta.prior = sum(dnormv(beta[,1], 0, 1000, log=TRUE),
  dnorm(beta[,-1], beta[,-Data$T], matrix(tau, Data$nKs,
  Data$T-1), log=TRUE))
  muKs.prior = mvtnorm::dmvnorm(muKs, rep(0,Data$S), covKs[ , , 1], log=TRUE)
  phi.prior = sum(dhalfnorm(phi[1], sqrt(1000), log=TRUE),
  dtrunc(phi[-1], "norm", a=0, b=Inf, mean=phi[-Data$T],
  sd=sigma[2], log=TRUE))
  kappa.prior = sum(dhalfnorm(kappa[1], sqrt(1000), log=TRUE),
  dtrunc(kappa[-1], "norm", a=0, b=Inf, mean=kappa[-Data$T],
  sd=sigma[3], log=TRUE))
  rho.t.prior = sum(dhalfnorm(rho.t[1], sqrt(1000), log=TRUE),
  dtrunc(rho.t[-1], "norm", a=0, b=Inf, mean=rho.t[-Data$T],
  sd=sigma[4], log=TRUE))
  sigma.prior = sum(dhalfcauchy(sigma, 25, log=TRUE))
  tau.prior = sum(dhalfcauchy(tau, 25, log=TRUE))
  ### Log-Likelihood
  mu = tcrossprod(Data$X, t(beta))
  Theta = matrix(muKs, Data$S, Data$T)
  for (t in 2:Data$T) {
  for (s in 1:Data$S) {
  Theta[s,t] = sum(covKs[,s,t] / sum(covKs[,s,t]) * Theta[,t-1])}}
  mu = mu + Theta
  LL = sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP = LL + beta.prior + muKs.prior + sum(phi.prior) +
  sum(kappa.prior) + sum(rho.t.prior) + sigma.prior + tau.prior
  Modelout = list(LP=LP, Dev=-2*LL, Monitor=LP,
  yhat=rnorm(prod(dim(mu)))* sigma[1] + mu, parm=parm)
  return(Modelout)
}

Initial.Values = c(rep(0,S), rep(c(mean(Y),rep(0,K-1)),T), rep(1,T),
rep(1,T), rep(1,T), rep(1,4), rep(1,K))


f = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )

f = LaplacesDemon(Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=1 )
f = VariationalBayes(Model, Data=Data, parm=as.initial.values(f), Iterations=10000, Samples=1000, CPUs=5 )
f = IterativeQuadrature(Model, Data=Data, parm=as.initial.values(f), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )

