# OR2_simulation2

rweibull.point.process = function(end.time, k, theta, delta = end.time/1000){ # Simulates a stationary speciation point process where the hazard rate depends solely on the age of the parental species (only one speciation event is allowed per time bins of legth 'delta')
  event.times = NULL
  for(t in seq(0,end.time, by=delta)){
    p = 1 - exp(-(((t+delta)/theta)^k - (t/theta)^k))  # probability of an event occuring between time t and t+delta
    z = rbinom(1,1,p)
    if(z==1) event.times = c(event.times, t+delta/2)
  }
  return(event.times)
}

generate.tree = function(theta,k,end.time,mean.lifespan, paleoPhylo=T){ # Generates a phylogenetic tree using 'rweibull.point.process' as the age-dependent speciation process and drawing species life-spans from an exponential distribution (i.e., constant extinction rates)
  t2 = rexp(1, 1/mean.lifespan)
  t2 = ifelse(t2<end.time, t2, end.time)
  S = data.frame(sp=1, gen=1, parent=NA, t1=0, t2=t2)
  sp = 1 # species number index
  g = 1  # generation number index
  while(any(S$t1[S$gen==g] < end.time)){
    for(p in S$sp[S$gen == g]){
      speciation.times = S$t1[S$sp==p] + rweibull.point.process(S$t2[S$sp==p] - S$t1[S$sp==p], k, theta)
      t2 = speciation.times + rexp(length(speciation.times), 1/mean.lifespan)
      temp = data.frame(t1=speciation.times, t2=ifelse(t2<end.time,t2,end.time))
      temp2 = temp[temp$t1<end.time,]
      n = nrow(temp2)
      if(n>0)
        new.species = cbind(sp=(sp+1):(sp+n), gen=rep(g+1,n), parent=rep(p,n), temp2)
      else
        new.species = NULL
      S = rbind(S, new.species)
      sp = sp + n
    }
    g = g+1
  }
  if(paleoPhylo) S = with(S, as.paleoPhylo(nm=sp, pn=parent, st=-(t1-end.time), en=-(t2-end.time), label=sp))
  return(S)
}

nll.const = function(Theta,DD){
  theta = exp(Theta)
  t1 = DD$bin.start
  t2 = c(DD$bin.start[-1], DD$bin.start[nrow(DD)]+(DD$bin.start[2]-DD$bin.start[1]))
  lambda = (t2 - t1)/theta
  lambda = ifelse(lambda>1e-16,lambda,1e-16) # For numerical robustness (avoiding numerical log(0) at extreme values)
  -sum(DD$X*log(DD$N*lambda) - DD$N*lambda)
}

nll.weibull = function(Theta,DD){
  theta = exp(Theta[1]) # using a log-link to get positive values of theta (scale parameter)
  k = exp(Theta[2])  # using a log-link to get positive values of k (shape parameter)
  t1 = DD$bin.start
  t2 = c(DD$bin.start[-1], DD$bin.start[nrow(DD)]+(DD$bin.start[2]-DD$bin.start[1]))
  lambda = (t2^k - t1^k)/(theta^k)
  lambda = ifelse(lambda>1e-16,lambda,1e-16) # For numerical robustness (avoiding numerical log(0) at extreme values)
  -sum(DD$X*log(DD$N*lambda) - DD$N*lambda)
}

nll.cpw = function(Theta,DD){ # constant + Weibull
  lam = exp(Theta[1])
  theta = exp(Theta[2]) # using a log-link to get positive values of theta (scale parameter)
  k = exp(Theta[3])  # using a log-link to get positive values of k (shape parameter)
  t1 = DD$bin.start
  t2 = c(DD$bin.start[-1], DD$bin.start[nrow(DD)]+(DD$bin.start[2]-DD$bin.start[1]))
  lambda = lam*(t2-t1) + (t2^k - t1^k)/(theta^k)
  lambda = ifelse(lambda>1e-16,lambda,1e-16) # For numerical robustness (avoiding numerical log(0) at extreme values)
  -sum(DD$X*log(DD$N*lambda) - DD$N*lambda)
}

# Gompertz-Makeham
nll.GM = function(Theta,DD){
  alpha = exp(Theta[1])
  beta = exp(Theta[2])
  lam = exp(Theta[3])
  t1 = DD$bin.start
  t2 = c(DD$bin.start[-1], DD$bin.start[nrow(DD)]+(DD$bin.start[2]-DD$bin.start[1]))
  lambda = lam*(t2-t1) + (alpha/beta)*(exp(beta*t2)-exp(beta*t1))
  lambda = ifelse(lambda>1e-16,lambda,1e-16) # For numerical robustness (avoiding numerical log(0) at extreme values)
  -sum(DD$X*log(DD$N*lambda) - DD$N*lambda)
}

# Gompert-Makeham (with beta allowed to be negative)
nll.GM2 = function(Theta,DD){
  alpha = exp(Theta[1])
  beta = Theta[2]
  lam = exp(Theta[3])
  t1 = DD$bin.start
  t2 = c(DD$bin.start[-1], DD$bin.start[nrow(DD)]+(DD$bin.start[2]-DD$bin.start[1]))
  lambda = lam*(t2-t1) + (alpha/beta)*(exp(beta*t2)-exp(beta*t1))
  lambda = ifelse(lambda>1e-16,lambda,1e-16) # For numerical robustness (avoiding numerical log(0) at extreme values)
  -sum(DD$X*log(DD$N*lambda) - DD$N*lambda)
}

nll.doubleWeibull = function(Theta,DD){
  theta1 = exp(Theta[1]) # using a log-link to get positive values of theta (scale parameter)
  k1 = exp(Theta[2])  # using a log-link to get positive values of k (shape parameter)
  theta2 = exp(Theta[3]) # using a log-link to get positive values of theta (scale parameter)
  k2 = exp(Theta[4])  # using a log-link to get positive values of k (shape parameter)
  t1 = DD$bin.start
  t2 = c(DD$bin.start[-1], DD$bin.start[nrow(DD)]+(DD$bin.start[2]-DD$bin.start[1]))
  lambda = (t2^k1 - t1^k1)/(theta1^k1) + (t2^k2 - t1^k2)/(theta2^k2)
  lambda = ifelse(lambda>1e-16,lambda,1e-16) # For numerical robustness (avoiding numerical log(0) at extreme values)
  -sum(DD$X*log(DD$N*lambda) - DD$N*lambda)
}

nll.IDB = function(Theta,DD){
  delta = exp(Theta[1])
  alpha = exp(Theta[2])
  beta = exp(Theta[3])
  t1 = DD$bin.start
  t2 = c(DD$bin.start[-1], DD$bin.start[nrow(DD)]+(DD$bin.start[2]-DD$bin.start[1]))
  lambda = (1/2)*(delta*(t2^2-t1^2) + t1^2 -t2^2) + (alpha/beta)*(log(t2*beta+1) - log(t1*beta+1))
  lambda = ifelse(lambda>1e-16,lambda,1e-16) # For numerical robustness (avoiding numerical log(0) at extreme values)
  -sum(DD$X*log(DD$N*lambda) - DD$N*lambda)
}

nll.Chen = function(Theta,DD){
  lam = exp(Theta[1])
  beta = exp(Theta[2])
  t1 = DD$bin.start
  t2 = c(DD$bin.start[-1], DD$bin.start[nrow(DD)]+(DD$bin.start[2]-DD$bin.start[1]))
  lambda = lam*(exp(t2^beta)-exp(t1^beta))
  lambda = ifelse(lambda>1e-16,lambda,1e-16) # For numerical robustness (avoiding numerical log(0) at extreme values)
  -sum(DD$X*log(DD$N*lambda) - DD$N*lambda)
}

nll.Chen.scaled = function(Theta,DD){
  lam = exp(Theta[1])
  beta = exp(Theta[2])
  theta = exp(Theta[3])
  t1 = DD$bin.start
  t2 = c(DD$bin.start[-1], DD$bin.start[nrow(DD)]+(DD$bin.start[2]-DD$bin.start[1]))
  lambda = lam*theta*(exp((t2^beta)*(theta^-beta))-exp((t1^beta)*(theta^-beta)))
  lambda = ifelse(lambda>1e-16,lambda,1e-16) # For numerical robustness (avoiding numerical log(0) at extreme values)
  -sum(DD$X*log(DD$N*lambda) - DD$N*lambda)
}

fit.tree = function(S, bin.length = 0.005, nll, dev=F, start, ...){
  # Check if all ID's are unique
  if (length(S$nm) != length(unique(S$nm)))
    stop(paste("ID codes are not unique, specifically", S$nm[duplicated(S$nm)]))
  age.of.death = -(S$en - S$st)
  age.of.parent.at.birth = -(S$st - S$st[match(S$pn,S$nm)])
  bins.start = seq(0, max(age.of.parent.at.birth, na.rm=T), bin.length)
  DD = data.frame(bin.start = bins.start)
  DD$X = DD$N = NA
  i=0
  for(s in bins.start){
    i=i+1
    DD$X[i] = sum(age.of.parent.at.birth > s & age.of.parent.at.birth <= (s+bin.length), na.rm=T)
    DD$N[i] = sum(age.of.death > s, na.rm=T)
  }
  fit = optim(start, nll, DD=DD, hessian=T, ...)
  aic = 2*length(fit$par) + 2*fit$value
  if(dev){
    s = rep(sum(DD$X)/nrow(DD), nrow(DD))
    dev = 2*(optim(s, nll.const, DD=DD, hessian=T, ...)$value - fit$value)
    df = nrow(DD) - length(fit$par)
  }
  else{dev = NA; df = NA}
  vcov = solve(fit$hessian)
  se = sqrt(diag(vcov))
  estimates = data.frame(par = fit$par, se = se, estimate = exp(fit$par), lower = exp(fit$par - 2*se), upper = exp(fit$par + 2*se))
  out = list(fit, list(ages.at.death = age.of.death, ages.at.births = age.of.parent.at.birth, vcov=vcov, se=se, estimates=estimates, deviance = dev, df = df, AIC=aic))
  return(unlist(out, recursive=F))
}

# For predictions (and SE) of the hazard function (eq. 1)
f.weibull = deriv(~logk + (exp(logk)-1)*log(x) - exp(logk)*logtheta, c("logtheta","logk"), function(logtheta, logk,x){})
f.doubleWeibull = deriv(~log((exp(logk1)/exp(logtheta1))*(x/exp(logtheta1))^(exp(logk1)-1) + (exp(logk2)/exp(logtheta2))*(x/exp(logtheta2))^(exp(logk2)-1)), c("logtheta1","logk1","logtheta2","logk2"), function(logtheta1,logk1,logtheta2,logk2,x){})
f.GM = deriv(~log(exp(logalpha)*exp(exp(logbeta)*x) + exp(loglambda)), c("logalpha","logbeta","loglambda"), function(logalpha,logbeta,loglambda,x){})
f.GM2 = deriv(~log(exp(logalpha)*exp(beta*x) + exp(loglambda)), c("logalpha","beta","loglambda"), function(logalpha,beta,loglambda,x){})
f.Chen.scaled = deriv(~llambda + lbeta + (exp(lbeta)-1)*(log(x)-ltheta) + (x/ltheta)^(exp(lbeta)), c("llambda","lbeta","ltheta"), function(llambda,lbeta,ltheta,x){})
f.IDB = deriv(~log((exp(logdelta)-1)*x + exp(logalpha)/(1+exp(logbeta)*x)), c("logdelta","logalpha","logbeta"), function(logdelta,logalpha,logbeta,x){})

prediction.const = function(fit,...){  # using the delta-method
  logpred = rep(-fit$par,length(x))
  SElogpred = rep(fit$se, length(x))
  list(
    pred = exp(logpred),
    lower = exp(logpred - 2*SElogpred),   # assuming ~normal on log-scale
    upper = exp(logpred + 2*SElogpred)
  )
}

prediction.weibull = function(fit,func,...){  # using the delta-method
  logpred = f.weibull(fit$par[1], fit$par[2],...)
  A = attr(logpred,"gradient")
  SElogpred = sqrt(diag(A %*% fit$vcov %*% t(A)))
  list(
    pred = exp(logpred),
    lower = exp(logpred - 2*SElogpred),   # assuming ~normal on log-scale
    upper = exp(logpred + 2*SElogpred)
  )
}

prediction.doubleWeibull = function(fit,func,...){  # using the delta-method
  logpred = f.doubleWeibull(fit$par[1], fit$par[2], fit$par[3], fit$par[4], ...)
  A = attr(logpred,"gradient")
  SElogpred = sqrt(diag(A %*% fit$vcov %*% t(A)))
  list(
    pred = exp(logpred),
    lower = exp(logpred - 2*SElogpred),   # assuming ~normal on log-scale
    upper = exp(logpred + 2*SElogpred)
  )
}

prediction.GM = function(fit,func,...){  # using the delta-method
  logpred = f.GM(fit$par[1], fit$par[2], fit$par[3],...)
  A = attr(logpred,"gradient")
  SElogpred = sqrt(diag(A %*% fit$vcov %*% t(A)))
  list(
    pred = exp(logpred),
    lower = exp(logpred - 2*SElogpred),   # assuming ~normal on log-scale
    upper = exp(logpred + 2*SElogpred)
  )
}

prediction.GM2 = function(fit,func,...){  # using the delta-method
  logpred = f.GM2(fit$par[1], fit$par[2], fit$par[3],...)
  A = attr(logpred,"gradient")
  SElogpred = sqrt(diag(A %*% fit$vcov %*% t(A)))
  list(
    pred = exp(logpred),
    lower = exp(logpred - 2*SElogpred),   # assuming ~normal on log-scale
    upper = exp(logpred + 2*SElogpred)
  )
}

prediction.Chen.scaled = function(fit,func,...){  # using the delta-method
  logpred = f.Chen.scaled(fit$par[1], fit$par[2], fit$par[3],...)
  A = attr(logpred,"gradient")
  SElogpred = sqrt(diag(A %*% fit$vcov %*% t(A)))
  list(
    pred = exp(logpred),
    lower = exp(logpred - 2*SElogpred),   # assuming ~normal on log-scale
    upper = exp(logpred + 2*SElogpred)
  )
}

prediction.IDB = function(fit,func,...){  # using the delta-method
  logpred = f.IDB(fit$par[1], fit$par[2], fit$par[3],...)
  A = attr(logpred,"gradient")
  SElogpred = sqrt(diag(A %*% fit$vcov %*% t(A)))
  list(
    pred = exp(logpred),
    lower = exp(logpred - 2*SElogpred),   # assuming ~normal on log-scale
    upper = exp(logpred + 2*SElogpred)
  )
}


### EXAMPLE ###

library(ape) # on CRAN
library(paleoPhylo) #download from http://r-forge.r-project.org/R/?group_id=266

# Generate a random tree
S = generate.tree(theta=4, k=1.5, end.time=40, mean.lifespan=4) # NB! Random tree has random size and may not be sufficient for model fitting
 # repeat until getting sufficiently large tree for model fitting

# Draw the tree
drawPhylo(S, addTimeLine="none", l2r=TRUE, cexTime=2,cexLab=1, dumpLast=TRUE )

# Rescaling time such that maximum lifespan is about 1 to avoid numerical problems in some cases
S2 = S
scale = max(-(S2$en-S2$st))
S2$st = S2$st/scale
S2$en = S2$en/scale

# Fitting WEIBULL model on rescaled time data
st = c(0,0) # starting parameter values (log scale) for numerical optimization
names(st) = c("theta","k")
(fit.Weibull = fit.tree(S2, nll=nll.weibull, start = st))

# Back-transforming estimate of theta to original scale
fit.Weibull$estimates["theta",3:5]*scale

# Plotting fitted predictions of age-dependant speciation rates
x = seq(0.01, max(-(S$en - S$st), na.rm=T), length.out=100) # x-values to plot
pred.Weibull = prediction.weibull(fit.Weibull, x=x/scale)

Ylim = c(0,max(pred.Weibull$upper/scale))
plot(x, pred.Weibull$pred/scale, ylim=Ylim, type="l", xlab="Age (MY)", ylab = "Speziation hazard rate")
lines(x, pred.Weibull$lower/scale,lty=2)
lines(x, pred.Weibull$upper/scale,lty=2)
axis(1, at=scale*fit.Weibull$ages.at.death, lab=F, tcl=0.5, col="red")   # Ages at death shown at the bottom
axis(3, at=scale*fit.Weibull$ages.at.birth, lab=F, tcl=0.5, col="blue")  # Parental ages at birth shown at the top

