### Metropolis sampler for analysis of Flint data

n <- 271
ybar <- 1.40
s.2 <- 1.684


mu0 <- 1.10
tau.2.0 <- sigma.2.0 <- 1.17
nu0 <- 1


dinvchisq <- function(x, df, scalesq){  # inverse chi-square density
  if(x>0)
    ((df/2)^(df/2) / gamma(df/2)) * sqrt(scalesq)^df * x^(-(df/2+1)) *
      exp(-df*scalesq/(2*x))
  else 0
}

likelihood <- function(mu, sigma.2){
  sigma.2^(-n/2) * exp(-(n-1)*s.2/(2*sigma.2)) *
    exp(-n*(mu-ybar)^2/(2*sigma.2))
}

ratio <- function(mu.prop, sigma.2.prop, mu.old, sigma.2.old){
  dnorm(mu.prop,mu0,sqrt(tau.2.0)) * dinvchisq(sigma.2.prop,nu0,sigma.2.0) *
    likelihood(mu.prop,sigma.2.prop)  /
  (dnorm(mu.old,mu0,sqrt(tau.2.0)) * dinvchisq(sigma.2.old,nu0,sigma.2.0) *
    likelihood(mu.old,sigma.2.old))
}


n.sim <- 10000
mu.sim <- numeric(n.sim)
sigma.2.sim <- numeric(n.sim)
accept.prob <- numeric(n.sim-1)

# rho <- 0.01

mu.sim[1] <- 1.4         # starting value
sigma.2.sim[1] <- 1.7    # starting value

for(t in 2:n.sim){
  mu.prop <- rnorm(1, mu.sim[t-1], sqrt(rho))
  sigma.2.prop <- rnorm(1, sigma.2.sim[t-1], sqrt(rho))

  accept.prob[t-1] <-
    min(ratio(mu.prop,sigma.2.prop,mu.sim[t-1],sigma.2.sim[t-1]), 1)
  if(runif(1) < accept.prob[t-1]){
    mu.sim[t] <- mu.prop
    sigma.2.sim[t] <- sigma.2.prop
  }else{
    mu.sim[t] <- mu.sim[t-1]
    sigma.2.sim[t] <- sigma.2.sim[t-1]
  }
}


# mean(accept.prob)


