### Gibbs sampler for (partially conjugate) analysis of Flint data

n <- 271
ybar <- 1.40
s.2 <- 1.684


mu0 <- 1.10
tau.2.0 <- sigma.2.0 <- 1.17
nu0 <- 1


mun <- function(sigma.2){
  (mu0/tau.2.0 + n*ybar/sigma.2) / (1/tau.2.0 + n/sigma.2)
}

tau.2.n <- function(sigma.2){
  1 / (1/tau.2.0 + n/sigma.2)
}

sigma.2.n <- function(mu){
  (nu0*sigma.2.0 + (n-1)*s.2 + n*(mu-ybar)^2) / (nu0 + n)
}


n.sim <- 10000
mu.sim <- numeric(n.sim)
sigma.2.sim <- numeric(n.sim)

mu.sim[1] <- 1.4         # starting value
sigma.2.sim[1] <- 1.7    # starting value

for(t in 2:n.sim){
  mu.sim[t] <- rnorm(1, mun(sigma.2.sim[t-1]),
                     sqrt(tau.2.n(sigma.2.sim[t-1])))
  
  sigma.2.sim[t] <- 1 / rgamma(1, (nu0+n)/2,
                               (nu0+n)*sigma.2.n(mu.sim[t])/2)
}


