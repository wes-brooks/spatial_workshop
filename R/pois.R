library("spBayes")
library("fields")
library("MBA")
library('deldir')
library('MASS')
library('scales')

# simulate abundance data
n <- 50
coords <- data.frame(x=runif(n, 0, 1), y=runif(n, 0, 1))
phi <- 6
sigma.sq <- 2
R <- exp(-phi * iDist(coords))
w <- mvrnorm(1, rep(0, n), sigma.sq * R)
w <- w - mean(w)
beta.0 <- 0.1
beta.1 <- 0.1
y <- rpois(n, exp(beta.0 + beta.1*coords$x + w))
x <- coords$x - mean(coords$x)

# plot the simulated counts
dat <- coords
dat$ct <- y
symbols(x=dat$x[dat$ct!=0], y=dat$y[dat$ct!=0], circles=dat$ct[dat$ct!=0]/2000, bg="#00ff0088",
        xlab='x', ylab='y', main='circle size is proportional to abundance', bty='n')
points(x=dat$x[dat$ct==0], y=dat$y[dat$ct==0], pch=4, col='red')

# plot the underlying random field
tri <- deldir(coords)
fill <- col_numeric("Reds", domain=NULL)(w)
plot(tile.list(tri), fillcol=fill)

# begin with coefficients estimated by simple GLM
pois.nonsp <- glm(y ~ x, family = "poisson")
beta.starting <- coefficients(pois.nonsp)
beta.tuning <- t(chol(vcov(pois.nonsp)))


# set up an MCMC chain
n.batch <- 500
batch.length <- 50
n.samples <- n.batch * batch.length
pois.sp.chain.1 <- spGLM(y ~ x, family = "poisson", coords = as.matrix(coords),
                         starting = list(beta = beta.starting, phi = 6, sigma.sq = 1, w = 0),
                         tuning = list(beta = c(0.1, 0.1), phi = 0.5, sigma.sq = 0.1, w = 0.1),
                         priors = list("beta.Flat", phi.Unif = c(3, 30), sigma.sq.IG = c(2, 1)),
                         amcmc = list(n.batch = n.batch, batch.length = batch.length, accept.rate = 0.43),
                         cov.model = "exponential", verbose = TRUE, n.report = 500)


pois.sp.chain.2 <- spGLM(y ~ 1, family = "poisson", coords = as.matrix(coords),
                         starting = list(beta = beta.starting[1], phi = 6, sigma.sq = 1, w = 0),
                         tuning = list(beta = 0.1, phi = 0.5, sigma.sq = 0.1, w = 0.1),
                         priors = list("beta.Flat", phi.Unif = c(3, 30), sigma.sq.IG = c(2, 1)),
                         amcmc = list(n.batch = n.batch, batch.length = batch.length, accept.rate = 0.43),
                         cov.model = "exponential", verbose = TRUE, n.report = 500)



# observe the Markov chains
# samps <- mcmc.list(pois.sp.chain.1$p.beta.theta.samples)#, pois.sp.chain.2$p.beta.theta.samples)#, pois.sp.chain.3$p.beta.theta.samples)
plot(mcmc.list(pois.sp.chain.1$p.beta.theta.samples))
plot(mcmc.list(pois.sp.chain.2$p.beta.theta.samples))

# observe some diagnostics
print(gelman.diag(samps))
gelman.plot(samps)
burn.in <- 15000
print(round(summary(window(samps, start = burn.in))$quantiles[, c(3, 1, 5)], 2))
