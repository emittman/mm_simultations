library(cudarpackage)

#simulate data
sim_data <- function(G, X, group_n, mu_beta, sd_beta, corr=NULL, mu_dispersion, sd_dispersion){
  V <- length(group_n)
  N <- sum(group_n)
  group_id <- unlist(sapply(1:V, function(g) rep(g, group_n[g])))
  X <- X[group_id,]
  if(!is.null(corr)){
    ch <- chol(corpcor::rebuild.cov(corr, sd_beta^2))
    beta <- sapply(1:G, function(g) t(rnorm(V)) %*% ch + mu_beta)
  } else{
    beta <- sapply(1:G, function(g) rnorm(V, mu_beta, sd_beta))
  }
  dispersion <- rlnorm(G, mu_dispersion, sd_dispersion)
  y <- t(sapply(1:G, function(g) rnbinom(N, mu = exp(X %*% beta[,g]), size = 1/dispersion[g])))
  out <- list(formatData(y, X), truth = list(beta = beta))
  out
}  

G <- 10000
X <- diag(3)

m <- matrix(c(3, 1, 1,
              1, 3, 1,
              1, 1, 3), 3, 3)
m

prior_mean <- c(0,0,0)

cr <- cov2cor(m)
prior_sd <- sqrt(diag(m))
d <- sim_data(G, X, c(3,2,2), prior_mean, prior_sd, cr, -2, .5)
d

K <- 2500

estimates <- indEstimates(d[[1]])
priors <- formatPriors(K, prior_mean, prior_sd)
#chain <- initChain(priors, G, estimates)
system.time(
s <- mcmc(d[[1]], priors, methodPi = "stickBreaking", n_iter=500000, idx_save=c(0,499,999), thin=1,
          n_save_P=100, alpha_fixed=F, verbose = 0, warmup=10000, slice_width=1, max_steps=100, estimates=estimates)
)
saveRDS(s, "samples_stick.rds")
saveRDS(d[[2]], "truth.rds")

system.time(
s2 <- mcmc(d[[1]], priors, methodPi = "symmDirichlet", n_iter=500000, idx_save=c(0,499,999), thin=1,
          n_save_P=100, alpha_fixed=F, verbose = 0, warmup=10000, slice_width=1, max_steps=100, estimates=estimates)
)

saveRDS(s2, "samples_SD.rds")
