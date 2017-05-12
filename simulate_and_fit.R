library(cudarpackage)

set.seed(134312517)

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
  out <- list(formatData(y, X), raw_counts = y, truth = beta)
  out
}  

G <- 10000
X <- matrix(c(1, -.5, 0,
              1, .5, 0,
              1, 0, 1), 3, 3, byrow=T)

m <- matrix(c(3, 1, -.5,
              1, 2, .5,
              -.5, .5, 1), 3, 3)
m

prior_mean <- c(3,0,0)

cr <- cov2cor(m)
prior_sd <- sqrt(diag(m))
d <- sim_data(G, X, c(4,4,6), prior_mean, prior_sd, cr, -2, .5)
d

K <- 3000

estimates <- indEstimates(d[[1]])
counts <- d[[2]]
truth <- d[[3]]
saveRDS(estimates, "estimates.rds")
saveRDS(counts, "counts.rds")
saveRDS(truth, "truth.rds")

prior_from_est <- informPriors(estimates)
priors <- formatPriors(K = K, 
                       prior_mean = prior_mean,
                       prior_sd = prior_sd,
                       a = prior_from_est$a,
                       b = prior_from_est$b,
                       A=5, B=.5) #vague prior on alpha
#chain <- initChain(priors, G, estimates)
system.time(
s <- mcmc(d[[1]], priors, methodPi = "stickBreaking", n_iter=100000, idx_save=c(0,499,999), thin=1,
          n_save_P=100, alpha_fixed=F, verbose = 0, warmup=30000, estimates=estimates)
)
saveRDS(s, "samples_stick.rds")

system.time(
s2 <- mcmc(d[[1]], priors, methodPi = "symmDirichlet", n_iter=100000, idx_save=c(0,499,999), thin=1,
          n_save_P=100, alpha_fixed=F, verbose = 0, warmup=30000, slice_width=1, max_steps=100, estimates=estimates)
)

saveRDS(s2, "samples_SD.rds")
