library(cudarpackage)

#simulate data
sim_data <- function(G, X, group_n, mu_beta, sd_beta, corr=NULL, mu_dispersion, sd_dispersion){
  V <- length(group_n)
  group_id <- unlist(sapply(1:V, function(g) rep(g, group_n[g])))
  X <- X[group_id,]
  if(!is.null(corr)){
    ch <- chol(corpcor::rebuild.cov(corr, sd_beta^2))
    beta <- sapply(1:G, function(g) t(rnorm(V)) %*% ch + mu_beta)
  } else{
    beta <- sapply(1:G, function(g) rnorm(V, mu_beta, sd_beta))
  }
}  
M
G <- 1000
group_n <- c(3, 3, 4)
V <- length(group_n)
group_id <- unlist(sapply(1:V, function(g) rep(g, group_n[g])))
X <- matrix(c(1, -1, 0,
              1,  1, 0,
              1,  0, 1), V, V,
              byrow=T)[group_id,]

mu_beta <- c(4, 0, 0)
sd_beta <- c(1.5, .2, .4)
mu_disper <- function(beta) pmin(-beta[1,]^.25, 1)
sd_disper <- 1

beta <- sapply(1:G, function(g) rnorm(V, mu_beta, sd_beta))
disper <- rlnorm(G, mu_disper(beta), sd_disper)

y <- matrix(rnbinom(G*sum(group_n), mu = exp(X %*% beta), size = 1/disper),
            G, sum(group_n), byrow=T)

ord <- order(beta[1,])
head(cbind(y[ord,], beta1 = beta[1,ord], beta2 = beta[2,ord], beta3 = beta[3,ord], disper = disper[ord]))

#fit model (no voom)




