#!/usr/bin/Rscript
#  dev/kalman_filt_mult.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.20.2018

## Kalman filter in multiple dimensions, isotropic error
require(mvtnorm)

sigma_x <- 0.5
sigma_z <- 0.1
phi <- 0.8
Y <- 100#Number of time points
P <- 2# Dimensionality of space.
Ms <- rpois(Y,100) + 1#Number of observations per time point
mu <- 0

# Transition matrix A
set.seed(123)
A <- matrix(rnorm(P*P), ncol = P)
A <- A %*% t(A)
A <- 0.8 * A / max(eigen(A)$values)

# Carry out latent state
X <- matrix(NA, nrow = Y+1, ncol = P)
X[1,] <- rnorm(P, 0, sigma_x)
for (y in 1:Y) {
    X[y+1,] <- A %*% X[y,] + rnorm(P, 0, sigma_x)
}
X <- X[-1,]

# Observe data
Z <- list()
for (y in 1:Y) {
    Z[[y]] <- matrix(rnorm(Ms[y]*P, 0, sigma_z), nrow = Ms[y], ncol = P)
    Z[[y]] <- t(apply(Z[[y]], 1, function(z) z + X[y,]))
}

## E step/likelihood
#run_est <- function(sigma_x) {
run_est <- function(a, sigma_z, sigma_x) {
    A <- matrix(a, nrow = P)
    # Storage
    X_hats <- matrix(NA, nrow = Y+1, ncol = P)
    X_hats[1,] <- rep(0,P)
    X_sig2s <- array(NA, dim = c(P, P, Y+1))
    X_sig2s[,,1] <- diag(sigma_x^2, P)

    log_lik <- 0
    for (y in 1:Y) {
        # TODO: only store inverses.
        X_hats_var <- solve(A %*% X_sig2s[,,y] %*% t(A) + diag(sigma_x^2, P))
        Z_bar_var <- diag(Ms[y] / sigma_z^2, P)

        X_sig2s[,,y+1] <- solve(X_hats_var + Z_bar_var)
        X_hats[y+1,] <- X_sig2s[,,y+1] %*% (X_hats_var %*% A %*% X_hats[y,] + 
                                            Z_bar_var %*% colMeans(Z[[y]]))
        for (m in 1:Ms[y]) {
            log_lik <- log_lik + dmvnorm(Z[[y]][m,], X_hats[y+1,], 
                                 X_sig2s[,,y+1] + diag(sigma_z^2, P), log = TRUE)
        }
    }
    X_hats <- X_hats[-1,]
    return(list(nll = -log_lik, X_hats = X_hats))
}
#n_log_lik <- function(param) run_est(param)$nll
n_log_lik <- function(params) run_est(params[1:P^2], params[P^2+1], params[P^2+2])$nll

set.seed(123)
minvar <- 1e-6
res <- optim(c(as.numeric(A), sigma_z, sigma_x), n_log_lik, method = 'L-BFGS-B',
             lower = c(rep(-Inf,P^2), minvar, minvar), 
             upper = c(rep(Inf, P^2), Inf, Inf))
est <- matrix(res$par[1:P^2], ncol = P)
(est + t(est)) / 2
A
res$par[P^2+1]
res$par[P^2+2]

xs <- seq(0.001, 1, length.out = 1e2)
plot(xs, sapply(xs, n_log_lik))
