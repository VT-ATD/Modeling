#!/usr/bin/Rscript
#  dev/kalman_filt_mult.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.20.2018

## Kalman filter in multiple dimensions, isotropic error
require(mvtnorm)

sigma_x <- 0.1
sigma_z <- 0.1
phi <- 0.8
Y <- 70#Number of time points
P <- 2# Dimensionality of space.
Ms <- rpois(Y,15) + 1#Number of observations per time point
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
run_est <- function(a) {
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
n_log_lik <- function(a) run_est(a)$nll

set.seed(123)
res <- optim(as.numeric(A), n_log_lik, method = 'L-BFGS-B')
est <- matrix(res$par, ncol = P)
(est + t(est)) / 2
A
