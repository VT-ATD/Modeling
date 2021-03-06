#!/usr/bin/Rscript
#  dev/kalman_filt_again.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.08.2018

## Get Kalman Filtering up and running for later integration with the EM algo.

#TODO: We may be assuming mean zero.
###### For starters, 1D latent process and observed data together with multiple observations per time point.
set.seed(12345)
# Params
sigma_x <- 0.5
sigma_z <- 0.1
phi <- 0.8
Y <- 50#Number of time points
Ms <- rpois(Y,50) + 1#Number of observations per time point
mu <- 0

## Simulate
# Do latent process
x_0 <- rnorm(1, 0, sigma_x)
x <- rep(NA, Y)
x <- c(x_0, x)
for (y in 1:Y) {
    x[y+1] <- phi * x[y] + rnorm(1, 0, sigma_x) 
}

# Observe data
z <- list()
for (y in 1:Y) {
    z[[y]] <- rnorm(Ms[y], mu + x[y+1], sigma_z)
}

## E step/likelihood
run_est <- function(phi, sigma_z, sigma_x) {
    # Storage
    x_hats <- c(0, rep(NA, Y))
    x_sig2s <- c(sigma_x^2, rep(NA, Y))

    log_lik <- 0
    for (y in 1:Y) {
        x_hats_var <- (phi^2 * x_sig2s[y] + sigma_x^2)
        z_bar_var <- (sigma_z^2 / Ms[y])

        x_sig2s[y+1] <- 1/(1/x_hats_var + 1/z_bar_var)
        x_hats[y+1] <- x_sig2s[y+1] * (phi*x_hats[y] / x_hats_var  + mean(z[[y]] - mu) / z_bar_var)
        for (m in 1:Ms[y]) {
            log_lik <- log_lik + dnorm(z[[y]][m], x_hats[y+1] + mu, sqrt(x_sig2s[y+1] + sigma_z^2), log = TRUE)
        }
    }

    return(list(ll = -log_lik, x_hats = x_hats))
}
n_log_lik <- function(params) run_est(params[1], params[2], params[3])$ll

minvar <- 1e-5
optim(c(phi, sigma_z, sigma_x),  n_log_lik, method = "L-BFGS-B", 
      lower = c(0, minvar, minvar), upper = c(1, Inf, Inf))

phis <- seq(0, 1, length.out = 500)
liks <- sapply(phis, n_log_lik)
plot(phis, liks)

x_hats <- run_est(phi, sigma_z, sigma_x)$x_hats
all <- c(x, x_hats)
plot(NA, NA, xlim = c(1, Y), ylim = c(min(all), max(all)))
points(x_hats, col = 'red')
points(x, col = 'blue')
points(2:(Y+1), sapply(z, mean), col = 'green')
