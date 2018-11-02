#!/usr/bin/Rscript
#  dev/kalman_filt.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.31.2018

## Get Kalman Filtering up and running for later integration with the EM algo.

#TODO: We may be assuming mean zero.
###### For starters, 1D latent process and observed data together with multiple observations per time point.
set.seed(123)
# Params
sigma_x <- 0.1
sigma_z <- 0.5
phi <- 0.8
Y <- 5000#Number of time points
Ms <- rpois(Y,300) + 1#Number of observations per time point

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
    z[[y]] <- rnorm(Ms[y], x[y], sigma_z)
}

## E step/likelihood
e_loglik <- function(z, Ms, phi, sigma_x, sigma_z) {
    # Storage
    x_hats <- c(0, rep(NA, Y))
    x_sig2s <- c(sigma_x^2, rep(NA, Y))

    log_lik <- 0
    for (y in 1:Y) {
        x_sig2s[y+1] <- 1/(1/(phi^2 * x_sig2s[y] + sigma_x^2) + Ms[y]/sigma_z^2)
        x_hats[y+1] <- x_sig2s[y+1] * (phi*x_hats[y] / (phi^2 * x_sig2s[y] + sigma_x^2) + mean(z[[y]]) / (sigma_z^2 / Ms[y]))
        for (m in 1:Ms[y]) {
            log_lik <- log_lik + dnorm(z[[y]][m], x_hats[y+1], sqrt(x_sig2s[y+1] + sigma_z^2), log = TRUE)
        }
    }

    return(list(log_lik, x_hats))
}

x_hats <- e_loglik(z, Ms, phi, sigma_x, sigma_z)[[2]]


all <- c(x, x_hats)
plot(NA, NA, xlim = c(1, Y), ylim = c(min(all), max(all)))
points(x_hats, col = 'red')
points(x, col = 'blue')
points(sapply(z, mean), col = 'green')

## Max Lik Estimation 
cost <- function(params) {
    phi <- params[1]
    sigma_x <- params[2]
    sigma_z <- params[3]

    return(-e_loglik(z, Ms, phi, sigma_x, sigma_z)[[1]])
}

min_sd <- 1e-4
optim(c(phi, sigma_x, sigma_z), cost, method = 'L-BFGS-B', lower = c(0,min_sd,min_sd))
