#!/usr/bin/Rscript
#  R/em_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.11.2018

## The negative log posterior on a single alpha given omega and the E step's Z's.
quick_alpha_cost <- function(alpha, Z_sums, N, omega, sigma_theta) {
    alpha <- c(alpha, 0)
    omega <- c(omega, 0)
    denom <- log(sum(exp(alpha)))
    logprior <- -1/(2*sigma_theta) * sum((alpha - omega)^2)

    loglik <- 0
    for (k in 1:K) {
        loglik <- loglik + alpha[k] * Z_sums[k] 
    }
    loglik <- loglik - N * denom

    return(-(loglik + logprior))
}

est_alpha <- function(Z_hat, omega, sigma_theta) {
    N <- nrow(Z_hat)
    Z_sums <- colSums(Z_hat)
    est <- optim(omega, quick_alpha_cost, Z_sums = Z_sums, N = N, omega = omega,
                 sigma_theta = sigma_theta, method = 'L-BFGS-B')
    #est <- optim(omega, alpha_cost, Z_hat = Z_hat, omega = omega, 
    #             sigma_theta = sigma_theta)
    return(est$par)
}

est_ALPHAs <- function(Zs, OMEGAs, sigma_theta) {
    ALPHAs <- list()
    for (y in 1:Y) {
        ALPHAs[[y]] <- list()
        for (l in 1:L) {
            ALPHAs[[y]][[l]] <- list()
            if (Ms[y,l] > 0) {
                ALPHA <- matrix(NA, nrow = Ms[y,l], ncol = K-1)
                for (m in 1:Ms[y,l]) {
                    Z_hat <- Zs[[y]][[l]][[m]]
                    omega <- OMEGAs[[y]][l,]
                    ALPHA[m,] <- est_alpha(Z_hat, omega, sigma_theta)
                }
                ALPHAs[[y]][[l]] <- ALPHA
            }
        }
    }

    return(ALPHAs)
}

e_step <- function(THETAs, PHI) {
    Zs <- list()# A list of length Y, each sublist of length L, each sublist of that of length Ms[y,l], finally, containing a matrix with as many rows as that doc has words, and a column for each topic.

    # Go across every document, get the expectation of each word's topic.
    for (y in 1:Y) {
        Zs[[y]] <- list()
        for (l in 1:L) {
            Zs[[y]][[l]] <- list()
            if (Ms[y,l] > 0) {

                for (m in 1:Ms[y,l]) {
                    N <- Ns[[y]][[l]][[m]]
                    Z <- matrix(NA, nrow = N, ncol = K)
                    if (N > 0) {

                        Z <- matrix(NA, nrow = N, ncol = K)
                        for (n in 1:N) {
                            w <- docs[[y]][[l]][[m]][[n]]
                            for (k in 1:K) {
                                Z[n,k] <- THETAs[[y]][[l]][m,k] * PHI[k,w]
                            }
                        }
                        Z <- Z / rowSums(Z)
                    }
                    Zs[[y]][[l]][[m]] <- Z
                }
            }
        }
    }

    return(Zs)
}
