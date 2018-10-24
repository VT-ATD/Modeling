#!/usr/bin/Rscript
#  R/naive_inf.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.24.2018

## Perform naive inference on the sttopic model
source('Modeling/R/stlda2_gen.R')
source('R/lib.R')

K <- 3
M_mu <- 100
V <- 3
N_mu <- 1000
L <- 3
Y <- 2

p = 2
eta = rep(1, V)
sigma_om = 1
sigma_b = 1
l_kern = 1
sigma_theta = 1

ret <- gen_stlda(K, M_mu, V, N_mu, L, Y, p = p, eta = eta, sigma_om = sigma_om, sigma_theta = sigma_theta, sigma_b = sigma_b, l_kern = l_kern)

# Inference step 1: Assume all but THETA is known
attach(ret)

## Random inits
ALPHAs_hat <- list()
for (y in 1:Y) {
    ALPHAs_hat[[y]] <- list()
    for (l in 1:L) {
        ALPHA <- matrix(NA, nrow = Ms[y,l], ncol = K-1)
        if (Ms[y,l] > 0) {
            for (m in 1:Ms[y,l]) {
                ALPHA[m,] <- rnorm((K-1))
            }
        }
        ALPHAs_hat[[y]][[l]] <- ALPHA
    }
}

OMEGAs_hat <- list()
for (y in 1:Y) {
    OMEGAs_hat[[y]] <- matrix(rnorm(L*(K-1)), nrow = L, ncol = K-1)
}

sigma_theta_hat <- 1
sigma_om_hat <- 1
######

#' Unfolds parameters into a vector for use with built in optimization routines
#' vec goes like this:
#' vec = [(ALPHAs by year, then by location, then by column), (OMEGAs by year then by column)]
pars2vec <- function(ALPHAs, OMEGAs, sigma_theta, sigma_om_hat) {
    c(as.numeric(unlist(unlist(ALPHAs))), 
      as.numeric(unlist(OMEGAs)),
      sigma_theta,
      sigma_om_hat)
}

#' Bring us back to our original representation
vec2pars <- function(vec, Ms) {
    used_params <- 0

    # Chunk up the ALPHAs
    ALPHAs_hat <- list()
    for (y in 1:Y) {
        ALPHAs_hat[[y]] <- list()
        for (l in 1:L) {
            if (Ms[y,l] > 0) {
                from <- used_params + 1
                to <- used_params + Ms[y,l] * (K-1)
                ALPHAs_hat[[y]][[l]] <- matrix(vec[from:to], nrow = Ms[y,l], ncol = K-1)
                used_params <- to
            } else {
                ALPHAs_hat[[y]][[l]] <- matrix(NA, nrow = 0, 0)
            }
        }
    }

    # OMEGAs
    OMEGAs <- list()
    for (y in 1:Y) {
        from <- used_params + 1
        to <- used_params + L * (K-1)
        OMEGAs[[y]] <- matrix(vec[from:to], nrow = L, ncol = K-1)
        used_params <- to
    }

    # Random Effects Variance
    sigma_theta <- vec[used_params + 1]
    used_params <- used_params + 1
    sigma_om <- vec[used_params + 1]
    used_params <- used_params + 1

    return(list(ALPHAs_hat = ALPHAs_hat, OMEGAs_hat = OMEGAs, sigma_theta_hat = sigma_theta, 
                sigma_om_hat = sigma_om))
}

ff <- vec2pars(pars2vec(ALPHAs_hat, OMEGAs_hat, sigma_theta_hat, sigma_om_hat), Ms)
ff$ALPHAs_hat
ALPHAs_hat
ff$OMEGAs_hat
OMEGAs_hat
ff$sigma_theta_hat
sigma_theta_hat
ff$sigma_om_hat
sigma_om_hat

#' The negative log posterior of the data
nlp <- function(ALPHAs_hat, OMEGAs_hat, sigma_theta_hat, sigma_om_hat) {
    lp <- 0

    ## Do the OMEGA AR process
    # Omega 0
    for (l in 1:L) {
        for (k in 1:(K-1)) {
            lp <- lp + dnorm(OMEGAs_hat[[1]][l,k], 0, sigma_om)
        }
    }
    ## Omegas afterwards
    BOMEGAs <- lapply(OMEGAs_hat, function(OMEGA) B %*% OMEGA)
    if (Y > 1) {
        for (y in 2:Y) {
            for (l in 1:L) {
                for (k in 1:(K-1)) {
                    lp <- lp + dnorm(OMEGAs_hat[[y]][l,k], BOMEGAs[[y-1]][l,k], sigma_om)
                }
            }
        }
    }

    ## Generate ALPHAs from OMEGA and softmax to get THETAs
    THETAs_hat <- list()
    for (y in 1:Y) {
        THETAs_hat[[y]] <- list()
        for (l in 1:L) {
            THETAs_hat[[y]][[l]] <- matrix(NA, nrow = Ms[y,l], ncol = K)
            if (Ms[y,l] > 0) {
                for (m in 1:Ms[y, l]) {
                    for (k in 1:(K-1)) {
                        lp <- lp + dnorm(ALPHAs_hat[[y]][[l]][m,k], OMEGAs_hat[[y]][l,k], sigma_theta_hat, log = TRUE)
                    }
                    THETAs_hat[[y]][[l]][m,] <- softmax(c(ALPHAs_hat[[y]][[l]][m,], 0))
                }
            }
        }
    }

    # Likelihood
    for (y in 1:Y) {
        for (l in 1:L) {
            if (Ms[y,l] > 0) {
                for (m in 1:Ms[y, l]) {
                    doc <- docs[[y]][[l]][[m]]
                    N <- length(doc)
                    for (n in 1:N) {
                        ll <- 0
                        for (k in 1:K) {
                            ll <- ll + THETAs_hat[[y]][[l]][m,k] * PHI[k,doc[n]]
                        }
                        lp <- lp + log(ll)
                    }
                }
            }
        }
    }

    return(-lp)
}

#ALPHAs_true <- lapply(THETAs, function(THETAY) lapply(THETAY, function(THETA) t(apply(THETA, 1, function(x) inv_softmax(x)[1:(K-1)]))))
#nlp(ALPHAs_hat, OMEGAs_hat)

init <- pars2vec(ALPHAs_hat, OMEGAs_hat, sigma_theta_hat, sigma_om_hat)
#init <- pars2vec(ALPHAs_true, OMEGAs)
cost_wrap <- function(vec) {
    do.call(nlp, vec2pars(vec, Ms))
}

minvar <- 1e-4
lower <- c(rep(-Inf, length(init)-2), minvar, minvar)
upper <- rep(Inf, length(init))
system.time(ret <- optim(init, cost_wrap, method = 'L-BFGS-B', lower = lower, upper = upper))
print(ret)

# Check for bound collisions
if (all.equal(ret$par[last_2_inds], c(minvar, minvar))) {
    warning("Lower bound reached in variance optimization.")
}

### Eyeball norm:
# See how we did with ALPHAs
ALPHA_est <- vec2pars(ret$par, Ms)$ALPHAs_hat
cat('est:\n')
lapply(ALPHA_est, function(ALPHAY) lapply(ALPHAY, function(ALPHA) t(apply(ALPHA, 1, function(row) softmax(c(row, 0))))))[[1]][[1]]
#print(t(apply(ALPHA_est[[1]][[1]], 1, function(row) softmax(c(row, 0)))))
cat('true:\n')
print(THETAs[[1]][[1]])

# See how we did with OMEGAs
OMEGA_est <- vec2pars(ret$par, Ms)$OMEGAs_hat
OMEGA_est[[1]]
OMEGAs[[1]]
