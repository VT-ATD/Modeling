#!/usr/bin/Rscript
#  R/em_algo.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.11.2018

## Perform naive inference on the sttopic model
source('Modeling/R/stlda2_gen.R')
source('R/lib.R')
source('R/em_funcs.R')

K <- 3
M_mu <- 200
V <- 3
N_mu <- 3000
L <- 1
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

ALPHAs_2_THETAs <- function(ALPHAs) {
    lapply(ALPHAs, function(ALPHAY) lapply(ALPHAY, function(ALPHA) {
               t(apply(ALPHA, 1, function(x) softmax(c(x,0))))
    }))
}

# Do our EM algo
max_iters <- 300
diff <- Inf
thresh <- 1e-3

# A random init
THETAs_hat <- lapply(THETAs, function(THETAY) lapply(THETAY, function(theta) {
        matrix(rgamma(nrow(theta)*ncol(theta), 1, 1), nrow = nrow(theta), ncol =ncol(theta))
}))
OMEGAs_hat <- lapply(1:Y, function(y) matrix(rnorm(L*(K-1)), nrow = L))
# Init at truth
#THETAs_hat <- THETAs
#OMEGAs_hat <- OMEGAs

iter <- 1
while (diff > thresh && iter < max_iters) {
    print(iter)
    print(diff)
    
    iter <- iter + 1
    Zs <- e_step(THETAs_hat, PHI)

    ALPHAs_hat <- est_ALPHAs(Zs, OMEGAs_hat, sigma_theta)
    THETAs_hat <- ALPHAs_2_THETAs(ALPHAs_hat)

    OMEGAs_new <- list()
    for (y in 1:Y) {
        OMEGA <- matrix(NA, nrow = L, ncol = K-1)
        for (l in 1:L) {
            OMEGA[l,] <- colMeans(ALPHAs_hat[[y]][[l]])
        }
        OMEGAs_new[[y]] <- OMEGA
    }

    diffs <- sapply(1:Y, function(y) norm(OMEGAs_new[[y]] - OMEGAs_hat[[y]]))
    diff <- mean(diffs)
    OMEGAs_hat <- OMEGAs_new
}
