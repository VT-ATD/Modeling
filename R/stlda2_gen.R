#!/usr/bin/Rscript
#  stlda2_gen.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.11.2018

# A different Spatio-Temporal LDA model

# Notes:
## Increasing Y to a large number seems to increase OMEGA's norm, that is, our limiting process blows up the Dirichlet hyperparam. Probably not desireable. Leads to numerical issues during exponentiation.

#' The Spatio-Temporal LDA Generative Model
#' 
#' Generate parameters and documents from the stLDA generative model. Documents are simply lists of integers between 1 and V.
#'
#' @param K Number of topics, a scalar integer.
#' @param M_mu Mean number of documents at each location-time combination
#' @param V The size of the vocabulary
#' @param N_mu Mean length of a document
#' @param L Number of spatial locations
#' @param Y Number of discrete time points
#' @param p Number of spatial dimensions
#' @param eta exchangible dirichlet prior topics, should be of length V and strictly positive.
#' @param sigma_om Anisotropic error variance for time component (kinda like a nugget).
#' @param sigma_theta Anisotropic random effect variance for topics of each doc
#' @param sigma_b Covariance parameter in Gaussian kernel for spatial closeness.
#' @param Lengthscale parameter in in Gaussian kernel for spatial closeness.
#' @return A list with named components giving each parameter.
gen_stlda <- function(K, M_mu, V, N_mu, L, Y, p = 2, eta = rep(1, V), sigma_om = 1, sigma_theta = 1, sigma_b = 1, l_kern = 1) {

    #TODO: Input checks

    # Get how many docs at each location
    Ms <- matrix(rpois(L*Y, M_mu), nrow = Y, ncol = L)

    # Generate spatial locations
    X <- matrix(rnorm(L*p), nrow = L, ncol = p)

    # Generate the OMEGA coefficient matrix
    B <- exp(-as.matrix(dist(X)^2 / l_kern))
    B <- sigma_b * B 

    # Generate OMEGAs for each time period
    OMEGAs <- lapply(1:Y, function(y) matrix(NA, nrow = L, ncol = K))

    # Generate initial OMEGA from iid normals
    OMEGAs[[1]] <- matrix(rnorm(L*(K-1)), nrow = L, ncol = K-1)

    # Step through other OMEGAs
    if (Y > 1) {
        for (y in 2:Y) {
            OMEGAs[[y]] <- B %*% OMEGAs[[y-1]] + matrix(rnorm(L*(K-1)), nrow = L, ncol = K-1)
        }
    }

    # Calculate ALPHAs
    ALPHAs <- lapply(OMEGAs, function(OMEGA) exp(OMEGA))

    # Topics
    PHI <- matrix(NA, nrow = K, ncol = V)
    for (k in 1:K) {
        PHI[k,] <- rdirich(eta)
    }

    # THETAs is a list of lists of numeric matrices whos rows live on a simplex.
    # The index of the first list gives the time, while that of the second gives the location.
    # The numeric matrix is of dimension Ms[time, location] by K, each row represents a document at this year location, while each column represents a topic.
    THETAs <- list()
    for (y in 1:Y) {
        THETAs[[y]] <- list()
        for (l in 1:L) {
            THETA <- matrix(NA, nrow = Ms[y,l], ncol = K)
            if (Ms[y,l] > 0) {
                for (m in 1:Ms[y,l]) {
                    #THETA[m,] <- rdirich(ALPHAs[[y]][l,])
                    ALPHA <- OMEGAs[[y]][l,] + matrix(rnorm((K-1), 0, sigma_theta))
                    THETA[m,] <- softmax(c(ALPHA, 0))
                }
            }
            THETAs[[y]][[l]] <- THETA
        }
    }

    # Similarly, docs is a list of lists, the outer representing time and the inner location, but this time, the innermost element is yet another list, each element of which is an integer vector representing a document.
    docs <- list()
    Ns <- list()
    for (y in 1:Y) {
        Ns[[y]] <- list()
        docs[[y]] <- list()
        for (l in 1:L) {
            docs[[y]][[l]] <- list()
            Ns[[y]][[l]] <- list()
            if (Ms[y,l] > 0) {
                for (m in 1:Ms[y,l]) {
                    N <- rpois(1, N_mu)
                    docs[[y]][[l]][[m]] <- rep(NA, N)
                    Ns[[y]][[l]][[m]] <- N
                    for (n in 1:N) {
                        z <- sample(1:K, 1, prob = THETAs[[y]][[l]][m,])
                        docs[[y]][[l]][[m]][n] <- sample(1:V, 1, prob = PHI[z,])
                    }
                }
            }
        }
    }

    return(list(X = X, Ms = Ms, OMEGAs = OMEGAs, PHI = PHI, THETAs = THETAs, 
                docs = docs, B = B, Ns = Ns))
}
