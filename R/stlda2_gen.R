#!/usr/bin/Rscript
#  stlda2_gen.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.11.2018

# A different Spatio-Temporal LDA model

# Notes:
## Increasing Y to a large number seems to increase OMEGA's norm, that is, our limiting process blows up the Dirichlet hyperparam. Probably not desireable. Leads to numerical issues during exponentiation.

rdirich <- function(alpha) {
    x <- rgamma(length(alpha), alpha, 1)
    x / sum(x)
}

# Corpus Characteristics
K <- 5
M_mu <- 2 #Mean documents for each time - location combination
V <- 3
N_mu <- 20 #Mean words per document
L <- 3 # Locations
Y <- 4 #Time points (discrete time) TODO: Algo will break if T < 2
p <- 2 #Spatial dimensions

# Hyperparams
sigma_om <- 1 # Error variance in OMEGA AR process
l_kern <- 1 # Length scale in kernel process
eta <- rep(1, V)

# Get how many docs at each location
Ms <- matrix(rpois(L*Y, M_mu), nrow = Y, ncol = L)

# Generate spatial locations
X <- matrix(rnorm(L*p), nrow = L, ncol = p)

# Generate the OMEGA coefficient matrix
B <- exp(-as.matrix(dist(X)^2 / l_kern))
B <- B / sum(diag(B))

# Generate OMEGAs for each time period
OMEGAs <- lapply(1:Y, function(y) matrix(NA, nrow = L, ncol = K))

# Generate initial OMEGA from iid normals
OMEGAs[[1]] <- matrix(rnorm(L*K), nrow = L, ncol = K)

# Step through other OMEGAs
for (y in 2:Y) {
    OMEGAs[[y]] <- B %*% OMEGAs[[y-1]] + matrix(rnorm(L*K), nrow = L, ncol = K)
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
                THETA[m,] <- rdirich(ALPHAs[[y]][l,])
            }
        }
        THETAs[[y]][[l]] <- THETA
    }
}

# Similarly, docs is a list of lists, the outer representing time and the inner location, but this time, the innermost element is yet another list, each element of which is an integer vector representing a document.
docs <- list()
for (y in 1:Y) {
    docs[[y]] <- list()
    for (l in 1:L) {
        docs[[y]][[l]] <- list()
        if (Ms[y,l] > 0) {
            for (m in 1:Ms[y,l]) {
                N <- rpois(1, N_mu)
                docs[[y]][[l]][[m]] <- rep(NA, N)
                for (n in 1:N) {
                    z <- sample(1:K, 1, prob = THETAs[[y]][[l]][m,])
                    docs[[y]][[l]][[m]][n] <- sample(1:V, 1, prob = PHI[z,])
                }
            }
        }
    }
}
