#!/usr/bin/Rscript
#  stlda/R/em_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.24.2018

## The negative log posterior on a single alpha given omega and the E step's Z's.
quick_alpha_cost <- function(alpha, Z_sums, N, omega, sigma_theta, K) {
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

est_alpha <- function(Z_hat, omega, sigma_theta, K) {
    N <- nrow(Z_hat)
    Z_sums <- colSums(Z_hat)
    est <- optim(omega, quick_alpha_cost, Z_sums = Z_sums, N = N, omega = omega,
                 K = K, sigma_theta = sigma_theta, method = 'L-BFGS-B')
    #est <- optim(omega, alpha_cost, Z_hat = Z_hat, omega = omega,
    #             sigma_theta = sigma_theta)
    return(est$par)
}

est_ALPHAs <- function(Zs, OMEGAs, sigma_theta, Ms) {
    Y <- length(OMEGAs)
    K <- ncol(OMEGAs[[1]]) + 1
    L <- nrow(OMEGAs[[1]])

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
                    ALPHA[m,] <- est_alpha(Z_hat, omega, sigma_theta, K)
                }
                ALPHAs[[y]][[l]] <- ALPHA
            }
        }
    }

    return(ALPHAs)
}

e_step <- function(THETAs, PHI, Ms, Ns) {
    K <- nrow(PHI)
    Y <- nrow(Ms)
    L <- ncol(Ms)

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
                            w <- docs[[y]][[l]][[m]][n]
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

est_PHI <- function(docs, Zs, THETAs_hat, Ms, Ns, K, eta) {
    Y <- nrow(Ms)
    L <- ncol(Ms)

    PHI_est <- matrix(rep(eta, K), nrow = K, byrow = TRUE)
    for (y in 1:Y) {
        for (l in 1:L) {
            if (Ms[y,l] > 0) {
                for (m in 1:Ms[y,l]) {
                    for (n in 1:Ns[[y]][[l]][[m]]) {
                        w <- docs[[y]][[l]][[m]][n]
                        for (k in 1:K) {
                            PHI_est[k,w] <- PHI_est[k,w] + Zs[[y]][[l]][[m]][n,k] * THETAs_hat[[y]][[l]][m,k]
                        }
                    }
                }
            }
        }
    }
    PHI_est <- PHI_est / rowSums(PHI_est)
    return(PHI_est)
}

#' SpatioTemporal LDA estimation via Expecation Maximization
#'
#'
#' @param docs docs
#' @param K K
#' @param V V
#' @param THETAs_init Initial value for document-topic prevalences. Initial values are important due to the nonconvex nature of the posterior.
#' @param OMEGAs_init Initial value for space-time-topic in R^d.
#' @param PHI_init Initial value for topic-word matrix
#' @param max_iters Max outer EM iterations
#' @param thresh Mean l2 norm on OMEGA matrices to stop iteration.
#' @param eta A vector of length V, the prior on PHI
#' @param sigma_om The variance for OMEGA, for now fixed
#' @param sigma_theta Random effect variance for each document
#' @param sigma_b Covariance for spatial kernel.
#' @param l_kern Lengthscale for spatial kernel.
#' @export
stlda_em <- function(docs, K, vocab, THETAs_init, OMEGAs_init, PHI_init,
                     max_iters = 1000, thresh = 1e-3, eta = rep(1, V), sigma_om = 1,
                     sigma_theta = 1, sigma_b = 1, l_kern = 1) {

    # A utility function to transition between the d-1 simplex and R^d, assuming the last element to be 0.
    V <- length(vocab)
    ALPHAs_2_THETAs <- function(ALPHAs) {
        lapply(ALPHAs, function(ALPHAY) lapply(ALPHAY, function(ALPHA) {

                   t(apply(ALPHA, 1, function(x) softmax(c(x,0))))
        }))
    }

    # Get documents in each location-time.
    Y <- length(docs)
    L <- length(docs[[1]])
    Ms <- matrix(NA, nrow = Y, ncol = L)
    for (y in 1:Y) {
        for (l in 1:L) {
            Ms[y,l] <- length(docs[[y]][[l]])
        }
    }

    # Compute document lengths
    Ns <- list()
    for (y in 1:Y) {
        Ns[[y]] <- list()
        for (l in 1:L) {
            Ns[[y]][[l]] <- list()
            for (m in 1:Ms[y,l]) {
                Ns[[y]][[l]][[m]] <- length(docs[[y]][[l]][[m]])
            }
        }
    }

    ## Initialize parameters
    if (missing(THETAs_init)) {
        THETAs_hat <- list()
        for (y in 1:Y) {
            THETAs_hat[[y]] <- list()
            for (l in 1:L) {
                if (Ms[y,l] > 0) {
                    THETA <- matrix(rgamma(Ms[y,l] * K, 1, 1), nrow = Ms[y,l], ncol = K)
                    THETA <- THETA / rowSums(THETA)
                } else {
                    THETA <- matrix(NA, nrow = Ms[y,l], ncol = K)
                }
                THETAs_hat[[y]][[l]] <- THETA
            }
        }
    } else {
        THETAs_hat <- THETAs_init
    }
    if (missing(OMEGAs_init)) {
        OMEGAs_hat <- lapply(1:Y, function(y) matrix(rnorm(L*(K-1)), nrow = L))
    } else {
        OMEGAs_hat <- OMEGAs_init
    }
    if (missing(PHI_init)) {
        PHI_hat <- matrix(rgamma(K*V, 1,1), nrow = K)
        PHI_hat <- PHI_hat / rowSums(PHI_hat)
    } else {
        PHI_hat <- PHI_init
    }

    # An EM Algorithm for maximum a posteriori inference.
    diff <- Inf
    iter <- 1
    while (diff > thresh && iter < max_iters) {
        print(iter)
        print(diff)

        ## E Step
        iter <- iter + 1
        Zs <- e_step(THETAs_hat, PHI_hat, Ms, Ns)

        ## M Step
        # Topic-word matrix (closed form update)
        PHI_hat <- est_PHI(docs, Zs, THETAs_hat, Ms, Ns, K, eta)

        # Document-topic prevalences
        ALPHAs_hat <- est_ALPHAs(Zs, OMEGAs_hat, sigma_theta, Ms)
        THETAs_hat <- ALPHAs_2_THETAs(ALPHAs_hat)

        # Location-topic prevalences
        OMEGAs_new <- list()
        for (y in 1:Y) {
            OMEGA <- matrix(NA, nrow = L, ncol = K-1)
            for (l in 1:L) {
                OMEGA[l,] <- colMeans(ALPHAs_hat[[y]][[l]])
            }
            OMEGAs_new[[y]] <- OMEGA
        }

        ## Convergence Checking via average l2 on Omegas
        diffs <- sapply(1:Y, function(y) norm(OMEGAs_new[[y]] - OMEGAs_hat[[y]]))
        diff <- mean(diffs)
        OMEGAs_hat <- OMEGAs_new
    }

    # Calculate most prevalent topics
    top_pop <- rep(0, K)
    for (y in 1:Y) {
        for (l in 1:L) {
                top_pop <- top_pop + colSums(THETAs_hat[[y]][[l]])
        }
    }

    # Create return object
    ret <- list(THETAs_hat = THETAs_hat, OMEGAs_hat = OMEGAs_hat, PHI_hat = PHI_hat)
    ret$vocab = vocab
    ret$Y = Y
    ret$L = L
    ret$Ms = Ms
    ret$Ns = Ns
    ret$V = V
    ret$K = K
    ret$top_pop <- top_pop
    class(ret) <- "stlda"
    return(ret)
}

#' Summary of an stLDA model fit object
#' 
#'
#'
#' @param fit An stlda object resulting from calling stlda_em.
#' @param top_topics How many topics to display.
#' @param n_top_words How many words to display from each topic.
#' @return No return value; used for side effect of printing information to console.
#' @export
summary.stlda <- function(fit, n_top_topics = 3, n_top_words = 5) {
    top_topics <- (1:fit$K)[order(-fit$top_pop)][1:n_top_topics]

    # Calculate most important words for each topic
    top_words <- t(sapply(top_topics, function(k) vocab[order(-fit$PHI_hat[k,])][1:n_top_words]))

    cat(paste(c(rep('_', Sys.getenv("COLUMNS")), '\n'), collapse = ''))
    cat("Spatio-Temporal Topic Model Summary\n")
    cat(paste(c(rep('-', Sys.getenv("COLUMNS")), '\n'), collapse = ''))
    cat("Corpus over ", fit$Y, " time points and ", fit$L, " locations with ", fit$K, " topics.\n", sep = '')
    cat(paste(sum(fit$Ms), " total documents with ", fit$V, " unique words and average document size of ", round(mean(unlist(fit$Ns)), 1), ".\n", sep = ""))
    cat(paste(c(rep('-', Sys.getenv("COLUMNS")), '\n'), collapse = ''))
    cat(paste("Top ", n_top_words, " words of ", n_top_topics, " Most Prevalent Topics:\n", sep = ""))
    for (ki in 1:n_top_topics) {
        cat(paste("Topic ", top_topics[ki], ": ", paste(top_words[ki,], collapse = ", "), "\n", sep = ""))
    }
    cat(paste(c(rep('_', Sys.getenv("COLUMNS")), '\n'), collapse = ''))
}

#' Print stlda Object
#'
#'
#' @param fit An stlda object.
#' @export
print.stlda <- function(fit) {
    cat("stlda object with", fit$Y, "time points,", fit$L, "locations and", fit$K, "topics.\n")
}
