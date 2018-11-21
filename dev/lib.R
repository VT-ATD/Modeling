#!/usr/bin/Rscript
#  R/lib.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.02.2018

#A numerically stable softmax function
softmax <- function(x) {
    exp(x - max(x) - log(sum(exp(x - max(x)))))
}

rdirich <- function(alpha) {
    x <- rgamma(length(alpha), alpha, 1)
    x / sum(x)
}

#An inverse softmax function which assumes that the last argument of the softmax was zero.
inv_softmax <- function(rho) {
    a <- log(rho)
    a <- a - log(rho[length(rho)])
    return(a)
}

