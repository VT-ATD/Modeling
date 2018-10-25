#!/usr/bin/Rscript
#  dev/tests.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.24.2018

## Test em_funcs.R


# Test PHI estimation
require(stlda)
source('R/stlda2_gen.R')


K <- 3
M_mu <- 200
V <- 3
N_mu <- 300
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

Zs <- e_step(THETAs, PHI)
norm(est_PHI(docs, Zs, THETAs) - PHI, '2') < 0.01

## Test the whole enchilada, initializing at the truth.
K <- 2
M_mu <- 200
V <- 2
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
attach(ret)

fit <- stlda_em(docs, K, V, THETAs, OMEGAs, PHI)

fit <- stlda_em(docs, K, V)
