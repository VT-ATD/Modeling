#!/usr/bin/Rscript
#  dev/1k_app.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.24.2018

## Run our model on the 1 thousand lines.
require(stlda)

# Load data
load("~/atd_data/1klines.RData")

# Fit model
system.time(fit <- stlda_em(docs = docs, K = 6, vocab = vocab, max_iters = 10))

summary(fit)
