#!/usr/bin/Rscript
#  dev/proc_dat.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.24.2018

## Turn JSON data into something the model can read.
require(ndjson)
dat <- as.data.frame(stream_in("~/atd_data/1000Lines_decahose-2015-10-01-00-12-24_us_eng_prcsd"))

# Filter data (hopefully will be a preprocessing step someday)
dat <- dat[!is.na(dat$gnip.profileLocations.0.address.region),]# Nonnull region
dat <- dat[sapply(strsplit(dat$prcsd_body, ' '), length) > 0,]# Nonempty body

text <- strsplit(dat$prcsd_body, ' ')
vocab <- unique(unlist(text))
V <- length(vocab)

# Store docs as a triple nested list giving time/location/doc/<word vector>
inds <- lapply(text, function(tt) match(tt, vocab))

locations <- unique(dat$gnip.profileLocations.0.address.region)
docs <- list()# Only 1 time point for now.
y <- 1
docs[[y]] <- list()
for (l in 1:length(locations)) {
    docs[[y]][[l]] <- inds[dat$gnip.profileLocations.0.address.region == locations[l]]
}

save(docs, vocab, file = "~/atd_data/1klines.RData")
