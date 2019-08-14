suppressWarnings(suppressMessages(library("chron")))
#setwd("~/ATD Group Spring 2019")

tform <- function(vec) {
  sfvec <- tolower(vec)
  sfvec <- removeNumbers(sfvec)
  sfvec <- removeWords(sfvec, stopwords("english"))
  sfvec <- removePunctuation(sfvec)
  sfvec <- stripWhitespace(sfvec)
  
  return (sfvec)
}

convertStable <- function(textdist,eps=2e-308) {
  # Takes in a topics vector from the lda model
  indices <- which(textdist$beta==0)
  textdist$beta[indices] <- eps
  numtopics <- max(textdist$topic)
  # Now I need to renormalize
  for(i in 1:numtopics) {
    inds <- which(textdist$topic==i)
    textdist$beta[inds] <- textdist$beta[inds]/sum(textdist$beta[inds])
  }
  return (textdist)
}

comp.Dists <- function(d1,d2,eps=2e-308) {
  d1 <- convertStable(d1$topics)
  d2 <- convertStable(d2$topics)
  print(dim(d1))
  print(dim(d2))
  # These two are the distributions we are comparing
  setseq <- seq(1,max(d2$topic),by=1)
  
  w1 <- setdiff(unique(d2$term),unique(d1$term))
  if(identical(w1,character(0))) {
    enddf1 <- d1
  }
  else {
    add.df1 <- data.frame(expand.grid(setseq,w1),2e-308)
    colnames(add.df1) <- colnames(d1)
    enddf1 <- suppressWarnings(bind_rows(d1,add.df1))
  }
  w2 <- setdiff(unique(d1$term),unique(d2$term))
  if(identical(w2,character(0))) {
    enddf2 <- d2
  }
  else {
    add.df2 <- data.frame(expand.grid(setseq,w2),2e-308)
    colnames(add.df2) <- colnames(d2)
    enddf2 <- suppressWarnings(bind_rows(d2,add.df2)) # you get a coertion warning, not needed
  }
  
  #print(w2)
  # w1 is words in d1 but not in d2
  # w2 is words in d2 but not in d1
  return (list(d1=enddf1,d2=enddf2))
}