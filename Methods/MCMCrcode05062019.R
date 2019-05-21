### MCMC code
# Shane Bookhultz
# Will implement the Gibbs, Collapsed Gibbs, Variational EM, and Variational EM algorithms
# To compare and contrast to one another. 

eps <- sqrt(.Machine$double.eps)


get.TF <- function(x, type=c("text", "url", "file"), 
                   lang="english", excludeWords=NULL, 
                   textStemming=FALSE,  colorPalette="Dark2",
                   min.freq=3, max.words=200)
{ 
  library("tm")
  library("SnowballC")
  library("wordcloud")
  library("RColorBrewer") 
  
  if(type[1]=="file") text <- readLines(x)
  else if(type[1]=="url") text <- html_to_text(x)
  else if(type[1]=="text") text <- x
  
  # Load the text as a corpus
  docs <- VCorpus(VectorSource(text))
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove stopwords for the language 
  docs <- tm_map(docs, removeWords, stopwords(lang))
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  # Remove your own stopwords
  if(!is.null(excludeWords)) 
    docs <- tm_map(docs, removeWords, excludeWords)
  
  # Text stemming
  if(textStemming) docs <- tm_map(docs, stemDocument)
  # Create term-document matrix
  tdm <- TermDocumentMatrix(docs)
  dtm <- DocumentTermMatrix(docs)
  m <- as.matrix(tdm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  # check the color palette name 
  # if(!colorPalette %in% rownames(brewer.pal.info)) colors = colorPalette
  # else colors = brewer.pal(8, colorPalette) 
  # # Plot the word cloud
  # set.seed(1234)
  # wordcloud(d$word,d$freq, min.freq=min.freq, max.words=max.words,
  #           random.order=FALSE, rot.per=0.35, 
  #           use.r.layout=FALSE, colors=colors)
  # 
  # invisible(list(tdm=tdm, freqTable = d))
  return (list(tdm,dtm))
}



##############
# Libraries
##############
library(DirichletReg)
# Will use document topic matrix

# So what I'll do is convert to TF then go through this

KLdiverge <- function(d1,d2) {
  return (-sum(d1*log(d2/d1)))
} # Not even sure if I need this

normalize <- function(vec) {
  vec/sum(vec)
}
# Do not pay attention to this function
#gen_lda <- function(K,D,topic.p=NULL,word.p=NULL) {
  # K topics
  # M number of documents
  # D list of documents
  # W is total unique words in corpus
  M <- length(D)
  alldocs <- unlist(D)
  uniquewords <- unique(alldocs)
  W <- length(uniquewords)
  if(is.null(topic.p)) topic.p <- rep((1/K),K) # Prior on topics
  if(is.null(word.p)) word.p <- rep((1/W),W) # Prior on words
  
  phi <- matrix(NA,nrow=M,ncol=length(topic.p)) #M by K
  theta <- matrix(NA,nrow=length(word.p),ncol=K) # V by K
  
  
  for(i in 1:M) { # Document-Topic loop
    
    # For each document, generate a probability of topics in doc
    phi[i,] <- rdirichlet(1,topic.p)
  }
  
  Nj <- c() # Lengths of each documents
  
  for(j in 1:K) { # Topic-word loop
    # Sample out of a dirichlet for each document
    # This shows each document is a mixture of topics, initially
    theta[,j] <- rdirichlet(1,word.p)
    Nj[j] <- length(D[[j]])
  }
  
  z <- list()
  w <- list()
  for(ii in 1:M) { # For each document
    z[[ii]] <- rep(NA,times=Nj[ii])
    w[[ii]] <- z[[ii]]
    for(jj in 1:Nj[ii]) { # For each word in each document
      # Sample a random topic assignment, from topic proportions theta
      # Sample a random word out of that topic
      # The z_ii,jj indicates which topic was selected for what word
      z[[ii]][jj] <- sample(1:K,size=1,prob=phi[ii,])
      w[[ii]][jj] <- sample(1:Nj[ii],size=1,prob=theta[,z[[ii]][jj]])
    }
  }
  
  return (list(phi=phi, theta=theta, topics.z=z, words.w = w,topic.p=topic.p, word.p= word.p))
}


gen_lda_tf <- function(K,D,topic.p=NULL,word.p=NULL) {
  # K topics
  # M number of documents
  # D TF object of documents
  # W is total unique words in corpus
  M <- D$nrow
  #alldocs <- unlist(D)
  #uniquewords <- unique(alldocs)
  W <- D$ncol
  if(is.null(topic.p)) topic.p <- rep((1),K) # Prior on topics, doesn't need to sum to 1, this gives uniform prior
  if(is.null(word.p)) word.p <- rep((1/W),W) # Prior on words
  
  phi <- matrix(NA,nrow=M,ncol=length(topic.p)) #M by K
  theta <- matrix(NA,nrow=length(word.p),ncol=K) # V by K
  
  Nj <- c() # Lengths of each documents
  
  for(i in 1:M) { # Document-Topic loop
    
    # For each document, generate a probability of topics in doc
    phi[i,] <- rdirichlet(1,topic.p)
    Nj[i] <- length(D$v[D$i==i])  
  }
  
  
  for(j in 1:K) { # Topic-word loop
    # Sample out of a dirichlet for each document
    # This shows each document is a mixture of topics, initially
    theta[,j] <- rdirichlet(1,word.p)
  }
  
  #theta <- theta/rowSums(theta) # Renormalize
  # Don't renormalize because with such sparse words, the rdirichlet will give 0 probability to words
  #print(theta[,1])
  # w and z are V by M
  z <- matrix(NA,nrow=D$ncol,ncol=D$nrow)
  w <- z
  for(ii in 1:M) { # Choosing a topic
      wordnum <- D$j[D$i==ii] # This will find the word (k)
      z[wordnum,ii] <- sample(1:K,size=length(wordnum),prob=phi[ii,],replace=TRUE)
      #print(z[wordnum,ii])
      # This statement above is fine, sampling for all j's
      for(jj in 1:Nj[ii]) { # Wordnum is all the indices, out of all words
        w[wordnum[jj],ii] <- sample(1:W,size=1,prob=theta[,z[wordnum[jj],ii]])
        #print(w[wordnum[jj],ii])
      }
  }
  
  
  # for(ii in 1:M) { # For each document
  #   z[[ii]] <- rep(NA,times=Nj[ii])
  #   w[[ii]] <- z[[ii]]
  #   for(jj in 1:Nj[ii]) { # For each word in each document
  #     # Sample a random topic assignment, from topic proportions theta
  #     # Sample a random word out of that topic
  #     # The z_ii,jj indicates which topic was selected for what word
  #     z[[ii]][jj] <- sample(1:K,size=1,prob=phi[ii,])
  #     w[[ii]][jj] <- sample(1:Nj[ii],size=1,prob=theta[,z[[ii]][jj]])
  #   }
  # }
  
  return (list(phi=phi, theta=theta, topics.z=z, words.w = w,topic.p=topic.p, word.p= word.p,Ns=Nj))
}

library(topicmodels)
Assoc.gen <- gen_lda_tf(K=2,D=AssociatedPress)
# Test the generative model
# Do I enter i
# Should I convert to document term matrix?


# Will need to implement Collapsed Gibbs here
# Because I already know the distributions I don't need a metropolis rule
# I already have the direct distribution, I just need to sample out of it

collap.gibbs <- function(K, D, phi, theta, Z, W, priortopic,priorword,Ns,iters=1000) {
  # Need to sample for z and w
  # W doesn't really matter since it will be for a i j combo
  # z sample the probabilities from a multinomial
  # Need three for loops, one for all the topics, one for all the documents, and one for all the words
  W <- D$ncol
  #alldocs <- unlist(D)
  #uniquewords <- unique(alldocs)
  M <- D$nrow
  alpha <- priortopic
  beta <- priorword
  # Will have Z in here
  
  
  postz <- array(dim=c(K,M,W))
  # z matrix is M by W full with topics
  num.j.m <- postz
  num.j.v <- postz
  summed.j.m <- array(dim=c(K,M,W,K))
  summed.j.v <- summed.j.m
  # I need to keep track of the counts for each topic in words
  
  
  # These for loops are for obtaining the posterior distribution for each topic
    
  for(j in 1:K) { # For each topic
    for(m in 1:M) { # For each document to pick the word
      #wordnum <- D$i[c(D$j==m)] # This will find the docs (k)
      # All words in doc m
      wordnum <- D$j[D$i==m]
      for(n in 1:Ns[m]) { # This will be the nth 
        # sample the topic and word combination pertaining to topic j
        # We have identified that we are in doc m and word n
        #print(wordnum)
        v <- wordnum[n] # Current word
        #print(v)
        # Find all places in the matrix with the index of the topic
        ind.j.m <- which(Z==j,arr.ind=TRUE)  # Get indices of the topic
        #print(dim(ind.j.m))
        #print((ind.j.m[1:15,2]))
        #print(ind.j.m)
        if(is.null(dim(ind.j.m))) ind.j.m <- ind.j.m[which(ind.j.m[2]==m)] # Vectorized notation
        else ind.j.m <- ind.j.m[which(ind.j.m[,2]==m),]
        #print(dim(ind.j.m))
        #print(ind.j.m)
        # Now I need to remove the nth word
        if(is.null(dim(ind.j.m))) {
          #print(ind.j.m)
          ind.j.m <- ind.j.m[which(ind.j.m[1]!=v)]
        }
        else ind.j.m <- ind.j.m[which(ind.j.m[,1]!=v),]
        #print(dim(ind.j.m))
        #print(ind.j.m)
        # Now get the vocab for these indices
        
        #if(is.null(dim(ind.j.m))) ind.j.m <- as.matrix(ind.j.m,nrow=1,ncol=2)
        
        if(is.null(dim(ind.j.m))) num.j.m[j,m,wordnum[n]] <- sum(D$v[ind.j.m[1]],na.rm=T)
        else num.j.m[j,m,wordnum[n]] <- sum(D$v[ind.j.m[,1]],na.rm=T)
        
        #print(ind.j.m[,1])
        #num.j.m <- sum(D$v[ind.j.m[,1]])/
        
        #print(num.j.m)
        
        #print(paste("Hello",num.j.m))
        
        ind.j.v <- which(Z==j,arr.ind = T)
        
        if(is.null(dim(ind.j.v))) ind.j.v <- ind.j.v[which(ind.j.v[1]==v)]
        else ind.j.v <- ind.j.v[which(ind.j.v[,1]==v),]
        
        if(is.null(dim(ind.j.v))) ind.j.v <- ind.j.v[which(ind.j.v[2]!=m)]
        else ind.j.v <- ind.j.v[which(ind.j.v[,2]!=m),]
        
        if(is.null(dim(ind.j.v))) num.j.v[j,m,wordnum[n]] <- sum(D$v[ind.j.v[1]],na.rm=T)
        else num.j.v[j,m,wordnum[n]] <- sum(D$v[ind.j.v[,1]],na.rm=T)
        
        # Just checking for vectorized
        
        # Now I have the number of words in doc m with word n
        # And number of words in all docs, not including selecting doc
        
        inter.summed.j.m <- c()
        inter.summed.j.v <- c()
        for(jj in 1:K) {
          # s.j.m is sum over j and m
          s.j.m <- which(Z==jj,arr.ind=TRUE)
          #print(dim(s.j.m))
          if(is.null(dim(s.j.m))) s.j.m <- s.j.m[which(s.j.m[2]==m)]
          else s.j.m <- s.j.m[which(s.j.m[,2]==m),]
          #print(dim(s.j.m))
          
          if(is.null(dim(s.j.m))) {
            #print(s.j.m)
            s.j.m <- s.j.m[which(s.j.m[1]!=v)]
          }
          else s.j.m <- s.j.m[which(s.j.m[,1]!=v),]
          #print("SJM")
          #print(dim(s.j.m))
          #print(s.j.m)
          if(is.null(dim(s.j.m))) inter.summed.j.m[jj] <- sum(D$v[s.j.m[1]],na.rm=T) + 0
          else inter.summed.j.m[jj] <- sum(D$v[s.j.m[,1]],na.rm=T) + 0
          
          
          s.j.v <- which(Z==jj,arr.ind = T)
          
          if(is.null(dim(s.j.v))) s.j.v <- s.j.v[which(s.j.v[1]==v)] 
          else s.j.v <- s.j.v[which(s.j.v[,1]==v),]
          
          if(is.null(dim(s.j.v))) s.j.v <- s.j.v[which(s.j.v[2]!=m)]
          else s.j.v <- s.j.v[which(s.j.v[,2]!=m),]
          #print("SJV")
          #print(dim(s.j.v))
          #print(head(s.j.v))
          if(is.null(dim(s.j.v))) inter.summed.j.v[jj] <- sum(D$v[s.j.v[2]],na.rm=T) + 0
          else inter.summed.j.v[jj] <- sum(D$v[s.j.v[,2]],na.rm=T) + 0
          #print(summed.j.m[jj])
          #print(summed.j.v[jj])
        }
        ######################################################
        # Assignment section to record all numbers
        ######################################################
        #print("Hi")
        #print(inter.summed.j.m)
        #print(inter.summed.j.v)
        #readline(prompt="Press [enter] to continue")
        
        summed.j.m[j,m,wordnum[n],] <- inter.summed.j.m
        summed.j.v[j,m,wordnum[n],] <- inter.summed.j.v
        
        
        
        topicprob <-(num.j.m[j,m,wordnum[n]] + alpha[j])/(sum(summed.j.m[j,m,wordnum[n],] + alpha,na.rm=T))
        #print(topicprob)
        topicword <-(num.j.v[j,m,wordnum[n]] + beta[n])/(sum(summed.j.v[j,m,,j] + beta,na.rm=T))
        #print(topicword)
        # This collapsed isn't sampling, it's just calculating
        #print(n)
        postz[j,m,wordnum[n]] <- topicprob*topicword
        #print(postz[j,m,wordnum[n]])
        #readline(prompt="Press [enter] to continue")
        
        #print("Probs")
        #print(postz[j,m,wordnum[n]])
      }
      #(postz[j,2,])
      print(m)
    }
    # Renormalize per topic
    #postz[j,,] <- postz[j,,]/sum(postz[j,,],na.rm=T) # Need to remove NAs because of the Nj[m], words in each doc
    # We don't really need to normalize since sample will take any set of probabilities
    # Num rows will be the max of num words in a doc
  }
  
  # This is to get the posterior with integrated out alpha beta
  
  
  ###########################################3
  # Now that I have the posterior distribution
  # Get collapsed Gibbs samples out of this
  # Collapsed Gibbs steps
  ############################################
  
  # Need to create an array for each time step of the Collapsed Gibbs
  
  
  print("Gibbs part")
  cGibbsPost <- array(dim=c(iters,K,M,W))
  #postz <- postz[which(is.na(postz))] <- eps
  cGibbsPost[1,,,] <- postz
  #print(postz[,2,])
  
  for(t in 2:iters) {
    cGibbsPost[t,,,] <- cGibbsPost[(t-1),,,] # So we're updating the current Gibbs post
    print(t)
    for(m in 1:M) {
      wordnum <- D$j[D$i==m]
      #print(wordnum)
      for(n in 1:Ns[m]) {
        # We are updating the counts here for whichever k we choose
        
        ###################################
        # Sampling a topic for m, n combo
        ###################################
        #print(cGibbsPost[t,,m,wordnum[n]])
        #print(postz[,m,wordnum[n]])
        #print(cGibbsPost[t,,m,wordnum[n]])
        samp <- sample(1:K, size=1, prob = cGibbsPost[t,,m,wordnum[n]])
        # Samp is a value from 1 to K
        # Once we get this sample from 1:K, then add the word combo and recalculate the posterior for that topic
        # Might need to subtract here, but not sure how?
        # Update the numeric counts
        # How many words to add?
        wordadd <- wordnum[n]
        wordfreq <- D$v[wordadd] # frequency to add to the total of it
        
        ind.j.v <- which(Z==samp,arr.ind = T)
        
        # I think this is right, but I'm not sure
        
        num.j.m[samp,m,wordnum[n]] <- num.j.m[samp,m,wordnum[n]] + wordfreq
        num.j.v[samp,m,wordnum[n]] <- num.j.v[samp,m,wordnum[n]] + wordfreq
        
        summed.j.m[samp,m,wordnum[n],samp] <- summed.j.m[samp,m,wordnum[n],samp] + wordfreq
        summed.j.v[samp,m,wordnum[n],samp] <- summed.j.v[samp,m,wordnum[n],samp] + wordfreq
        
        # How do I update these counts?
        
        # Need to check for if equality, if so subtract counts
        # Will need a better accounting of the Gibbs sampler, keep track of a matrix
        
        topicprob <-(num.j.m[samp,m,wordnum[n]] + alpha[j])/(sum(summed.j.m[samp,m,wordnum[n],] + alpha, na.rm=T))
        #print(topicprob)
        topicword <-(num.j.v[samp,m,wordnum[n]] + beta[n])/(sum(summed.j.v[samp,m,,samp] + beta, na.rm = T))
        #print(topicword)
        # This collapsed isn't sampling, it's just calculating
        cGibbsPost[t,samp,m,wordnum[n]] <- topicprob*topicword
        #print(cGibbsPost[t,samp,m,wordnum[n]])
        
        # Need to update after each n is chosen, update posterior probability
      }
    }
  }
  #cGibbsPost <- apply()
  # Well this takes 5 hours apparently
  # Think some way later about cacheing the n's
  # Later try to make this faster
  # And then I need to sample out of this? 
  return (list(post.K.M.W=postz,Gibbspost=cGibbsPost))
  # This z is multinomial, maybe try to get the theta and phi from this as well
  # Will still need estimates of each theta and phi
  # Actually these theta and phi are found in the topicprob and topicword for each
  # Need to get better caching
  
}

# So I'm doing Variational EM and collapsed Gibbs
#vEM

# Read in the file, has AllContent

library(dplyr)

newdf <- complete.cases(AllContent2019_03_02$content)
interdf <- AllContent2019_03_02[newdf,]
#interdf <- distinct(interdf)

# Now do the LDA on the short content

mytf <- get.TF(interdf$content[1:10])[[2]]
content.gen <- gen_lda_tf(K=2,D=mytf)  
  
system.time(cGibbsinfer.content <- collap.gibbs(K=2,D=mytf,
                phi = content.gen$phi, theta=content.gen$theta, Z=content.gen$topics.z, W=content.gen$words.w, 
                priortopic=content.gen$topic.p, priorword = content.gen$word.p, Ns = content.gen$Ns))
Assoc.gen <- gen_lda_tf(K=2,D=AssociatedPress)

anothertest <- collap.gibbs(K=2,D=AssociatedPress,phi = Assoc.gen$phi, theta= Assoc.gen$theta, Z=Assoc.gen$topics.z, W=Assoc.gen$words.w,
                            priortopic = Assoc.gen$topic.p, priorword = Assoc.gen$word.p, Ns= Assoc.gen$Ns,iters=5)
# The algorithm that I saw from other shows that all words have some probability of showing up

cont400 <- get.TF(interdf$content[1:400])[[2]]
cont400.gen <- gen_lda_tf(K=2,D=cont400)

system.time(anothertest2 <- collap.gibbs(K=2,D=cont400,phi = cont400.gen$phi, theta= cont400.gen$theta, Z=cont400.gen$topics.z, W=cont400.gen$words.w,
                            priortopic = cont400.gen$topic.p, priorword = cont400.gen$word.p, Ns= cont400.gen$Ns,iters=100))

cont400.gen4 <- gen_lda_tf(K=4,D=cont400)
system.time(anothertest3 <- collap.gibbs(K=4,D=cont400,phi = cont400.gen4$phi, theta= cont400.gen4$theta, Z=cont400.gen4$topics.z, W=cont400.gen4$words.w,
                            priortopic = cont400.gen4$topic.p, priorword = cont400.gen4$word.p, Ns= cont400.gen4$Ns,iters=100))

# Figure out how to compare these? Maybe by the perplexity or something like that
loglikeli <- function(gam,phi,alpha,beta,wordindices) {
  # gam is a vector from i=1:K, where K is number of topics
  # phi is a N_d by K matrix (N_d is number of words in current doc) there will be some NAs for the W!=N_d
  # alpha is an updated parameter on prior on topics, K dimensional
  # beta is a prior on the words within a topic, K by V (W) matrix
  # words is the wordvector using here, length W
  # wordindices are the indices for the phi needed over all the i's (words in doc D)
  #print(length(gam))
  #print(dim(phi))
  #print(dim(phi[wordindices,]))
  #print(length(alpha))
  #print(dim(beta))
  parts <- c()
  parts[1] <- log(gamma(sum(alpha))) - sum(log(gamma(alpha))) + sum( (alpha-1)*(digamma(gam)-digamma(sum(gam))) )
  print(parts[1])
  parts[2] <- sum(phi[wordindices,]%*%(digamma(gam)-digamma(sum(gam))),na.rm=T)
  print(parts[2])
  word <- rep(0,ncol(beta))
  word[wordindices] <- 1
  #print("hi")
  #print(dim(t(word%*%t(log(beta)))))
  #print(t(word%*%t(log(beta))))
  #print(word)
  #print(word%*%t(log(beta)))
  #print(beta)
  #print(rowSums(beta))
  #print(log(beta))
  logbeta <- log(beta)
  print(logbeta)
  logbeta[is.infinite(logbeta)] <- -600
  print("i")
  print(logbeta)
  #print(logbeta)
  #print(phi[wordindices,])
  #print(word%*%t(log(beta)))
  #print(word%*%t(logbeta))
  #print(phi[wordindices,])
  #print(log(phi[wordindices,]))
  logphi <- log(phi[wordindices,])
  logphi[is.infinite(logphi)] <- log(eps)
  parts[3] <- sum(phi[wordindices,]%*%t(word%*%t(logbeta)))
  print(parts[3])
  parts[4] <- -log(gamma(sum(gam))) + sum(log(gamma(gam))) - sum( (gam-1)* (digamma(gam) - digamma(sum(gam))) )
  print(parts[4])
  parts[5] <- -sum(phi[wordindices,]%*%t(logphi))
  print(parts[5])
  
  # Parts 3 and 5 are failing
  
  return (sum(parts))
}


vEM <- function(genldaobj,D,TT=100,iters=500,tol=0.01) {
  #############################################################################################
  # This function will use the variational EM algorithm as discussed by Blei et al. 2003
  ############################################################################################
  
  
  
  M <- nrow(genldaobj$phi) #docs
  K <- ncol(genldaobj$phi) #topics
  W <- D$ncol
  #####################
  saveSkv <- array(dim=c(K,W,iters))
  saveSa <- matrix(NA,nrow=K,ncol=iters) # Not even sure if I need these anymore
  #####################
  alphasave <- matrix(NA,nrow=K,ncol=iters)
  alphasave[,1] <- genldaobj$topic.p
  stop <- FALSE
  betasave <- array(dim=c(K,W,iters))
  betasave[,,1] <- t(genldaobj$theta) # Because theta is W by K I need it to be K by W
  loglikeli.save <- matrix(NA,nrow=TT,ncol=iters)
  save.gam.d.k <- array(dim=c(M,K,iters))
  for(vv in 1:iters) {
    gam.d.k <- array(dim=c(M,K,TT))
    #gam.d.k <- matrix(NA,nrow=M,ncol=K)
    #phi.d.i.k <- array(dim=c(M,W,K))
    phi.d.i.k <- array(dim=c(M,W,K,TT))
    #Skv <- matrix(NA,nrow=K,ncol=W)
    Sdik <- phi.d.i.k
    wordnums <- list()
    for(d in 1:M) {
      wordnums[[vv]] <- D$i[D$j==d] # Gives correct word indices
      for(k in 1:K) {
        gam.d.k[d,k,1] <- genldaobj$Ns[d]/K # Initializations
        #print(gam.d.k)
      }
      #Sa <- matrix(NA,nrow=M,ncol=TT)
      for(t in 2:TT) { 
        stop <- FALSE #E step convergence iterations
        for(i in 1:genldaobj$Ns[d]) {
          for(k2 in 1:K) {
            phi.d.i.k[d,wordnums[[vv]][i],k2,t] <- betasave[k2,wordnums[[vv]][i],vv]*exp(digamma(gam.d.k[d,k2,(t-1)])) #genldaobj$theta[wordnums[i],k2]*exp(digamma(gam.d.k[d,k2,(t-1)]))
          }
          # Now normalize phi
          phi.d.i.k[d,wordnums[[vv]][i],,t] <- phi.d.i.k[d,wordnums[[vv]][i],,t]/sum(phi.d.i.k[d,wordnums[[vv]][i],,t],na.rm=T)
        }
        for(k3 in 1:K) {
          gam.d.k[d,k3,t] <- alphasave[k3,vv] + sum(phi.d.i.k[d,,k3,t],na.rm=T)
          #print(gam.d.k[d,k3,t])
        }
        #print(phi.d.i.k[d,,1,t])
        #print(sum(phi.d.i.k[d,,1,t],na.rm=T))
        #print(gam.d.k[d,1,t])
      # if(t > 2) {
      #   print(KLdiverge(phi.d.i.k[d,wordnums,,(t-1),phi.d.i.k[d,wordnums,,t]]))
      #   if(KLdiverge(phi.d.i.k[d,wordnums,,(t-1)],phi.d.i.k[d,wordnums,,(t)]) < 0.1 && KLdiverge(gam.d.k[d,,(t-1)],gam.d.k[d,,(t)]) < 0.1) 
      #     break # Look for KL convergence on other gam as well
      # }
        loglikeli.save[t,vv] <- loglikeli(gam.d.k[d,,t],phi.d.i.k[d,,,t],alphasave[,vv],betasave[,,vv],wordindices = wordnums[[vv]])
        print(loglikeli.save[t,vv])
        if(t > 2) { # Just make sure there are multiple time points
          if(abs(loglikeli.save[t,vv]-loglikeli.save[(t-1),vv]) < tol) {
            stop <- TRUE
            break
          }
        }
        
        if(stop==T){break}
        
      # Need to think about the convergence criterion  
        
      # for suff stats
      
      # Maybe I don't need a for loop for the Ns
      # Look at words in doc d
      }
      indswords <- D$i[D$j==d]
      indvec <- rep(0,W)
      wghts <- D$v[indswords]
      indvec[indswords] <- wghts
      #print(wghts)
      #print(phi.d.i.k[d,,,t])
      Sdik[d,,,t]<-indvec*phi.d.i.k[d,,,t]
      save.gam.d.k[d,,vv] <- gam.d.k[d,,t]
      # 
      # for(ii in 1:genldaobj$Ns[d]) {
      #   #Skv[,wordnums[ii]] <- which(genldaobj$words.w[,ii]==wordnums[ii])*phi.d.i.k[d,wordnums[ii],]
      #   # Instead of multiplying, I'll just using the parts of phi where the indices are equal
      #   # Convert to indicator
      #   inds <- rep(0,W)
      #   # Want to find alll words in this dth
      #   print(wordnums[ii])
      #   inds[wordnums[ii]] <- 1
      #   Skv[d,wordnums[ii]] <- inds*phi.d.i.k[d,wordnums[ii],]
      # }
    
      #Sa[d] <- sum(digamma(gam.d.k) - K*digamma(sum(gam.d.k[d,])))
      #print(dim(gam.d.k))
      #print(gam.d.k[,,t])
      #print(dim(digamma(gam.d.k[,,t])))
      #print(dim(digamma(gam.d.k[d,,t])))
      #print(sum(digamma(gam.d.k[,,t])-K*digamma(gam.d.k[d,,t]),na.rm=T))
      #print(dim(Sa))
      #Sa[d,t] <- sum(digamma(gam.d.k[,,t])-K*digamma(gam.d.k[d,,t]),na.rm=T)
      #print(Sa[d,t])
      #print(t)
      #print(Sa)
      }
    #print(d)  
    Skv <- apply(Sdik[,,,t],c(2:3),sum,na.rm=T) #To be W by K, maybe this'll still work
    Skv <- t(Skv) # To be K and W (v)
    # Might need to renormalize
    #Skv <- apply(Skv,2,normalize)
  
    
    #SaSave <- sum(Sa[,t]) # At the best element, the latest time step for each doc
    
    # Now this is the m step
    saveSkv[,,vv] <- Skv
    # saveSa[vv] <- SaSave
    # # Need to change the theta and topic.p to beta and alpha respectively
    # print(vv)
    # # Now to implement the M step
    # 
    
    #########################
    # M step, update beta
    #########################
    
    zz <- M*trigamma(sum(alphasave[,vv]))
    gk <- M*digamma(sum(alphasave[,vv])) - M*digamma(alphasave[,vv]) + colSums(digamma(save.gam.d.k[,,vv]) - digamma(rowSums(save.gam.d.k[,,vv])))
    # Colsums over K, Rowsums over D
    hkk <- M*(trigamma(alphasave[,vv])-trigamma(sum(alphasave[,vv])))
    const <- (sum(gk/hkk)) / (1/zz + sum(1/hkk))
    print("this")
    print(zz)
    print(gk)
    print(hkk)
    print(const)
    if(vv!=iters) {
      # I guess at the most recent iteration of it
      #print("kkkkkk")
      #print(saveSkv[,,vv])
      #print("T1")
      #print(apply(saveSkv[,,vv],1,normalize))
      #print(colSums(apply(saveSkv[,,vv],1,normalize)))
      #print("T2")
      #print(apply(saveSkv[,,vv],2,normalize))
      #print(rowSums(apply(saveSkv[,,vv],2,normalize)))
      betasave[,,(vv+1)] <- t(apply(saveSkv[,,vv],1,normalize)) 
      alphasave[,(vv+1)] <- alphasave[,vv] - ((gk-const)/hkk)
      print(betasave[,,(vv+1)])
      print(alphasave[,(vv+1)])
      print("Dpme")
      
      # Ask about these zeros in the beta matrix and look up more tomorrow.
      
    }
    
  }
    # if(vv!=iters){
    #   betasave[,,(vv+1)] <- apply(saveSkv[,,vv],2,normalize)
    #   print(betasave[,,vv])
    #   alphasave[(vv+1)] <- alphasave[vv]
    #   }
    # }
  #return (list(Skv=Skv,Sa=Sa,gam=gam.d.k,phi=phi.d.i.k))
  
  #return (list(Skv=saveSkv,Sa=saveSa,beta=betasave,alpha=alphasave))
  
  return (list(alpha=alphasave,beta=betasave,loglikelihood=loglikeli.save,phi=saveSkv,gam=save.gam.d.k))
    
  # Might have to change around some variables to be able to model them
  # Should make topic p one dimensional by t dimensional until convergence, then get a cached alpha
  # Get a beta matrix then have it updated by time for the probabilities of the words
}



# Look up newton raphson tomrrow for updates on alpha
# Look up updates on collapsed gibbs and see why it's running so slow
# See what I need to see in the variational EM, maybe look at Bleis stuff



tryloglikeli <- function(X) {
  # gam is a vector from i=1:K, where K is number of topics
  # phi is a N_d by K matrix (N_d is number of words in current doc) there will be some NAs for the W!=N_d
  # alpha is an updated parameter on prior on topics, K dimensional
  # beta is a prior on the words within a topic, K by V (W) matrix
  # words is the wordvector using here, length W
  # wordindices are the indices for the phi needed over all the i's (words in doc D)
  
  # X is 2*Nd + 2 by K
  # xvector
  X <- matrix(as.vector(t(X)),ncol=2,byrow=T)
  
  dims <- dim(X)
  Nd <- (dim(X)[1]-2)/2
  
  
  gam <- X[1,]
  alpha <- X[2,]
  phi <- X[3:(2+Nd),]
  beta <- X[(3+Nd):nrow(X),]
  beta <- t(beta)
  
  
  
  parts <- c()
  parts[1] <- log(gamma(sum(alpha))) - sum(log(gamma(alpha))) + sum( (alpha-1)*(digamma(gam)-digamma(sum(gam))) )
  print(parts[1])
  parts[2] <- sum(phi%*%(digamma(gam)-digamma(sum(gam))),na.rm=T)
  print(parts[2])
  word <- rep(1,ncol(beta))
  logbeta <- log(beta)
  print(logbeta)
  logbeta[is.infinite(logbeta)] <- -600
  print("i")
  print(logbeta)

  
  logphi <- log(phi)
  logphi[is.infinite(logphi)] <- log(eps)
  parts[3] <- sum(phi%*%t(word%*%t(logbeta)))
  print(parts[3])
  parts[4] <- -log(gamma(sum(gam))) + sum(log(gamma(gam))) - sum( (gam-1)* (digamma(gam) - digamma(sum(gam))) )
  print(parts[4])
  parts[5] <- -sum(phi%*%t(logphi))
  print(parts[5])
  
  # Parts 3 and 5 are failing
  
  return (sum(parts))
}

testgam <- c(23,23)/2
testalpha <- content.gen$topic.p
testbeta <- content.gen$theta[mytf$i[mytf$j==1],]
testbeta <- testbeta/rowSums(testbeta)
testphi <- matrix(NA,nrow=23,ncol=2)
for(k2 in 1:2) {
    testphi[,k2] <- testbeta[,k2]*exp(digamma(testgam[k2])) #genldaobj$theta[wordnums[i],k2]*exp(digamma(gam.d.k[d,k2,(t-1)]))
}
newX <- rbind(matrix(testgam,nrow=1),matrix(testalpha,nrow=1),testbeta,testphi)

input <- as.vector(t(newX))
#matrix(as.vector(t(newX)),ncol=2,byrow=T)
mytest1EM <- optim(input,tryloglikeli)

#############################################################
# This is the section for making plots for the report
#############################################################
wordd <- 3207

par(mfcol=c(1,2)) # 21
plot(apply(anothertest2$Gibbspost[,1,,wordd],1,mean,na.rm=T),type="l",ylim=c(0,0.25),main="Trace plot of average probability of word \"points\" in each topic",ylab="Probability")
lines(apply(anothertest2$Gibbspost[,2,,wordd],1,mean,na.rm=T),col=2)

plot(apply(anothertest3$Gibbspost[,1,,wordd],1,mean,na.rm=T),type="l",ylim=c(0,0.25),main="Trace plot of probability of word \"points\" in each topic",ylab="Probability")
lines(apply(anothertest3$Gibbspost[,2,,wordd],1,mean,na.rm=T),col=2)
lines(apply(anothertest3$Gibbspost[,3,,wordd],1,mean,na.rm=T),col=3)
lines(apply(anothertest3$Gibbspost[,4,,wordd],1,mean,na.rm=T),col=4)

wordd <- 3207

wordd <- which(cont400$dimnames$Terms=="points")

par(mfcol=c(1,2)) # 21
plot(apply(anothertest2$Gibbspost[,1,,wordd],1,mean,na.rm=T),type="l",ylim=c(0,0.25),main="Average probability of word \"points\" in each topic",ylab="Probability")
lines(apply(anothertest2$Gibbspost[,2,,wordd],1,mean,na.rm=T),col=2)
legend("topright",c("Topic 1","Topic 2"),col=c(1,2),lwd=1)


plot(apply(anothertest3$Gibbspost[,1,,wordd],1,mean,na.rm=T),type="l",ylim=c(0,0.25),main="Average probability of word \"points\" in each topic",ylab="Probability")
lines(apply(anothertest3$Gibbspost[,2,,wordd],1,mean,na.rm=T),col=2)
lines(apply(anothertest3$Gibbspost[,3,,wordd],1,mean,na.rm=T),col=3)
lines(apply(anothertest3$Gibbspost[,4,,wordd],1,mean,na.rm=T),col=4)
legend("topright",c("Topic 1","Topic 2","Topic 3","Topic 4"),col=seq(1,4,by=1),lwd=1)


wordd <- which(cont400$dimnames$Terms=="chars")

par(mfcol=c(1,2)) # 21
plot(apply(anothertest2$Gibbspost[,1,,wordd],1,mean,na.rm=T),type="l",ylim=c(0,0.25),main="Average probability of word \"chars\" in each topic",ylab="Probability")
lines(apply(anothertest2$Gibbspost[,2,,wordd],1,mean,na.rm=T),col=2)
legend("topright",c("Topic 1","Topic 2"),col=c(1,2),lwd=1)


plot(apply(anothertest3$Gibbspost[,1,,wordd],1,mean,na.rm=T),type="l",ylim=c(0,0.25),main="Average probability of word \"chars\" in each topic",ylab="Probability")
lines(apply(anothertest3$Gibbspost[,2,,wordd],1,mean,na.rm=T),col=2)
lines(apply(anothertest3$Gibbspost[,3,,wordd],1,mean,na.rm=T),col=3)
lines(apply(anothertest3$Gibbspost[,4,,wordd],1,mean,na.rm=T),col=4)
legend("topright",c("Topic 1","Topic 2","Topic 3","Topic 4"),col=seq(1,4,by=1),lwd=1)


wordd <- which(cont400$dimnames$Terms=="points")

par(mfcol=c(1,2)) # 21
plot(apply(anothertest2$Gibbspost[,1,,wordd],1,mean,na.rm=T),type="l",ylim=c(0,0.15),main=paste0("Average probability of word \"",cont400$dimnames$Terms[wordd],"\" in each topic"),ylab="Probability")
lines(apply(anothertest2$Gibbspost[,2,,wordd],1,mean,na.rm=T),col=2)
legend("topright",c("Topic 1","Topic 2"),col=c(1,2),lwd=1)


plot(apply(anothertest3$Gibbspost[,1,,wordd],1,mean,na.rm=T),type="l",ylim=c(0,0.15),main=paste0("Average probability of word \"",cont400$dimnames$Terms[wordd],"\" in each topic"),ylab="Probability")
lines(apply(anothertest3$Gibbspost[,2,,wordd],1,mean,na.rm=T),col=2)
lines(apply(anothertest3$Gibbspost[,3,,wordd],1,mean,na.rm=T),col=3)
lines(apply(anothertest3$Gibbspost[,4,,wordd],1,mean,na.rm=T),col=4)
legend("topright",c("Topic 1","Topic 2","Topic 3","Topic 4"),col=seq(1,4,by=1),lwd=1)

wordd <- which(cont400$dimnames$Terms=="trump")
#wordd <- 2500

par(mfcol=c(1,2)) # 21
plot(apply(anothertest2$Gibbspost[,1,,wordd],1,mean,na.rm=T),type="l",ylim=c(0,0.15),main=paste0("Average probability of word \"",cont400$dimnames$Terms[wordd],"\" in each topic"),ylab="Probability")
lines(apply(anothertest2$Gibbspost[,2,,wordd],1,mean,na.rm=T),col=2)
legend("topright",c("Topic 1","Topic 2"),col=c(1,2),lwd=1)

plot(apply(anothertest3$Gibbspost[,1,,wordd],1,mean,na.rm=T),type="l",ylim=c(0,0.15),main=paste0("Average probability of word \"",cont400$dimnames$Terms[wordd],"\" in each topic"),ylab="Probability")
lines(apply(anothertest3$Gibbspost[,2,,wordd],1,mean,na.rm=T),col=2)
lines(apply(anothertest3$Gibbspost[,3,,wordd],1,mean,na.rm=T),col=3)
lines(apply(anothertest3$Gibbspost[,4,,wordd],1,mean,na.rm=T),col=4)
legend("topright",c("Topic 1","Topic 2","Topic 3","Topic 4"),col=seq(1,4,by=1),lwd=1)




plot(anothertest2$Gibbspost[,1,21,2355],type="l",ylim=c(0,0.25),main="Trace plot of average probability of word \"points\" in each topic",ylab="Probability")
lines(anothertest2$Gibbspost[,2,21,2355],col=2)

plot(anothertest3$Gibbspost[,1,21,2355],type="l",ylim=c(0,0.25),main="Trace plot of probability of word \"points\" in each topic",ylab="Probability")
lines(anothertest3$Gibbspost[,2,21,2355],col=2)
lines(anothertest3$Gibbspost[,3,21,2355],col=3)
lines(anothertest3$Gibbspost[,4,21,2355],col=4)






