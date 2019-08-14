suppressWarnings(suppressMessages(library("topicmodels")))
suppressWarnings(suppressMessages(library("tm")))
suppressWarnings(suppressMessages(library("SnowballC")))
suppressWarnings(suppressMessages(library("wordcloud")))
suppressWarnings(suppressMessages(library("RColorBrewer")))
suppressWarnings(suppressMessages(library("RCurl")))
suppressWarnings(suppressMessages(library("XML")))
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("forcats")))
# For gini
suppressWarnings(suppressMessages(library("ineq")))
# For ldagraphic
suppressWarnings(suppressMessages(library("ggpubr")))
suppressWarnings(suppressMessages(library("knitr")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("plyr")))

# Website is http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need

#++++++++++++++++++++++++++++++++++
# rquery.wordcloud() : Word cloud generator
# - http://www.sthda.com
#+++++++++++++++++++++++++++++++++++
# x : character string (plain text, web url, txt file path)
# type : specify whether x is a plain text, a web page url or a file path
# lang : the language of the text
# excludeWords : a vector of words to exclude from the text
# textStemming : reduces words to their root form
# colorPalette : the name of color palette taken from RColorBrewer package, 
# or a color name, or a color code
# min.freq : words with frequency below min.freq will not be plotted
# max.words : Maximum number of words to be plotted. least frequent terms dropped
# value returned by the function : a list(tdm, freqTable)
rquery.wordcloud <- function(x, type=c("text", "url", "file"), 
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
  
  ###### I'm going to remove quotes
  # Somehow
  
  # Text stemming
  if(textStemming) docs <- tm_map(docs, stemDocument)
  # Create term-document matrix
  tdm <- TermDocumentMatrix(docs)
  m <- as.matrix(tdm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  # check the color palette name 
  if(!colorPalette %in% rownames(brewer.pal.info)) colors = colorPalette
  else colors = brewer.pal(8, colorPalette) 
  # Plot the word cloud
  set.seed(1234)
  wordcloud(d$word,d$freq, min.freq=min.freq, max.words=max.words,
            random.order=FALSE, rot.per=0.35, 
            use.r.layout=FALSE, colors=colors)
  
  invisible(list(tdm=tdm, freqTable = d))
}
#++++++++++++++++++++++
# Helper function
#++++++++++++++++++++++
# Download and parse webpage
html_to_text<-function(url){
  library(RCurl)
  library(XML)
  # download html
  html.doc <- getURL(url)  
  #convert to plain text
  doc = htmlParse(html.doc, asText=TRUE)
  # "//text()" returns all text outside of HTML tags.
  # We also do not want text such as style and script codes
  text <- xpathSApply(doc, "//text()[not(ancestor::script)][not(ancestor::style)][not(ancestor::noscript)][not(ancestor::form)]", xmlValue)
  # Format text vector into one character string
  return(paste(text, collapse = " "))
}

get.TF <- function(x, type=c("text", "url", "file"), 
                   lang="english", excludeWords=NULL, 
                   textStemming=FALSE,  colorPalette="Dark2",
                   min.freq=3, max.words=200)
{ 
  library("tm")
  library("SnowballC")
  library("wordcloud")
  library("RColorBrewer") 
  suppressWarnings(library("topicmodels"))
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
  return (tdm)
}

RenormGini <- function(probvec) {
  # This takes in a probability vector, and returns a gini of greater than uniform weighing
  gtindex <- probvec >= 1/length(probvec)
  renormvec <- probvec[gtindex]/sum(probvec[gtindex])
  ginivec <- ineq(renormvec,type="Gini")
  return (ginivec)
}

KLdiverge <- function(d1,d2,dir="Forward") {
  # This function takes two discrete distribution vectors
  # Returns KL divergence between the two vectors
  # Check if d1,d2 are dists, if not, normalize
  d1 <- unlist(d1)
  d2 <- unlist(d2)
  if(sum(d1)!=1) d1 <- d1/sum(d1)
  if(sum(d2)!=1) d2 <- d2/sum(d2)
  Dkl <- rep(NA,times=length(d1))
  if(dir=="Forward") {
    # Check if d1 > d2 at certain points, then do -log of one over the other
    ind1 <- which(d1 < d2)
    Dkl[ind1] <- d1[ind1]*(log(d2[ind1])-log(d1[ind1]))
    ind2 <- which(d1 >= d2)
    Dkl[ind2] <- -d1[ind2]*(log(d1[ind2])-log(d2[ind2]))
    Dkl <- -sum(Dkl)
    #Dkl <- -sum(d1*log(d2/d1))
  }
  else {
    ind1 <- which(d2 < d1)
    Dkl[ind1] <- d2[ind1]*log(d1[ind1]/d2[ind1])
    ind2 <- which(d2 >= d1)
    Dkl[ind2] <- -d2[ind2]*log(d2[ind2]/d1[ind2])
    Dkl <- -sum(Dkl)
    #Dkl <- -sum(d2*log(d1/d2))
  }
  return (Dkl)
} 

JSdiverge <- function(d1,d2, dir="Forward") {
  M <- 0.5*(d1+d2)
  M <- unlist(M)
  JSD <- 0.5*(KLdiverge(d1,M,dir=dir)+KLdiverge(d2,M,dir=dir))
  return (JSD)
} 

ShannonIC <- function(d,gtunif=FALSE) {
  # Takes in a vector, returns a log sum val
  if(sum(d)!=1) d <- d/sum(d)
  if(gtunif==TRUE) {
    index <- d >= 1/length(d)
    d <- d[index]/sum(d[index])
  }
  return (-sum(d*log(d)))
}
# Could run some optim wrapper around this to get the optimial size of the vector

LDAgraphic <- function(textvec,numtopics,seedset=150,stopwordlist=c("said","will"),topn=10) {
  # home -> #setwd("~/ATD Group Spring 2019/Analysis for Feb 13")
  #setwd("~/")
  suppressWarnings(suppressMessages(library("topicmodels")))
  suppressWarnings(suppressMessages(library("tm")))
  suppressWarnings(suppressMessages(library("SnowballC")))
  suppressWarnings(suppressMessages(library("wordcloud")))
  suppressWarnings(suppressMessages(library("RColorBrewer")))
  suppressWarnings(suppressMessages(library("RCurl")))
  suppressWarnings(suppressMessages(library("XML")))
  suppressWarnings(suppressMessages(library("ggplot2")))
  suppressWarnings(suppressMessages(library("dplyr")))
  suppressWarnings(suppressMessages(library("forcats")))
  # For gini
  suppressWarnings(suppressMessages(library("ineq")))
  # For ldagraphic
  suppressWarnings(suppressMessages(library("ggpubr")))
  suppressWarnings(suppressMessages(library("knitr")))
  suppressWarnings(suppressMessages(library("dplyr")))
  suppressWarnings(suppressMessages(library("plyr")))
  textlist <- list()
  docs <- list()
  iterator <- 1
  #stopwordlist <- c("said", "will")
  # Will need to gsub those words out for blanks
  suppressWarnings(library("topicmodels"))
  for(iter in 1:length(textvec)) {
    textlist[[iter]] <- textvec[iter] #textlist[[iter]] <- readLines(textvec[iter])
    for(iter2 in 1:length(textlist[[iter]])) {
      docs[[iterator]] <- get.TF(textlist[[iter]][iter2],excludeWords = stopwordlist)
      iterator <- iterator + 1
    }
  }
  
  tlist <- unlist(textlist)
  tlist <- gsub("â€”","",tlist)
  tlist <- gsub("Â","",tlist)
  # This will show up rampantly if not
  tlist <- gsub("Let friends in your social network know what you are reading about","",tlist)
  tlist <- gsub("A link has been sent to your friend's email address.","",tlist) # Will need to check how to get rid of the additional information from dif sites
  #tlist <- gsub("it????Ts","",tlist)
  #tlist <- gsub("Trump????Ts","",tlist)
  #tlist <- c(textlist[[1]],textlist[[2]],textlist[[3]],textlist[[4]],textlist[[5]],textlist[[6]])
  tmatdoc <- as.DocumentTermMatrix(get.TF(tlist,excludeWords = stopwordlist))
  # Docs is now a list of 32 dataframe tf matrices
  LDAboi <- LDA(tmatdoc,k=numtopics,control=list(seed=seedset))
  suppressMessages(suppressWarnings(library(dplyr)))
  suppressMessages(suppressWarnings(library(tidytext)))
  topicss <- tidy(LDAboi, matrix="beta")
  
  
  #topn <- 10
  
  ap_top_terms <- topicss %>%
    group_by(topic) %>%
    top_n(topn, beta) %>%
    ungroup() %>%
    arrange(topic, -beta) 
  
  # if(numtopics==6) {
  # # Now make a dataframe for each top terms
  # graph.array <- list()
  # mygraphlist <- list()
  # collist <- rainbow(6)
  # collist[2] <- "#B5BC24"
  # collist[4] <- "#1BC2BD"
  # for(nn in 1:numtopics) {
  #     thisdf <- ap_top_terms[(((nn-1)*topn)+1):(nn*topn),]
  #     thisdf <- as.data.frame(thisdf[order(-thisdf$beta),])
  #     graph.array[[nn]] <- thisdf
  #     mygraphlist[[nn]] <- graph.array[[nn]] %>%
  #       ggplot(aes(x=factor(term,levels=rev(term)), y=beta,fill=max(beta))) +
  #       geom_col(show.legend = FALSE) +
  #       scale_fill_gradient(low=collist[nn],high=collist[nn]) +
  #       facet_wrap(~ topic, scales = "free") +
  #       labs(x="term",y="beta") +
  #       coord_flip() +
  #       ylim(0,0.03)
  # }
  # #Range=factor(Range, levels=new.levels)
  # # ap_top_terms <- transform(ap_top_terms, beta = reorder(beta, order(topic, beta, decreasing=T)))
  # # 
  # 
  # # mygraphic <- ap_top_terms %>%
  # #   mutate(term = reorder(term, beta)) %>%
  # #   ggplot(aes(x=factor(term,levels=rev(term)), y=beta, fill = factor(topic))) +
  # #   geom_col(show.legend = FALSE) +
  # #   facet_wrap(~ topic, scales = "free") +
  # #   labs(x="term",y="probability") +
  # #   coord_flip() +
  # #   ylim(0,0.025)
  # 
  # mygraphic <- ggarrange(mygraphlist[[1]],mygraphlist[[2]],mygraphlist[[3]],mygraphlist[[4]],mygraphlist[[5]],mygraphlist[[6]],
  #                        #labels=c("1","Topic 2","Topic 3","Topic 4","Topic 5","Topic 6"),
  #                        ncol=3,nrow=2)
  # 
  # }
  #else {
  
  ###############################################
  # Need to change this to have topics in order
  ###############################################
  
  
  
  mygraphic <- ap_top_terms %>%
    mutate(term = reorder(term, beta)) %>%
    ggplot(aes(x=term, y=beta, fill = factor(topic))) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ topic, scales = "free") +
    labs(x="term",y="probability") +
    coord_flip() +
    ylim(0,max(ap_top_terms$beta)) #}
  # Now to calculate Gini indexes
  # Find out how to get Gini's of the top 10 betas for the terms and calculate Gini over that
  Ginivectop <- c()
  ShannonICvectop <- c()
  Ginivec <- c()
  ShannonICvec <- c()
  GiniGTUnif <- c()
  ShannonGTUnif <- c()
  for(j in 1:numtopics) {
    Ginivectop[j] <- ineq(unlist(ap_top_terms[which(ap_top_terms$topic==j),3]),type="Gini")
    ShannonICvectop[j] <- ShannonIC(unlist(ap_top_terms[which(ap_top_terms$topic==j),3]))
    Ginivec[j] <- ineq(unlist(topicss[which(topicss$topic==j),3]),type="Gini")
    ShannonICvec[j] <- ShannonIC(unlist(topicss[which(topicss$topic==j),3]))
    GiniGTUnif[j] <- RenormGini(unlist(topicss[which(topicss$topic==j),3]))
    ShannonGTUnif[j] <- ShannonIC(unlist(topicss[which(topicss$topic==j),3]),gtunif = T)
  }
  # Ginivec <- c()
  # for(i in 1:numtopics) {
  #   Ginivec[i] <- ineq(unlist(topicss[which(topicss$topic==i),3]),type="Gini")
  #   
  # }
  # GiniGTUnif <- c()
  # for(ii in 1:numtopics) {
  #   GiniGTUnif[ii] <- RenormGini(unlist(topicss[which(topicss$topic==ii),3]))
  # }
  
  # Will need KL divergence 
  
  # Need to put a calculation in here for getting all the combinations of KL divergences, or KL between the number of topic models
  
  
  retlist <- list(graphic=mygraphic,LDAmodel=LDAboi,all.gini = Ginivec,top.gini=Ginivectop,gtunif.gini = GiniGTUnif,all.SIC=ShannonICvec,top.SIC=ShannonICvectop,gtunif.SIC=ShannonGTUnif,topics = topicss,top.terms = ap_top_terms)
  
  # Need to return the LDA model as well 
  
  return (retlist)
  
}


###########################
# Kendall's Tau measure
###########################


ktaub <- function(topwords, SFvec, LFvec) {
  # This function takes in 3 criteria
  # Topwords: A number > 10 that represents looking at top words in both distributions
  # SFvec, LFvec, short form and long form text vectors
  # Please note that this can be modified for any criteria, I will use pure wordfreq for now
  
  freq.SF <- sort(table(SFvec), decreasing = T)
  freq.LF <- sort(table(LFvec), decreasing = T)
  
  numconcord <- 0
  numdiscord <- 0
  numtieSF <- 0
  numtieLF <- 0
  numerrors <- 0
  for(i in 2:topwords) {
    for(j in 1:(i-1)) {
      #######################
      # Looking at sf -> lf
      #######################
      word1 <- names(freq.SF)[i]  
      word2 <- names(freq.LF)[j]
      if(word1 %in% names(freq.LF) && word2 %in% names(freq.SF) && word1 != word2) {
        # This makes sure that both words are in each set, can just set freqs to 0 if don't exist
        # Check with Peter tomorrow on seeing if this condition needs to be satisfied
        
        # Get freqs for ties
        freqSFA <- freq.SF[i]
        freqSFB <- freq.SF[which(names(freq.SF) == word2)]
        freqLFA <- freq.LF[which(names(freq.LF) == word1)]
        freqLFB <- freq.LF[j]
        
        # Get ranks for concordant, discordant
        rankSFA <- i
        rankSFB <- which(names(freq.SF) == word2)
        rankLFA <- which(names(freq.LF) == word1)
        rankLFB <- j
        
        # concordant
        if(rankSFA > rankSFB && rankLFA > rankLFB) { numconcord <- numconcord + 1 }
        else {numdiscord <- numdiscord + 1}
        # tieSF
        if(freqSFA == freqSFB && freqLFA != freqLFB) {
          numtieSF <- numtieSF + 1
          numdiscord <- numdiscord - 1
        }
        # tieLF
        if(freqSFA != freqSFB && freqLFA == freqLFB) {
          numtieLF <- numtieLF + 1
          numdiscord <- numdiscord - 1}
        
        # What do we do in case of the situation when both are ties? I won't worry about that for now
        # if(freqSFA == freqSFB && freqLFA == freqLFB) {
        #   numdiscord
        # }
        
      }
      else {
        numerrors <- numerrors + 1 
      }
    }
    # for(j in 1:(i-1)) {
    #   #######################
    #   # Looking at lf <- sf
    #   #######################
    #   
    #   # Ask Peter if I need to repeat this for hashtags not mentioned in the above
    #   
    #   word1 <- names(freq.SF)[j]
    #   word2 <- names(freq.LF)[i]
    #   if(word1 %in% names(freq.LF) && word2 %in% names(freq.SF)) {
    #     # This makes sure that both words are in each set, can just set freqs to 0 if don't exist
    #     # Check with Peter tomorrow on seeing if this condition needs to be satisfied
    #     
    #     # Get freqs for ties
    #     freqSFA <- freq.SF[i]
    #     freqSFB <- freq.SF[which(names(freq.SF) == word2)]
    #     freqLFA <- freq.LF[which(names(freq.LF) == word1)]
    #     freqLFB <- freq.LF[j]
    #     
    #     # Get ranks for concordant, discordant
    #     rankSFA <- i
    #     rankSFB <- which(names(freq.SF) == word2)
    #     rankLFA <- which(names(freq.LF) == word1)
    #     rankLFB <- j
    #     ###################################
    #     # NEED TO CHECK FOR REPEATS
    #     ###################################
    #     
    #     
    #     # concordant
    #     if(rankSFA > rankSFB && rankLFA > rankLFB) { numconcord <- numconcord + 1 }
    #     else {numdiscord <- numdiscord + 1}
    #     # tieSF
    #     if(freqSFA == freqSFB && freqLFA != freqLFB) {
    #       numtieSF <- numtieSF + 1
    #       numdiscord <- numdiscord - 1
    #     }
    #     # tieLF
    #     if(freqSFA != freqSFB && freqLFA == freqLFB) {
    #       numtieLF <- numtieLF + 1
    #       numdiscord <- numdiscord - 1}
    #     # What do we do in case of the situation when both are ties? I won't worry about that for now
    #     # if(freqSFA == freqSFB && freqLFA == freqLFB) {
    #     #   numdiscord
    #     # }
    #     
    #   }
    #  
    #}
  }
  # Calculation of ktaub
  taubnumerator <- numconcord - numdiscord
  taubdenominator <- sqrt((numconcord + numdiscord + numtieLF)*(numconcord + numdiscord + numtieSF))
  
  taub <- taubnumerator/taubdenominator
  
  return (list(taub=taub,errors=numerrors,concord=numconcord,discord = numdiscord,ties=(numtieLF+numtieSF)))
  
}

ktaubprob <- function(topwords, SFprob, LFprob) {
  # SFprob is an LDA matrix for one specific topic, on short form comparison
  # LFprob is an LDA matrix for one specific topic, on long form comparison
  
  
  
  prob.sf <- SFprob[order(SFprob$beta,decreasing = T),]
  prob.lf <- LFprob[order(LFprob$beta,decreasing = T),]
  
  numconcord <- 0
  numdiscord <- 0
  numtieSF <- 0 
  numtieLF <- 0
  numerrors <- 0
  for(i in 2:topwords) {
    for(j in 1:(i-1)) {
      word1 <- prob.sf$term[i]
      word2 <- prob.lf$term[j]
      if(word1 %in% prob.lf$term && word2 %in% prob.sf$term && word1 != word2) {
        probSFA <- prob.sf$beta[i]
        probSFB <- prob.sf$beta[which(prob.sf$term == word2)]
        probLFA <- prob.lf$beta[which(prob.lf$term == word1)]
        probLFB <- prob.lf$beta[j]
        
        rankSFA <- i
        rankSFB <- which(prob.sf$term == word2)
        rankLFA <- which(prob.lf$term == word1)
        rankLFB <- j
        
        # concordant
        if(rankSFA > rankSFB && rankLFA > rankLFB) { numconcord <- numconcord + 1 }
        else {numdiscord <- numdiscord + 1}
        # tieSF
        if(probSFA == probSFB && probLFA != probLFB) {
          numtieSF <- numtieSF + 1
          numdiscord <- numdiscord - 1
        }
        # tieLF
        if(probSFA != probSFB && probLFA == probLFB) {
          numtieLF <- numtieLF + 1
          numdiscord <- numdiscord - 1}
        
        # What do we do in case of the situation when both are ties? I won't worry about that for now
        # if(freqSFA == freqSFB && freqLFA == freqLFB) {
        #   numdiscord
        # }
        
      }
      else {
        numerrors <- numerrors + 1 
      }
    }
  }
  taubnumerator <- numconcord - numdiscord
  taubdenominator <- sqrt((numconcord + numdiscord + numtieLF)*(numconcord + numdiscord + numtieSF))
  
  taub <- taubnumerator/taubdenominator
  
  return (list(taub=taub,errors=numerrors,concord=numconcord,discord = numdiscord,ties=(numtieLF+numtieSF)))
  
}



