softmax <- function(vec) {
  v <- abs(vec)/sum(abs(vec))
  return (v)
}

simstudytm <- function(topiclist, signaltopics, nulltopics, drawtopwords, samptopwords, sampfulldocwords, numdocs, numtopics, prop.signal, prop.mix, stopwordlist,wholevocab=FALSE) {
  # This function will input a probabilistic document structure, generate documents, topic models
  # % null topic and other topics will seep in
  # topiclist is a full list of topics from a model
  # signaltopics is the topics from the model
  # drawtopwords is the number of words from the vocabulary I'll draw
  # samptopwords is amount of words in short form on average
  
  # Step 0: Get topic matrix correctly formatted (from LDA model)
  #signal.dist <- array(dim=c(length(signaltopics),drawtopwords,3))
  signal.dist <- list()
  #print(dim(signal.dist))
  for(i in 1:length(signaltopics)) {
    save1 <- topiclist[(topiclist$topic == signaltopics[i]),]
    #print(signal.dist[i,,])
    #signal.dist[i,,] <- save1[sort.int(save1$beta,decreasing = T,index.return = T)$ix[1:drawtopwords],]
    if(wholevocab==TRUE) {signal.dist[[i]] <- save1}
    else {signal.dist[[i]] <- save1[sort.int(save1$beta,decreasing = T,index.return = T)$ix[1:drawtopwords],]}
    
    #print(matrix(save1[sort.int(save1$beta,decreasing = T,index.return = T)$ix[1:drawtopwords],]),nrow=drawtopwords,ncol=3)
    #print(1)
    #print(signal.dist[i,,])
    #print(dim(signal.dist))
    #signal.dist[i,,3] <- softmax(signal.dist[i,,3])
    signal.dist[[i]][,3] <- softmax(signal.dist[[i]][,3])
    #print(2)
  }
  null.dist <- list()
  #null.dist <- array(dim=c(length(nulltopics),drawtopwords,3))
  for(j in 1:length(nulltopics)) {
    save2 <- topiclist[(topiclist$topic == nulltopics[j]),]
    #null.dist[j,,] <- save2[sort.int(save2$beta, decreasing = T, index.return = T)$ix[1:drawtopwords],]
    #null.dist[j,,] <- softmax(null.dist[j,,])
    if(wholevocab==TRUE) {null.dist[[j]] <- save2}
    else {null.dist[[j]] <- save2[sort.int(save2$beta, decreasing = T, index.return = T)$ix[1:drawtopwords],]}
    null.dist[[j]][,3] <- softmax(null.dist[[j]][,3])
  }
  
  # Step 1: Generate the documents
  longdocuments <- c()
  shortdocuments <- c()
  longwordvec <- c()
  shortwordvec <- c()
  for(n in 1:numdocs) {
    if(runif(1) < prop.signal) {
      # Draw from signal topic
      if(runif(1) < prop.mix) {
        # Mix from one or more of the signal topics
        
      }
      else {
        # Just draw words from a particular topic
        onetopic <- sample(1:length(signaltopics),1)
        # Then sample out words according to the betas respective to that topic
        # Decide the number of words in the document
        numwords <- rpois(1,lambda=sampfulldocwords)
        #numwords <- max(c(numwords,))
        #words <- sample(signal.dist[onetopic,,2],size=numwords,replace=T,prob=signal.dist[onetopic,,3])
        #TEST1 <<- signal.dist
        words <- sample(signal.dist[[onetopic]]$term,size=numwords,replace=T,prob=signal.dist[[onetopic]]$beta)
        numshortwords <- sample(seq(samptopwords,(samptopwords+10),by=1),1)
        shortwords <- words[1:numshortwords]
        
        longwordvec <- c(longwordvec,words)
        shortwordvec <- c(shortwordvec, shortwords)
      }
    }
    else {
      # Draw from null topic
      onetopic <- sample(1:length(nulltopics),1)
      
      numwords <- rnorm(1,mean=sampfulldocwords,sd=0.3*sampfulldocwords)
      numwords <- max(c(numwords, 2*samptopwords))
      #words <- sample(null.dist[onetopic,,2],size=numwords,replace=T,prob=null.dist[onetopic,,3])
      #TEST2 <<- null.dist
      words <- sample(null.dist[[onetopic]]$term,size=numwords,replace=T,prob=null.dist[[onetopic]]$beta)
      # Now need to get the short form content
      
      numshortwords <- sample(seq(samptopwords,(samptopwords+10),by=1),1)
      shortwords <- words[1:numshortwords]
      
      longwordvec <- c(longwordvec,words)
      shortwordvec <- c(shortwordvec, shortwords)
      
    }
    shortdocuments[n] <- paste(shortwords,collapse = " ")
    longdocuments[n] <- paste(words, collapse = " ")
    
    
  }
  
  # Now that the documents have been created, move to step 2 (done because of short words)
  
  # Step 3: Calculate Kendall's Tau between the documents and see how it appears
  
  ktaubrecord <- ktaub(topwords = 30,SFvec = shortwordvec, LFvec = longwordvec)$taub
  
  # Step 4: Create the topic distributions (5 topics)
  LDAnumtopic <- numtopics
  topicdist.ss <- LDAgraphic(textvec=shortdocuments,numtopics = LDAnumtopic, seedset = 5, stopwordlist = stopwordlist, topn = 10)
  topicdist.sl <- LDAgraphic(textvec=longdocuments, numtopics = LDAnumtopic, seedset = 5, stopwordlist = stopwordlist, topn = 10)
  
  # Step 5: Match the topic distributions via different methods
  
  # Match by ktaub, gini, and jsd
  iter <- 1
  
  taub.matrix <- matrix(NA,nrow=(LDAnumtopic^2),ncol=3)
  js.matrix <- matrix(NA,nrow=(LDAnumtopic^2),ncol=3)
  match.gini <- matrix(NA,nrow=LDAnumtopic,ncol=3)
  match.js <- match.gini
  match.taub <- match.gini
  
  
  indexs <- seq(1,LDAnumtopic,by=1)
  
  if(wholevocab==TRUE){compdists <- comp.Dists(topicdist.ss,topicdist.sl)}
  maxtomin.ss.gini <- order(topicdist.ss$gtunif.gini,decreasing=T) 
  maxtomin.sl.gini <- order(topicdist.sl$gtunif.gini, decreasing=T)
  ss.index <- indexs[which(topicdist.ss$gtunif.gini > 0.40)] 
  sl.index <- indexs[which(topicdist.sl$gtunif.gini > 0.52)] 
  
  for(i in 1:LDAnumtopic) {
    for(j in 1:LDAnumtopic) {
      ######################################################################
      # Match on taub
      if(wholevocab==TRUE) {taub.matrix[iter,] <- c(i,j,ktaubprob(30,compdists$d1[compdists$d1$topic==i,],
                                                                  compdists$d2[compdists$d2$topic==j,])$taub)}
      else{
        # taub.matrix[iter,] <- c(i,j,ktaubprob(30,compdists$d1[compdists$d1$topic==i,],
        #                                           compdists$d2[compdists$d2$topic==j,])$taub)
        
        taub.matrix[iter,] <- c(i,j,ktaubprob(30,topicdist.ss$topics[topicdist.ss$topics$topic==i,],
                                              topicdist.sl$topics[topicdist.sl$topics$topic==j,])$taub)
      }
      #print(77)
      ######################################################################
      
      ######################################################################
      # Match on JS divergence
      
      # js.matrix[iter,] <- c(i,j,JSdiverge(compdists$d1$beta[compdists$d1$topic==i],
      #                                         compdists$d2$beta[compdists$d2$topic==j]))
      if(wholevocab==TRUE) {js.matrix[iter,] <- c(i,j,JSdiverge(compdists$d1$beta[compdists$d1$topic==i], compdists$d2$beta[compdists$d2$topic==j]))}  
      else{  
        js.matrix[iter,] <- c(i,j,JSdiverge(topicdist.ss$topics$beta[topicdist.ss$topics$topic==i],
                                            topicdist.sl$topics$beta[topicdist.sl$topics$topic==j]))
      }
      iter <- iter + 1
      
    }
    if(wholevocab==TRUE) {match.gini[i,] <- c(maxtomin.ss.gini[i],maxtomin.sl.gini[i],JSdiverge(compdists$d1$beta[compdists$d1$topic==maxtomin.ss.gini[i]],compdists$d2$beta[compdists$d2$topic==maxtomin.sl.gini[i]]))}
    else {
      match.gini[i,] <- c(maxtomin.ss.gini[i],maxtomin.sl.gini[i],JSdiverge(topicdist.ss$topics$beta[topicdist.ss$topics$topic==maxtomin.ss.gini[i]],topicdist.sl$topics$beta[topicdist.sl$topics$topic==maxtomin.sl.gini[i]]))
    }
  }
  match.taub <- matchJS(taub.matrix,min=F)
  match.js <- matchJS(js.matrix,min=T)
  
  
  rawdocs <- list(shortdocs=shortdocuments,longdocs=longdocuments)
  topicdists <- list(ss=topicdist.ss,sl=topicdist.sl)
  matches <- list(taub=match.taub,gini=match.gini,jsd=match.js)
  
  # Step 6: Validate the process using the topic-document probabilities
  
  # Now need to run the build.contin
  
  doctopic.ss <- posterior(topicdist.ss$LDAmodel)$topics
  doctopic.sl <- posterior(topicdist.sl$LDAmodel)$topics
  #print("next")
  contin.taub <- build.contin(matched.df=match.taub, p1.doctopic = doctopic.ss, p2.doctopic = doctopic.sl) # no plurality, but threshold
  contin.gini <- build.contin(matched.df=match.gini, p1.doctopic = doctopic.ss, p2.doctopic = doctopic.sl)
  contin.jsd <- build.contin(matched.df=match.js, p1.doctopic = doctopic.ss, p2.doctopic = doctopic.sl)
  
  contintables <- list(taub=contin.taub, gini=contin.gini, jsd=contin.jsd)
  
  return (list(rawdocs=rawdocs, ktaub.pretopic=ktaubrecord, topicdists=topicdists, matches=matches, contin.tables = contintables))
  
}

# It's like 30-40 words compose short form content
# 400 words is the median of the data set for the full
# I'll use an sd of 0.3*max

test1 <- simstudytm(topiclist=LDAlist.sl[[15]]$topics,signaltopics=c(2,5,15,17,19),nulltopics=c(20), drawtopwords=40,
                    samptopwords = 39, sampfulldocwords = 400, numdocs=100, numtopics = 5, prop.signal = 0.8, prop.mix = 0, stopwordlist = swordlist)

test2 <- simstudytm(topiclist=LDAlist.sl[[15]]$topics,signaltopics=c(2,5,15,17,19),nulltopics=c(20), drawtopwords=40,
                    samptopwords = 39, sampfulldocwords = 400, numdocs=100, numtopics = 6, prop.signal = 0.8, prop.mix = 0, stopwordlist = swordlist)


test2 <- simstudytm(topiclist=LDAlist.sl[[15]]$topics,signaltopics=c(2,5,15,17,19),nulltopics=c(20), drawtopwords=40,
                    samptopwords = 39, sampfulldocwords = 400, numdocs=100, numtopics = 4, prop.signal = 0.8, prop.mix = 0, stopwordlist = swordlist)


test3 <- simstudytm(topiclist=LDAlist.sl[[15]]$topics,signaltopics=c(2,5,15,17,19),nulltopics=c(20), drawtopwords=40,
                    samptopwords = 39, sampfulldocwords = 400, numdocs=100, numtopics = 6, prop.signal = 0.8, prop.mix = 0, stopwordlist = swordlist)

#simstudytm <- function(topiclist, signaltopics, nulltopics, drawtopwords, samptopwords, sampfulldocwords, numdocs, numtopics, prop.signal, prop.mix, stopwordlist) {

test4 <- simstudytm(topiclist=LDAlist.sl[[15]]$topics,signaltopics=c(2,5,15,17,19),nulltopics=c(20), drawtopwords=40,
                    samptopwords = 39, sampfulldocwords = 400, numdocs=100, numtopics = 6, prop.signal = 0.8, prop.mix = 0, stopwordlist = swordlist,wholevocab = T)

test4 <- simstudytm(topiclist=LDAlist.sl[[15]]$topics,signaltopics=c(2,5,15,17,19),nulltopics=c(20), drawtopwords=9045,
                    samptopwords = 39, sampfulldocwords = 400, numdocs=100, numtopics = 6, prop.signal = 0.8, prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)


# Now that the wholetopic does indeed work, set up the sim study

noisevec <- c(0,0.1,0.2) # noise topic will always be number 20
# use full vocab
totaldocs <- 1000
topicvec <- c(3,5,10,19) # Noise topic again is always 20
longavglength <- c(50,100,400,1000)

simlist <- list()
iterate <-1

set.seed(5)
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      tic <- proc.time()
      simlist[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[15]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                       samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist,file="~/Documents/ATD group/LDAmodelsApril/simstudy07222019.Rdata")

# Need to do 6 topic model
mydate <- complete.fulldf %>%
  filter(date == "2019-04-15")
sixtopicmodel.ss <- LDAgraphic(tform(mydate$shortprocessed.content),numtopics=6,
                               stopwordlist = swordlist)

sixtopicmodel.sl <- LDAgraphic(tform(mydate$shortonlong.content),numtopics = 6,
                               stopwordlist = swordlist)

six.doctopic1 <- posterior(sixtopicmodel.ss$LDAmodel)$topics
six.doctopic2 <- posterior(sixtopicmodel.sl$LDAmodel)$topics

six.index <- seq(1,6,by=1)
six.savetaub <- matrix(NA,nrow=36,ncol=3)
six.savejsd <- six.savetaub
six.savejsdr <- six.savetaub

six.matchtaub <- matrix(NA,nrow=6,ncol=3)
six.matchgini <- matrix(NA,nrow=6,ncol=3)
six.matchjsd <- matrix(NA,nrow=6,ncol=3)
six.matchjsdr <- six.matchjsd

i1 <- 1 

six.compdists <- comp.Dists(sixtopicmodel.ss, sixtopicmodel.sl)
maxtomin.six.ss.gini <- order(sixtopicmodel.ss$gtunif.gini,decreasing=T) 
maxtomin.six.sl.gini <- order(sixtopicmodel.sl$gtunif.gini, decreasing=T)
ss.index <- six.index[which(sixtopicmodel.ss$gtunif.gini > 0.40)] 
sl.index <- six.index[which(sixtopicmodel.sl$gtunif.gini > 0.52)] 
for(i in 1:6) {
  for(j in 1:6) {
    ######################################################################
    # Match on taub
    six.savetaub[i1,] <- c(i,j,ktaubprob(30,six.compdists$d1[six.compdists$d1$topic==i,],
                                         six.compdists$d2[six.compdists$d2$topic==j,])$taub)
    ######################################################################
    
    ######################################################################
    # Match on JS divergence
    
    six.savejsd[i1,] <- c(i,j,JSdiverge(six.compdists$d1$beta[six.compdists$d1$topic==i],
                                        six.compdists$d2$beta[six.compdists$d2$topic==j]))
    ######################################################################
    
    ######################################################################
    # Match on JS divergence, remove non-informative topics
    
    
    
    if(is.element(i,ss.index) && is.element(j,sl.index)) {
      six.savejsdr[i1,] <- six.savejsd[i1,]
    }
    ######################################################################
    
    
    i1 <- i1 + 1
    #print(iter)
    
  }
  six.matchgini[i,] <- c(maxtomin.six.ss.gini[i],maxtomin.six.sl.gini[i],JSdiverge(six.compdists$d1$beta[six.compdists$d1$topic==maxtomin.six.ss.gini[i]],six.compdists$d2$beta[six.compdists$d2$topic==maxtomin.six.sl.gini[i]]))
}
six.matchtaub <- matchJS(six.savetaub,min=F)
six.matchjsd <- matchJS(six.savejsd,min=T)
six.matchjsdr <-  matchJS(six.savejsdr[!is.na(six.savejsdr[,1]),],min=T)

six.doublecount <- c()
# Contingency tables

sixtau <-  build.contin(matched.df=six.matchtaub, p1.doctopic=six.doctopic1, p2.doctopic=six.doctopic2, cutoff=FALSE, cutoff.num=5,
                        threshold=0.1, plurality=TRUE, plurality.tol=0.2) 

six.doublecount[1] <- length(sixtau$matches[rowSums(sixtau$matches) > 1,])

sixgini <- build.contin(matched.df=six.matchgini, p1.doctopic=six.doctopic1, p2.doctopic=six.doctopic2, cutoff=TRUE, cutoff.num=5,
                        threshold=0.1, plurality=TRUE, plurality.tol=0.2) 

six.doublecount[2] <- length(sixgini$matches[rowSums(sixgini$matches) > 1,])

sixjsd <-  build.contin(matched.df=six.matchjsd, p1.doctopic=six.doctopic1, p2.doctopic=six.doctopic2, cutoff=TRUE, cutoff.num=5,
                        threshold=0.1, plurality=TRUE, plurality.tol=0.2) 
six.doublecount[3] <- length(sixjsd$matches[rowSums(sixjsd$matches) > 1,])

# Check across all topics

apply(sixtau$table,c(1,2),sum)
apply(sixgini$table,c(1,2),sum)
apply(sixjsd$table,c(1,2),sum)

# Now do the plot for 1 of the graphics

# Run some more sims




data.ss1 <- data.frame(
  term=sixtopicmodel.ss$top.terms[1:10,2],
  beta=sixtopicmodel.ss$top.terms[1:10,3]
)

# Reorder following the value of another column:
ss.plot1 <- data.ss1 %>%
  mutate(term = fct_reorder(term, beta)) %>%
  ggplot( aes(x=term, y=beta,fill = rgb(197/255,23/255,1/255))) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme(legend.position = "none")

ez.plot <- function(indexes, colorscheme,tmodel,limity=0.035) {
  
  colnames(tmodel$top.terms) <- c(colnames(tmodel$top.terms)[c(1,2)],"probability") 
  
  data <- data.frame(
    term=tmodel$top.terms[indexes[1]:indexes[2],2],
    probability=tmodel$top.terms[indexes[1]:indexes[2],3]
  )
  
  plotthis <- data %>%
    mutate(term = fct_reorder(term, probability)) %>%
    ggplot( aes(x=term, y=probability)) +
    geom_bar(stat="identity",fill=colorscheme) +
    coord_flip() +
    theme(legend.position = "none") +
    theme(text = element_text(size=20),
          axis.text.x = element_text(angle=0, hjust=1),
          axis.text.y = element_text(angle=0)) +
    scale_y_continuous(limits=c(0,limity),breaks=c(0.01, 0.02,0.03)) +
    xlab(" ")
  
  return(plotthis)
}



ss.graphlist <- list()

ss.graphlist[[1]] <- ez.plot(indexes=c(1,10),colorscheme="#C51701",tmodel=sixtopicmodel.ss)
ss.graphlist[[2]] <- ez.plot(indexes=c(11,20),colorscheme="#960422",tmodel=sixtopicmodel.ss)
ss.graphlist[[3]] <- ez.plot(indexes=c(21,30),colorscheme="#FD8078",tmodel=sixtopicmodel.ss)
ss.graphlist[[4]] <- ez.plot(indexes=c(31,40),colorscheme="#D75F11",tmodel=sixtopicmodel.ss)
ss.graphlist[[5]] <- ez.plot(indexes=c(41,50),colorscheme="#B24E0F",tmodel=sixtopicmodel.ss)
ss.graphlist[[6]] <- ez.plot(indexes=c(51,60),colorscheme="#FDA31B",tmodel=sixtopicmodel.ss)




newplot.ss <- ggarrange(ss.graphlist[[1]],ss.graphlist[[2]],ss.graphlist[[3]],ss.graphlist[[4]],ss.graphlist[[5]],ss.graphlist[[6]],ncol=2,nrow=3,labels = c("T1","T2","T3","T4","T5","T6"))

sl.graphlist <- list()

sl.graphlist[[1]] <- ez.plot(indexes=c(1,10),colorscheme="#C51701",tmodel=sixtopicmodel.sl,limity=0.025)
sl.graphlist[[2]] <- ez.plot(indexes=c(11,20),colorscheme="#960422",tmodel=sixtopicmodel.sl,limity=0.025)
sl.graphlist[[3]] <- ez.plot(indexes=c(21,30),colorscheme="#FD8078",tmodel=sixtopicmodel.sl,limity=0.025)
sl.graphlist[[4]] <- ez.plot(indexes=c(31,40),colorscheme="#D75F11",tmodel=sixtopicmodel.sl,limity=0.025)
sl.graphlist[[5]] <- ez.plot(indexes=c(41,50),colorscheme="#B24E0F",tmodel=sixtopicmodel.sl,limity=0.025)
sl.graphlist[[6]] <- ez.plot(indexes=c(51,60),colorscheme="#FDA31B",tmodel=sixtopicmodel.sl,limity=0.025)

newplot.sl <- ggarrange(sl.graphlist[[1]],sl.graphlist[[2]],sl.graphlist[[3]],sl.graphlist[[4]],sl.graphlist[[5]],sl.graphlist[[6]],ncol=2,nrow=3,labels = c("T1","T2","T3","T4","T5","T6"))
