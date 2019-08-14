##############################
####### New Sim Study  #######
##############################

noisevec <- c(0,0.1,0.2) # noise topic will always be number 20
# use full vocab
totaldocs <- 1000
topicvec <- c(3,5,10,19) # Noise topic again is always 20
longavglength <- c(50,100,400,1000)

simlist0401 <- list()
iterate <-1

set.seed(10)
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      tic <- proc.time()
      simlist0401[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[1]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0401,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr01_07232019.Rdata")

simlist0402 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0402[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[2]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0402,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr02_07232019.Rdata")

simlist0403 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0403[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[3]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0403,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr03_07232019.Rdata")

########################
# 6 topic model to show
########################

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

sixgini <- build.contin(matched.df=six.matchgini, p1.doctopic=six.doctopic1, p2.doctopic=six.doctopic2, cutoff=FALSE, cutoff.num=5,
                        threshold=0.1, plurality=TRUE, plurality.tol=0.2) 

six.doublecount[2] <- length(sixgini$matches[rowSums(sixgini$matches) > 1,])

sixjsd <-  build.contin(matched.df=six.matchjsd, p1.doctopic=six.doctopic1, p2.doctopic=six.doctopic2, cutoff=FALSE, cutoff.num=5,
                        threshold=0.1, plurality=TRUE, plurality.tol=0.2) 
six.doublecount[3] <- length(sixjsd$matches[rowSums(sixjsd$matches) > 1,])

# Check across all topics

apply(sixtau$table,c(1,2),sum)
apply(sixgini$table,c(1,2),sum)
apply(sixjsd$table,c(1,2),sum)

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

#######################################
# Sim study for many days
#######################################

simlist0403 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0403[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[3]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0403,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr03_07232019.Rdata")

print("4/4")

simlist0404 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0404[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[4]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0404,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr04_07232019.Rdata")

print("4/5")

simlist0405 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0405[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[5]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0405,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr05_07232019.Rdata")

print("4/6")

simlist0406 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0406[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[6]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0406,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr06_07232019.Rdata")

print("4/7")

simlist0407 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0407[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[7]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0407,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr07_07232019.Rdata")

print("4/8")

simlist0408 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0408[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[8]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0408,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr08_08232019.Rdata")

print("4/8")

simlist0408 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0408[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[8]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0408,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr08_07232019.Rdata")

print("4/9")

simlist0409 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0409[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[9]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0409,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr09_07232019.Rdata")

print("4/10")

simlist0410 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0410[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[10]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0410,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr10_07232019.Rdata")

print("4/11")

simlist0411 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0411[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[11]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0411,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr11_07232019.Rdata")

print("4/12")

simlist0412 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0412[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[12]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}

save(simlist0412,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr12_07232019.Rdata")

#################################################
# Looking at one of the sims
###############################################
runningsum.taub <- matrix(0,nrow=2,ncol=2)
runningsum.gini <- matrix(0,nrow=2,ncol=2)
runningsum.jsd <- matrix(0,nrow=2,ncol=2)

for(b in 1:48) {
  runningsum.taub <- runningsum.taub + apply(simlist0401[[b]]$contin.tables$taub$table,c(1,2),sum)
  runningsum.gini <- runningsum.gini + apply(simlist0401[[b]]$contin.tables$gini$table,c(1,2),sum)
  runningsum.jsd <- runningsum.jsd + apply(simlist0401[[b]]$contin.tables$jsd$table,c(1,2),sum)
}

correctmatchraw.taub <- c()
correctmatchraw.gini <- c()
correctmatchraw.jsd <- c()

correctmatchpercent.taub <- c()
correctmatchpercent.gini <- c()
correctmatchpercent.jsd <- c()


for(d in 1:48) {
  correctmatchraw.taub[d] <- apply(simlist0401[[d]]$contin.tables$taub$table,c(1,2),sum)[1,1]
  correctmatchpercent.taub[d] <- (apply(simlist0401[[d]]$contin.tables$taub$table,c(1,2),sum)/sum(simlist0401[[d]]$contin.tables$taub$table))[1,1]
  correctmatchraw.gini[d] <-  apply(simlist0401[[d]]$contin.tables$gini$table,c(1,2),sum)[1,1]
  correctmatchpercent.gini[d] <- (apply(simlist0401[[d]]$contin.tables$gini$table,c(1,2),sum)/sum(simlist0401[[d]]$contin.tables$gini$table))[1,1]
  correctmatchraw.jsd[d] <-  apply(simlist0401[[d]]$contin.tables$jsd$table,c(1,2),sum)[1,1]
  correctmatchpercent.jsd[d] <- (apply(simlist0401[[d]]$contin.tables$jsd$table,c(1,2),sum)/sum(simlist0401[[d]]$contin.tables$jsd$table))[1,1]
}
boxplot(correctmatchpercent.taub,correctmatchpercent.gini,correctmatchpercent.jsd,main="Box plots ")
# For all current sim studies

# Which indices are for each topic

bigsimlist <- list(simlist0401,simlist0402,simlist0403,simlist0404,simlist0405,simlist0406,simlist0407)



correctmatchraw.taub <- matrix(NA,nrow=7,ncol=48)
correctmatchraw.gini <- matrix(NA,nrow=7,ncol=48)
correctmatchraw.jsd <- matrix(NA,nrow=7,ncol=48)

correctmatchpercent.taub <- matrix(NA,nrow=7,ncol=48)
correctmatchpercent.gini <- matrix(NA,nrow=7,ncol=48)
correctmatchpercent.jsd <- matrix(NA,nrow=7,ncol=48)

for(a in 1:7) {
  for(b in 1:48) {
    correctmatchraw.taub[a,b] <- apply(bigsimlist[[a]][[b]]$contin.tables$taub$table,c(1,2),sum)[1,1]
    correctmatchpercent.taub[a,b] <- (apply(bigsimlist[[a]][[b]]$contin.tables$taub$table,c(1,2),sum)/sum(bigsimlist[[a]][[b]]$contin.tables$taub$table))[1,1]
    correctmatchraw.gini[a,b] <-  apply(bigsimlist[[a]][[b]]$contin.tables$gini$table,c(1,2),sum)[1,1]
    correctmatchpercent.gini[a,b] <- (apply(bigsimlist[[a]][[b]]$contin.tables$gini$table,c(1,2),sum)/sum(bigsimlist[[a]][[b]]$contin.tables$gini$table))[1,1]
    correctmatchraw.jsd[a,b] <-  apply(bigsimlist[[a]][[b]]$contin.tables$jsd$table,c(1,2),sum)[1,1]
    correctmatchpercent.jsd[a,b] <- (apply(bigsimlist[[a]][[b]]$contin.tables$jsd$table,c(1,2),sum)/sum(bigsimlist[[a]][[b]]$contin.tables$jsd$table))[1,1]
    
  }
}

boxplot(correctmatchpercent.taub[1,],correctmatchpercent.gini[1,],correctmatchpercent.jsd[1,],correctmatchpercent.taub[2,],correctmatchpercent.gini[2,],correctmatchpercent.jsd[2,],correctmatchpercent.taub[3,],correctmatchpercent.gini[3,],correctmatchpercent.jsd[3,],correctmatchpercent.taub[4,],correctmatchpercent.gini[4,],correctmatchpercent.jsd[4,],main="Box plots 4/1 - 4/4")
boxplot(correctmatchpercent.taub[5,],correctmatchpercent.gini[5,],correctmatchpercent.jsd[5,],correctmatchpercent.taub[6,],correctmatchpercent.gini[6,],correctmatchpercent.jsd[6,],correctmatchpercent.taub[7,],correctmatchpercent.gini[7,],correctmatchpercent.jsd[7,],main="Box plots 4/5 - 4/7 ")

boxplot(correctmatchraw.taub[1,],correctmatchraw.gini[1,],correctmatchraw.jsd[1,],correctmatchraw.taub[2,],correctmatchraw.gini[2,],correctmatchraw.jsd[2,],correctmatchraw.taub[3,],correctmatchraw.gini[3,],correctmatchraw.jsd[3,],correctmatchraw.taub[4,],correctmatchraw.gini[4,],correctmatchraw.jsd[4,],main="Box plots 4/1 - 4/4")
boxplot(correctmatchraw.taub[5,],correctmatchraw.gini[5,],correctmatchraw.jsd[5,],correctmatchraw.taub[6,],correctmatchraw.gini[6,],correctmatchraw.jsd[6,],correctmatchraw.taub[7,],correctmatchraw.gini[7,],correctmatchraw.jsd[7,],main="Box plots 4/5 - 4/7 ")

# Before going tonight, run the Kendall's Tau-b for the random sample with a random sample of 40 words
# Run a 20 topic model, calculate all various needed information


datevec <- as.character(seq(as.Date("2019-04-01"),by="day",length.out=30))
savetime <- c()
#savearray1 <- array(dim=c(30,length(sepvec1),5))
# This is for the short vs short long analysis
set.seed(5)

for(jj in 1:length(datevec)) {
  
  subdf <- complete.fulldf %>% 
    filter(date %in% datevec[[jj]])
  tic <- proc.time()
  
  newlongfreq <- list()
  # Need to change a's 
  # Remove quotations
  for(i in 1:nrow(subdf)) {
    # 20 words
    long <- unlist(strsplit(tform(as.character(subdf$fullprocessed.content[i])),split=" "))
    if(length(long) < 30) next
    newlongfreq[[i]] <- sample(long,size=30)
    #print(i)
  }
  sf <- unlist(strsplit(tform(as.character(subdf$shortprocessed.content)), split = " "))
  for(j in 1:length(sepvecshyam)) {
    savemat[jj,j] <- ktaub(sepvecshyam[j],sf,unlist(newlongfreq))$taub
    #print(j)
  }
  toc <- proc.time()
  print(jj)
  print(toc-tic)
  savetime[jj] <- toc-tic
}
# needs to be from 18-30
# Times take between 500-700 seconds for weekends
# Times take between 800-1000 seconds for weekdays

# Plot the Savemat

plot(sepvecshyam,savemat[1,],type="l",main="Random sample of 30 words over 30 days",ylab="Kendall's Tau-b",xlab="N top words",ylim=c(-.7,.7))
for(i in 2:30) {
  lines(sepvecshyam,savemat[i,])
}

# Now run the twenty topic model

mydate <- complete.fulldf %>%
  filter(date == "2019-04-15")
twentytopicmodel.ss <- LDAgraphic(tform(mydate$shortprocessed.content),numtopics=20,
                                  stopwordlist = swordlist)

twentytopicmodel.sl <- LDAgraphic(tform(mydate$shortonlong.content),numtopics = 20,
                                  stopwordlist = swordlist)

twenty.doctopic1 <- posterior(twentytopicmodel.ss$LDAmodel)$topics
twenty.doctopic2 <- posterior(twentytopicmodel.sl$LDAmodel)$topics

twenty.index <- seq(1,20,by=1)
twenty.savetaub <- matrix(NA,nrow=(20*20),ncol=3)
twenty.savejsd <- twenty.savetaub
twenty.savejsdr <- twenty.savetaub

twenty.matchtaub <- matrix(NA,nrow=20,ncol=3)
twenty.matchgini <- matrix(NA,nrow=20,ncol=3)
twenty.matchjsd <- matrix(NA,nrow=20,ncol=3)
twenty.matchjsdr <- twenty.matchjsd

i1 <- 1 

twenty.compdists <- comp.Dists(twentytopicmodel.ss, twentytopicmodel.sl)
maxtomin.twenty.ss.gini <- order(twentytopicmodel.ss$gtunif.gini,decreasing=T) 
maxtomin.twenty.sl.gini <- order(twentytopicmodel.sl$gtunif.gini, decreasing=T)
ss.index <- twenty.index[which(twentytopicmodel.ss$gtunif.gini > 0.40)] 
sl.index <- twenty.index[which(twentytopicmodel.sl$gtunif.gini > 0.52)] 
for(i in 1:20) {
  for(j in 1:20) {
    ######################################################################
    # Match on taub
    twenty.savetaub[i1,] <- c(i,j,ktaubprob(30,twenty.compdists$d1[twenty.compdists$d1$topic==i,],
                                            twenty.compdists$d2[twenty.compdists$d2$topic==j,])$taub)
    ######################################################################
    
    ######################################################################
    # Match on JS divergence
    
    twenty.savejsd[i1,] <- c(i,j,JSdiverge(twenty.compdists$d1$beta[twenty.compdists$d1$topic==i],
                                           twenty.compdists$d2$beta[twenty.compdists$d2$topic==j]))
    ######################################################################
    
    ######################################################################
    # Match on JS divergence, remove non-informative topics
    
    
    
    if(is.element(i,ss.index) && is.element(j,sl.index)) {
      twenty.savejsdr[i1,] <- twenty.savejsd[i1,]
    }
    ######################################################################
    
    
    i1 <- i1 + 1
    #print(iter)
    
  }
  twenty.matchgini[i,] <- c(maxtomin.twenty.ss.gini[i],maxtomin.twenty.sl.gini[i],JSdiverge(twenty.compdists$d1$beta[twenty.compdists$d1$topic==maxtomin.twenty.ss.gini[i]],twenty.compdists$d2$beta[twenty.compdists$d2$topic==maxtomin.twenty.sl.gini[i]]))
}
twenty.matchtaub <- matchJS(twenty.savetaub,min=F)
twenty.matchjsd <- matchJS(twenty.savejsd,min=T)
twenty.matchjsdr <-  matchJS(twenty.savejsdr[!is.na(twenty.savejsdr[,1]),],min=T)

twenty.doublecount <- c()
# Contingency tables

twentytau <-  build.contin(matched.df=twenty.matchtaub, p1.doctopic=twenty.doctopic1, p2.doctopic=twenty.doctopic2, cutoff=FALSE, cutoff.num=5,
                           threshold=0.1, plurality=TRUE, plurality.tol=0.2) 

twenty.doublecount[1] <- length(twentytau$matches[rowSums(twentytau$matches) > 1,])

twentygini <- build.contin(matched.df=twenty.matchgini, p1.doctopic=twenty.doctopic1, p2.doctopic=twenty.doctopic2, cutoff=FALSE, cutoff.num=5,
                           threshold=0.1, plurality=TRUE, plurality.tol=0.2) 

twenty.doublecount[2] <- length(twentygini$matches[rowSums(twentygini$matches) > 1,])

twentyjsd <-  build.contin(matched.df=twenty.matchjsd, p1.doctopic=twenty.doctopic1, p2.doctopic=twenty.doctopic2, cutoff=FALSE, cutoff.num=5,
                           threshold=0.1, plurality=TRUE, plurality.tol=0.2) 
twenty.doublecount[3] <- length(twentyjsd$matches[rowSums(twentyjsd$matches) > 1,])

# Check across all topics

apply(twentytau$table,c(1,2),sum)
apply(twentygini$table,c(1,2),sum)
apply(twentyjsd$table,c(1,2),sum)


#simlist0409 <- list()
tictic <- proc.time()
# Register my backend for parallel
library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)

iterate <-1
totalsimlist <- list()
tic <- proc.time()
# library(foreach) 
#  <- foreach(ff=8:30) %dopar% {
#   iterate <- 1
#   
# 
# iterate <- 1  

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

totalsimlist <- foreach(i = 1:length(noisevec[1:2]), .combine='list', .init=list()) %:% foreach(j = 1:length(longavglength[1:2])) %:% foreach(k = 1:length(topicvec[1:2])) %dopar% {
  simstudytm(topiclist = LDAlist.sl[[8]]$topics, signaltopics = seq(1,topicvec[k],by=1), nulltopics = 20, drawtopword = 9000, 
             samptopwords = 40, sampfulldocwords = longavglength[j], numdocs = totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
  print(iterate)
  iterate <- iterate + 1
  print(iterate)
}

#simlist0409 <- list()
tictic <- proc.time()
# Register my backend for parallel
library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)

iterate <-1
totalsimlist <- list()
tic <- proc.time()

foreach(ff = 8:8) %dopar% {
  for(i in 1:length(noisevec)) {
    for(j in 1:length(longavglength)) {
      for(k in 1:length(topicvec)) {
        print("hi")
        totalsimlist[[ff]][[iterate]] <- simstudytm(topiclist = LDAlist.sl[[ff]]$topics, signaltopics = seq(1,topicvec[k],by=1), nulltopics = 20, drawtopword = 9000, 
                                                    samptopwords = 40, sampfulldocwords = longavglength[j], numdocs = totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist)
        print("bye")
        toc <- proc.time()
        print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
        iterate <- iterate + 1
      }
    }
  }
}

stopCluster(cl)

toctoc <- proc.time()

print(toctoc[3]-tictic[3])

# If I can't do the dopar I'm just gonna do the regular

# Now just run all the simulations
timemat <- c()

# Let's just save an array with the [1,1,] table object

sim.correctarray <- array(dim=c(2,3,30,20,length(noisevec),length(longavglength),length(topicvec)))

# 2 for table and table percent
# 3 for each method
# 20 for topics
print("4/8")

simlist0408 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0408[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[8]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,8,1:topicvec[k],i,j,k] <- simlist0408[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,8,1:topicvec[k],i,j,k] <- simlist0408[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,8,1:topicvec[k],i,j,k] <- simlist0408[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,8,1:topicvec[k],i,j,k] <- simlist0408[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,8,1:topicvec[k],i,j,k] <- simlist0408[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,8,1:topicvec[k],i,j,k] <- simlist0408[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[8] <- toc[3]-tic[3]
save(simlist0408,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr08_07252019.Rdata")
rm(simlist0408)
print("4/9")

simlist0409 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0409[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[9]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,9,1:topicvec[k],i,j,k] <- simlist0409[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,9,1:topicvec[k],i,j,k] <- simlist0409[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,9,1:topicvec[k],i,j,k] <- simlist0409[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,9,1:topicvec[k],i,j,k] <- simlist0409[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,9,1:topicvec[k],i,j,k] <- simlist0409[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,9,1:topicvec[k],i,j,k] <- simlist0409[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[9] <- toc[3]-tic[3]
save(simlist0409,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr09_07252019.Rdata")
rm(simlist0409)
print("4/10")

simlist0410 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0410[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[10]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,10,1:topicvec[k],i,j,k] <- simlist0410[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,10,1:topicvec[k],i,j,k] <- simlist0410[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,10,1:topicvec[k],i,j,k] <- simlist0410[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,10,1:topicvec[k],i,j,k] <- simlist0410[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,10,1:topicvec[k],i,j,k] <- simlist0410[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,10,1:topicvec[k],i,j,k] <- simlist0410[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[10] <- toc[3]-tic[3]
save(simlist0410,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr10_07252019.Rdata")
rm(simlist0410)
print("4/11")

simlist0411 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0411[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[11]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,11,1:topicvec[k],i,j,k] <- simlist0411[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,11,1:topicvec[k],i,j,k] <- simlist0411[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,11,1:topicvec[k],i,j,k] <- simlist0411[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,11,1:topicvec[k],i,j,k] <- simlist0411[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,11,1:topicvec[k],i,j,k] <- simlist0411[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,11,1:topicvec[k],i,j,k] <- simlist0411[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[11] <- toc[3]-tic[3]
save(simlist0411,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr11_07252019.Rdata")
rm(simlist0411)
print("4/12")

simlist0412 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0412[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[12]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,12,1:topicvec[k],i,j,k] <- simlist0412[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,12,1:topicvec[k],i,j,k] <- simlist0412[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,12,1:topicvec[k],i,j,k] <- simlist0412[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,12,1:topicvec[k],i,j,k] <- simlist0412[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,12,1:topicvec[k],i,j,k] <- simlist0412[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,12,1:topicvec[k],i,j,k] <- simlist0412[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[12] <- toc[3]-tic[3]
save(simlist0412,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr12_07252019.Rdata")
rm(simlist0412)

print("4/13")

simlist0413 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0413[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[13]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,13,1:topicvec[k],i,j,k] <- simlist0413[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,13,1:topicvec[k],i,j,k] <- simlist0413[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,13,1:topicvec[k],i,j,k] <- simlist0413[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,13,1:topicvec[k],i,j,k] <- simlist0413[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,13,1:topicvec[k],i,j,k] <- simlist0413[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,13,1:topicvec[k],i,j,k] <- simlist0413[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[13] <- toc[3]-tic[3]
save(simlist0413,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr13_07252019.Rdata")
rm(simlist0413)

print("4/14")

simlist0414 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0414[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[14]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,14,1:topicvec[k],i,j,k] <- simlist0414[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,14,1:topicvec[k],i,j,k] <- simlist0414[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,14,1:topicvec[k],i,j,k] <- simlist0414[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,14,1:topicvec[k],i,j,k] <- simlist0414[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,14,1:topicvec[k],i,j,k] <- simlist0414[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,14,1:topicvec[k],i,j,k] <- simlist0414[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[14] <- toc[3]-tic[3]
save(simlist0414,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr14_07252019.Rdata")
rm(simlist0414)
print("4/15")

simlist0415 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0415[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[15]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,15,1:topicvec[k],i,j,k] <- simlist0415[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,15,1:topicvec[k],i,j,k] <- simlist0415[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,15,1:topicvec[k],i,j,k] <- simlist0415[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,15,1:topicvec[k],i,j,k] <- simlist0415[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,15,1:topicvec[k],i,j,k] <- simlist0415[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,15,1:topicvec[k],i,j,k] <- simlist0415[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[15] <- toc[3]-tic[3]
save(simlist0415,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr15_07252019.Rdata")
rm(simlist0415)
print("4/16")

simlist0416 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0416[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[16]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,16,1:topicvec[k],i,j,k] <- simlist0416[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,16,1:topicvec[k],i,j,k] <- simlist0416[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,16,1:topicvec[k],i,j,k] <- simlist0416[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,16,1:topicvec[k],i,j,k] <- simlist0416[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,16,1:topicvec[k],i,j,k] <- simlist0416[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,16,1:topicvec[k],i,j,k] <- simlist0416[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[16] <- toc[3]-tic[3]
save(simlist0416,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr16_07252019.Rdata")
rm(simlist0416)
print("4/17")

simlist0417 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0417[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[17]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,17,1:topicvec[k],i,j,k] <- simlist0417[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,17,1:topicvec[k],i,j,k] <- simlist0417[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,17,1:topicvec[k],i,j,k] <- simlist0417[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,17,1:topicvec[k],i,j,k] <- simlist0417[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,17,1:topicvec[k],i,j,k] <- simlist0417[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,17,1:topicvec[k],i,j,k] <- simlist0417[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[17] <- toc[3]-tic[3]
save(simlist0417,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr17_07252019.Rdata")
rm(simlist0417)
print("4/18")

simlist0418 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0418[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[18]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,18,1:topicvec[k],i,j,k] <- simlist0418[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,18,1:topicvec[k],i,j,k] <- simlist0418[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,18,1:topicvec[k],i,j,k] <- simlist0418[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,18,1:topicvec[k],i,j,k] <- simlist0418[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,18,1:topicvec[k],i,j,k] <- simlist0418[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,18,1:topicvec[k],i,j,k] <- simlist0418[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[18] <- toc[3]-tic[3]
save(simlist0418,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr18_07252019.Rdata")
rm(simlist0418)
print("4/19")

simlist0419 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0419[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[19]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,19,1:topicvec[k],i,j,k] <- simlist0419[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,19,1:topicvec[k],i,j,k] <- simlist0419[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,19,1:topicvec[k],i,j,k] <- simlist0419[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,19,1:topicvec[k],i,j,k] <- simlist0419[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,19,1:topicvec[k],i,j,k] <- simlist0419[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,19,1:topicvec[k],i,j,k] <- simlist0419[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[19] <- toc[3]-tic[3]
save(simlist0419,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr19_07252019.Rdata")
rm(simlist0419)
print("4/20")

simlist0420 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0420[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[20]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,20,1:topicvec[k],i,j,k] <- simlist0420[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,20,1:topicvec[k],i,j,k] <- simlist0420[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,20,1:topicvec[k],i,j,k] <- simlist0420[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,20,1:topicvec[k],i,j,k] <- simlist0420[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,20,1:topicvec[k],i,j,k] <- simlist0420[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,20,1:topicvec[k],i,j,k] <- simlist0420[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[20] <- toc[3]-tic[3]
save(simlist0420,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr20_07252019.Rdata")
rm(simlist0420)
print("4/21")

simlist0421 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0421[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[21]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,21,1:topicvec[k],i,j,k] <- simlist0421[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,21,1:topicvec[k],i,j,k] <- simlist0421[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,21,1:topicvec[k],i,j,k] <- simlist0421[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,21,1:topicvec[k],i,j,k] <- simlist0421[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,21,1:topicvec[k],i,j,k] <- simlist0421[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,21,1:topicvec[k],i,j,k] <- simlist0421[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[21] <- toc[3]-tic[3]
save(simlist0421,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr21_07252019.Rdata")
rm(simlist0421)
print("4/22")

simlist0422 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0422[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[22]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,22,1:topicvec[k],i,j,k] <- simlist0422[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,22,1:topicvec[k],i,j,k] <- simlist0422[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,22,1:topicvec[k],i,j,k] <- simlist0422[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,22,1:topicvec[k],i,j,k] <- simlist0422[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,22,1:topicvec[k],i,j,k] <- simlist0422[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,22,1:topicvec[k],i,j,k] <- simlist0422[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[22] <- toc[3]-tic[3]
save(simlist0422,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr22_07252019.Rdata")
rm(simlist0422)
print("4/23")

simlist0423 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0423[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[23]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,23,1:topicvec[k],i,j,k] <- simlist0423[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,23,1:topicvec[k],i,j,k] <- simlist0423[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,23,1:topicvec[k],i,j,k] <- simlist0423[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,23,1:topicvec[k],i,j,k] <- simlist0423[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,23,1:topicvec[k],i,j,k] <- simlist0423[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,23,1:topicvec[k],i,j,k] <- simlist0423[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[23] <- toc[3]-tic[3]
save(simlist0423,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr23_07252019.Rdata")
rm(simlist0423)
print("4/24")

simlist0424 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0424[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[24]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,24,1:topicvec[k],i,j,k] <- simlist0424[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,24,1:topicvec[k],i,j,k] <- simlist0424[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,24,1:topicvec[k],i,j,k] <- simlist0424[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,24,1:topicvec[k],i,j,k] <- simlist0424[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,24,1:topicvec[k],i,j,k] <- simlist0424[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,24,1:topicvec[k],i,j,k] <- simlist0424[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[24] <- toc[3]-tic[3]
save(simlist0424,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr24_07252019.Rdata")
rm(simlist0424)
print("4/25")

simlist0425 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0425[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[25]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,25,1:topicvec[k],i,j,k] <- simlist0425[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,25,1:topicvec[k],i,j,k] <- simlist0425[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,25,1:topicvec[k],i,j,k] <- simlist0425[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,25,1:topicvec[k],i,j,k] <- simlist0425[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,25,1:topicvec[k],i,j,k] <- simlist0425[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,25,1:topicvec[k],i,j,k] <- simlist0425[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[25] <- toc[3]-tic[3]
save(simlist0425,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr25_07252019.Rdata")
rm(simlist0425)
print("4/26")

simlist0426 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0426[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[26]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,26,1:topicvec[k],i,j,k] <- simlist0426[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,26,1:topicvec[k],i,j,k] <- simlist0426[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,26,1:topicvec[k],i,j,k] <- simlist0426[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,26,1:topicvec[k],i,j,k] <- simlist0426[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,26,1:topicvec[k],i,j,k] <- simlist0426[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,26,1:topicvec[k],i,j,k] <- simlist0426[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[26] <- toc[3]-tic[3]
save(simlist0426,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr26_07252019.Rdata")
rm(simlist0426)

print("4/27")

simlist0427 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0427[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[27]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,27,1:topicvec[k],i,j,k] <- simlist0427[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,27,1:topicvec[k],i,j,k] <- simlist0427[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,27,1:topicvec[k],i,j,k] <- simlist0427[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,27,1:topicvec[k],i,j,k] <- simlist0427[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,27,1:topicvec[k],i,j,k] <- simlist0427[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,27,1:topicvec[k],i,j,k] <- simlist0427[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[27] <- toc[3]-tic[3]
save(simlist0427,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr27_07252019.Rdata")
rm(simlist0427)
print("4/28")

simlist0428 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0428[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[28]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,28,1:topicvec[k],i,j,k] <- simlist0428[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,28,1:topicvec[k],i,j,k] <- simlist0428[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,28,1:topicvec[k],i,j,k] <- simlist0428[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,28,1:topicvec[k],i,j,k] <- simlist0428[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,28,1:topicvec[k],i,j,k] <- simlist0428[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,28,1:topicvec[k],i,j,k] <- simlist0428[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[28] <- toc[3]-tic[3]
save(simlist0428,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr28_07252019.Rdata")
rm(simlist0428)
print("4/29")

simlist0429 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0429[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[29]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,29,1:topicvec[k],i,j,k] <- simlist0429[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,29,1:topicvec[k],i,j,k] <- simlist0429[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,29,1:topicvec[k],i,j,k] <- simlist0429[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,29,1:topicvec[k],i,j,k] <- simlist0429[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,29,1:topicvec[k],i,j,k] <- simlist0429[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,29,1:topicvec[k],i,j,k] <- simlist0429[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[29] <- toc[3]-tic[3]
save(simlist0429,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr29_07252019.Rdata")
rm(simlist0429)
print("4/30")

simlist0430 <- list()
iterate <-1

set.seed(10)
tic <- proc.time()
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      simlist0430[[iterate]] <- simstudytm(topiclist = LDAlist.sl[[30]]$topics, signaltopics = seq(1,topicvec[k],by=1),nulltopics=20,drawtopwords=9000,
                                           samptopwords=40, sampfulldocwords = longavglength[j], numdocs=totaldocs, numtopics = topicvec[k], prop.signal = (1-noisevec[i]), prop.mix = 0, stopwordlist = swordlist, wholevocab = TRUE)
      sim.correctarray[1,1,30,1:topicvec[k],i,j,k] <- simlist0430[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,30,1:topicvec[k],i,j,k] <- simlist0430[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,30,1:topicvec[k],i,j,k] <- simlist0430[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,30,1:topicvec[k],i,j,k] <- simlist0430[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,30,1:topicvec[k],i,j,k] <- simlist0430[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,30,1:topicvec[k],i,j,k] <- simlist0430[[iterate]]$contin.tables$taub$table.percent[1,1,]
      toc <- proc.time()
      print(paste0("Loop: ", iterate, " out of 48, Time elapsed: ", round(toc[3]-tic[3],3)))
      iterate <- iterate + 1
    }
  }
}
timemat[30] <- toc[3]-tic[3]
save(simlist0430,file="~/Documents/ATD group/LDAmodelsApril/simstudyApr30_07252019.Rdata")
rm(simlist0430)

iterate <- 1
for(i in 1:length(noisevec)) {
  for(j in 1:length(longavglength)) {
    for(k in 1:length(topicvec)) {
      sim.correctarray[1,1,1,1:topicvec[k],i,j,k] <- simlist0401[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,1,1:topicvec[k],i,j,k] <- simlist0401[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,1,1:topicvec[k],i,j,k] <- simlist0401[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,1,1:topicvec[k],i,j,k] <- simlist0401[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,1,1:topicvec[k],i,j,k] <- simlist0401[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,1,1:topicvec[k],i,j,k] <- simlist0401[[iterate]]$contin.tables$taub$table.percent[1,1,]
      
      sim.correctarray[1,1,2,1:topicvec[k],i,j,k] <- simlist0402[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,2,1:topicvec[k],i,j,k] <- simlist0402[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,2,1:topicvec[k],i,j,k] <- simlist0402[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,2,1:topicvec[k],i,j,k] <- simlist0402[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,2,1:topicvec[k],i,j,k] <- simlist0402[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,2,1:topicvec[k],i,j,k] <- simlist0402[[iterate]]$contin.tables$taub$table.percent[1,1,]
      
      sim.correctarray[1,1,3,1:topicvec[k],i,j,k] <- simlist0403[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,3,1:topicvec[k],i,j,k] <- simlist0403[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,3,1:topicvec[k],i,j,k] <- simlist0403[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,3,1:topicvec[k],i,j,k] <- simlist0403[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,3,1:topicvec[k],i,j,k] <- simlist0403[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,3,1:topicvec[k],i,j,k] <- simlist0403[[iterate]]$contin.tables$taub$table.percent[1,1,]
      
      sim.correctarray[1,1,4,1:topicvec[k],i,j,k] <- simlist0404[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,4,1:topicvec[k],i,j,k] <- simlist0404[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,4,1:topicvec[k],i,j,k] <- simlist0404[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,4,1:topicvec[k],i,j,k] <- simlist0404[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,4,1:topicvec[k],i,j,k] <- simlist0404[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,4,1:topicvec[k],i,j,k] <- simlist0404[[iterate]]$contin.tables$taub$table.percent[1,1,]
      
      sim.correctarray[1,1,5,1:topicvec[k],i,j,k] <- simlist0405[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,5,1:topicvec[k],i,j,k] <- simlist0405[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,5,1:topicvec[k],i,j,k] <- simlist0405[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,5,1:topicvec[k],i,j,k] <- simlist0405[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,5,1:topicvec[k],i,j,k] <- simlist0405[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,5,1:topicvec[k],i,j,k] <- simlist0405[[iterate]]$contin.tables$taub$table.percent[1,1,]
      
      sim.correctarray[1,1,6,1:topicvec[k],i,j,k] <- simlist0406[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,6,1:topicvec[k],i,j,k] <- simlist0406[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,6,1:topicvec[k],i,j,k] <- simlist0406[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,6,1:topicvec[k],i,j,k] <- simlist0406[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,6,1:topicvec[k],i,j,k] <- simlist0406[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,6,1:topicvec[k],i,j,k] <- simlist0406[[iterate]]$contin.tables$taub$table.percent[1,1,]
      
      sim.correctarray[1,1,7,1:topicvec[k],i,j,k] <- simlist0407[[iterate]]$contin.tables$gini$table[1,1,]
      sim.correctarray[2,1,7,1:topicvec[k],i,j,k] <- simlist0407[[iterate]]$contin.tables$gini$table.percent[1,1,]
      sim.correctarray[1,2,7,1:topicvec[k],i,j,k] <- simlist0407[[iterate]]$contin.tables$jsd$table[1,1,]
      sim.correctarray[2,2,7,1:topicvec[k],i,j,k] <- simlist0407[[iterate]]$contin.tables$jsd$table.percent[1,1,]
      sim.correctarray[1,3,7,1:topicvec[k],i,j,k] <- simlist0407[[iterate]]$contin.tables$taub$table[1,1,]
      sim.correctarray[2,3,7,1:topicvec[k],i,j,k] <- simlist0407[[iterate]]$contin.tables$taub$table.percent[1,1,]
      
      iterate <- iterate + 1
    }
  }
}


sim.correctarray[2,1,8,1:topicvec[k],i,j,k] <- simlist0408[[iterate]]$contin.tables$gini$table.percent[1,1,]


# Need to go through each list and extract the tables like I need
# Place them in a data frame

# making the plots
# library
library(ggplot2)
# It's either 20 or 19 I'm not sure 
# create a data frame
# Only will go up to the 12th
#variety=rep(LETTERS[1:7], each=40)
#treatment=rep(c("high","low"),each=20)
#note=seq(1:280)+sample(1:150, 280, replace=T)
#data=data.frame(variety, treatment ,  note)

tunables = c("Noise: 0%, Doc Length: 50", "Noise: 0%, Doc Length: 100", "Noise: 0%, Doc Length: 400", "Noise: 0%, Doc Length 1000",
             "Noise: 10%, Doc Length: 50", "Noise: 10%, Doc Length: 100", "Noise: 10%, Doc Length: 400", "Noise: 10%, Doc Length 1000",
             "Noise: 20%, Doc Length: 50", "Noise: 20%, Doc Length: 100", "Noise: 20%, Doc Length: 400", "Noise: 20%, Doc Length 1000")
method <- rep(c("Gini","JS","Tau-b"),times=length(tunables)*12*10)
tunables <- rep(tunables, times = 360)


start <- 1
for(n in 1:12) {
  for(i in 1:length(noisevec)) {
    for(j in 1:length(longavglength)) {
      for(k in 1:3) {
        end <- start + 9
        match.percent[start:end] <-       sim.correctarray[2,2,n,1:10,i,j,1]
      }
    }
  }
}

match.percent <- # put in information from the array here
  data = data.frame(rep(tunables,each=numdays*20), match.percent,method)


bp <- ggplot(data, aes(x=tunables, y=match.percent, fill=method)) + 
  geom_boxplot() +
  facet_wrap(~tunables, scale="free")

bp + scale_fill_manual(values=c("#8B1F41", "#E87722", "#D2B48C"))

avgalldaysarray.gini <- apply(sim.correctarray[1,1,,1:10,1,1,1],MARGIN = 3,FUN = mean)
giniarray <- array(dim=c(12,10,3,4,3))
jsdarray <- giniarray
taubarray <- giniarray
for(n in 1:12) {
  for(i in 1:3) {
    for(j in 1:4) {
      for(k in 1:3) {
        giniarray[n,,i,j,k] <- sim.correctarray[1,1,n,1:10,i,j,k]
        jsdarray[n,,i,j,k] <- sim.correctarray[1,2,n,1:10,i,j,k]
        taubarray[n,,i,j,k] <- sim.correctarray[1,3,n,1:10,i,j,k]
      }
    }
  }
}
newnew <- apply(giniarray[,,1,1,1],2,mean,na.rm=T)

gini.newarray <- array(dim=c(10,3,4,3))
jsd.newarray <- array(dim=c(10,3,4,3))
taub.newarray <- array(dim=c(10,3,4,3))

for(i in 1:3) {
  for(j in 1:4) {
    for(k in 1:3) {
      gini.newarray[,i,j,k] <- apply(giniarray[,,i,j,k],2,mean,na.rm=T)
      jsd.newarray[,i,j,k] <- apply(jsdarray[,,i,j,k],2,mean,na.rm=T)
      taub.newarray[,i,j,k] <- apply(taubarray[,,i,j,k],2,mean,na.rm=T)
    }
  }
}

# Now make a data frame out of it
tunables = c("Noise: 0%, Doc Length: 50", "Noise: 0%, Doc Length: 100", "Noise: 0%, Doc Length: 400", "Noise: 0%, Doc Length 1000",
             "Noise: 10%, Doc Length: 50", "Noise: 10%, Doc Length: 100", "Noise: 10%, Doc Length: 400", "Noise: 10%, Doc Length 1000",
             "Noise: 20%, Doc Length: 50", "Noise: 20%, Doc Length: 100", "Noise: 20%, Doc Length: 400", "Noise: 20%, Doc Length 1000")
method <- rep(c("Gini","JS","Tau-b"),times=1)

#newgrid <- expand.grid(tunables,method)


#tunables <- rep(tunables, times = 360)
topic <- c("3","5","10")
newgrid <- expand.grid(rep(1:10),tunables,topic,method)
newnewgrid <- expand.grid(tunables,newgrid,rep(1:10))

#match.percent <- # put in information from the array here
data = data.frame(newgrid)
colnames(data) <- c("matches","tunables","topics","method")
data$matches <- 
  #newdata <- do.call("rbind",replicate(10,data,simplify=FALSE))
  start <- 1
for(k in 1:3) {
  for(i in 1:3) {
    for(j in 1:4) {
      end <- start + 9
      data$matches[start:end] <- gini.newarray[,i,j,k]
      start <- end + 1
    }
  }
}
start <- 361
for(k in 1:3) {
  for(i in 1:3) {
    for(j in 1:4) {
      end <- start + 9
      data$matches[start:end] <- jsd.newarray[,i,j,k]
      start <- end + 1
    }
  }
}
start <- 721
for(k in 1:3) {
  for(i in 1:3) {
    for(j in 1:4) {
      end <- start + 9
      data$matches[start:end] <- taub.newarray[,i,j,k]
      start <- end + 1
    }
  }
}

# Data should be good now

#TOPIC <- c(rep("3",12960),rep("5",12960),rep("10",12960))

bp <- ggplot(data, aes(x=topics, y=matches, fill=method)) + 
  geom_boxplot() +
  facet_wrap(~tunables, scale="free") + 
  ggtitle("Amount of document-topic matches averaged across 4/1-4/12") +
  theme(plot.title = element_text(hjust = 0.5, size=16)) 

bp + scale_fill_manual(values=c("#8B1F41", "#E87722", "#D2B48C"))