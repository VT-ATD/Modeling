for(k in 1:30) {
  iter <- 1
  tic <- proc.time()
  shortcomp1 <- comp.Dists(LDAlist.ss[[k]],LDAlist.sl[[k]])
  maxtomin.shortshortgini <- order(LDAlist.ss[[k]]$gtunif.gini,decreasing=T) 
  maxtomin.shortlonggini <- order(LDAlist.sl[[k]]$gtunif.gini, decreasing=T)
  shortshortindex <- indexs[which(LDAlist.ss[[k]]$gtunif.gini > 0.40)] 
  shortlongindex <- indexs[which(LDAlist.sl[[k]]$gtunif.gini > 0.52)] 
  for(i in 1:20) {
    tic1 <- proc.time()
    for(j in 1:20) {
      ######################################################################
      # Match on taub
      saveTAUBshort[k,iter,] <- c(i,j,ktaubprob(30,shortcomp1$d1[shortcomp1$d1$topic==i,],
                                                shortcomp1$d2[shortcomp1$d2$topic==j,])$taub)
      ######################################################################
      
      ######################################################################
      # Match on JS divergence
      
      saveJSshort[k,iter,] <- c(i,j,JSdiverge(shortcomp1$d1$beta[shortcomp1$d1$topic==i],
                                              shortcomp1$d2$beta[shortcomp1$d2$topic==j]))
      ######################################################################
      
      ######################################################################
      # Match on JS divergence, remove non-informative topics
      
      
      
      if(is.element(i,shortshortindex) && is.element(j,shortlongindex)) {
        saveJSRemove[k,iter,] <- saveJSshort[k,iter,]
      }
      ######################################################################
      
      
      iter <- iter + 1
      #print(iter)
      
    }
    toc1 <- proc.time()
    print(paste0("April ",k," topic ",i," Time:",round(toc1[3]-tic1[3],2)))
    # Match on gini
    matchGINI[k,i,] <- c(maxtomin.shortshortgini[i],maxtomin.shortlonggini[i],JSdiverge(shortcomp1$d1$beta[shortcomp1$d1$topic==maxtomin.shortshortgini[i]],shortcomp1$d2$beta[shortcomp1$d2$topic==maxtomin.shortlonggini[i]]))
  }
  matchTAUB[k,,] <- matchJS(saveTAUBshort[k,,],min=F)
  matchJSdiverge[k,,] <- matchJS(saveJSshort[k,,],min=T)
  #inter <- saveJSRemove[k,1:(min(which(is.na(saveJSRemove[k,,])))-1),]
  matchJSdivergeRemove[[k]] <- matchJS(saveJSRemove[k,!is.na(saveJSRemove[k,,1]),],min=T)
  toc <- proc.time()
  print(paste0("One loop has finished running, loop:",k))
  print(round(toc[3]-tic[3],2))
  timesfordays[k] <- toc[3]-tic[3]
}

# This is returning me an error right now, need to just present something

# Takes around an hour to run




# Now I need to match these kendall taus

# double checking probabilities on articles

# Look up the articles from April 15th

april15th

a <- which(april15th$url=="https://www.nytimes.com/2019/04/15/podcasts/the-daily/julian-assange-wikileaks-arrest.html")
posterior(LDAlist.ss[[15]]$LDAmodel)$topics[a,]
posterior(LDAlist.sl[[15]]$LDAmodel)$topics[a,]*100

b <- which(april15th$url=="https://www.cnn.com/2019/04/15/us/micah-herndon-boston-marathon-crawl-finish-trnd/index.html")
posterior(LDAlist.ss[[15]]$LDAmodel)$topics[b,]*100
posterior(LDAlist.sl[[15]]$LDAmodel)$topics[b,]*100

c <- which(april15th$url == "http://www.msnbc.com/rachel-maddow-show/sanders-congress-not-smart-enough-understand-trumps-tax-returns")
posterior(LDAlist.ss[[15]]$LDAmodel)$topics[c,]*100
posterior(LDAlist.sl[[15]]$LDAmodel)$topics[c,]*100

d <- which(april15th$url == "https://www.breitbart.com/economy/2019/04/14/feds-company-faked-white-collar-jobs-1900-chinese-migrants/")
posterior(LDAlist.ss[[15]]$LDAmodel)$topics[d,]*100
posterior(LDAlist.sl[[15]]$LDAmodel)$topics[d,]*100

e <- which(april15th$url == "https://www.foxnews.com/world/american-witness-describes-notre-dame-burn-all-my-insides-just-fell-apart")
posterior(LDAlist.ss[[15]]$LDAmodel)$topics[e,]*100
posterior(LDAlist.sl[[15]]$LDAmodel)$topics[e,]*100

f <- which(april15th$url == "https://www.nytimes.com/2019/04/15/nyregion/newyorktoday/nyc-news-city-hall-station.html")
sort(posterior(LDAlist.ss[[15]]$LDAmodel)$topics[f,]*100,decreasing = T)
sort(posterior(LDAlist.sl[[15]]$LDAmodel)$topics[f,]*100,decreasing = T)

# Try to get the matches printed out, get efficient tables

# Now lists all of the arrays together

allmatches <- list(taub=matchTAUB,gini=matchGINI,jsd=matchJSdiverge,jsdr=matchJSdivergeRemove)

save(allmatches,file="~/Documents/ATD group/LDAmodelsApril/allmatches.Rdata")


######################################
# Creating a 2x2 contingency table
######################################

TAUB.contin <- array(dim=c(30,20,2,2))
GINI.contin <- array(dim=c(30,20,2,2))
JSd.contin <- array(dim=c(30,20,2,2))

TAUB.day.contin <- array(dim=c(30,2,2))
GINI.day.contin <- array(dim=c(30,2,2))
JSd.day.contin <- array(dim=c(30,2,2))

totarticles <- c()
for(d in 1:30) { # For each day
  doctopic1 <- posterior(LDAlist.ss[[d]]$LDAmodel)$topics
  doctopic2 <- posterior(LDAlist.sl[[d]]$LDAmodel)$topics
  totarticles <- dim(doctopic1)[1]
  tic <- proc.time()
  for(i in 1:20) { # For each combination of topics
    # Find out the which statement
    TAUB.contin[d,i,,] <- matrix(c(0,0,0,0),ncol=2,nrow=2)
    GINI.contin[d,i,,] <- TAUB.contin[d,i,,]
    JSd.contin[d,i,,] <- TAUB.contin[d,i,,]
    TAUB.contin[d,i,1,1] <- length(which(doctopic1[,allmatches$taub[d,i,1]] >= 0.5 & doctopic2[,allmatches$taub[d,i,2]] >= 0.5))
    TAUB.contin[d,i,1,2] <- length(which(doctopic1[,allmatches$taub[d,i,1]] < 0.5 & doctopic2[,allmatches$taub[d,i,2]] >= 0.5))
    TAUB.contin[d,i,2,1] <- length(which(doctopic1[,allmatches$taub[d,i,1]] >= 0.5 & doctopic2[,allmatches$taub[d,i,2]] < 0.5))
    TAUB.contin[d,i,2,2] <- length(which(doctopic1[,allmatches$taub[d,i,1]] < 0.5 & doctopic2[,allmatches$taub[d,i,2]] < 0.5))
    
    GINI.contin[d,i,1,1] <- length(which(doctopic1[,allmatches$gini[d,i,1]] >= 0.5 & doctopic2[,allmatches$gini[d,i,2]] >= 0.5))
    GINI.contin[d,i,1,2] <- length(which(doctopic1[,allmatches$gini[d,i,1]] < 0.5 & doctopic2[,allmatches$gini[d,i,2]] >= 0.5))
    GINI.contin[d,i,2,1] <- length(which(doctopic1[,allmatches$gini[d,i,1]] >= 0.5 & doctopic2[,allmatches$gini[d,i,2]] < 0.5))
    GINI.contin[d,i,2,2] <- length(which(doctopic1[,allmatches$gini[d,i,1]] < 0.5 & doctopic2[,allmatches$gini[d,i,2]] < 0.5))
    
    JSd.contin[d,i,1,1] <- length(which(doctopic1[,allmatches$jsd[d,i,1]] >= 0.5 & doctopic2[,allmatches$jsd[d,i,2]] >= 0.5))
    JSd.contin[d,i,1,2] <- length(which(doctopic1[,allmatches$jsd[d,i,1]] < 0.5 & doctopic2[,allmatches$jsd[d,i,2]] >= 0.5))
    JSd.contin[d,i,2,1] <- length(which(doctopic1[,allmatches$jsd[d,i,1]] >= 0.5 & doctopic2[,allmatches$jsd[d,i,2]] < 0.5))
    JSd.contin[d,i,2,2] <- length(which(doctopic1[,allmatches$jsd[d,i,1]] < 0.5 & doctopic2[,allmatches$jsd[d,i,2]] < 0.5))
    
  }
  toc <- proc.time()
  TAUB.day.contin[d,,] <- apply(TAUB.contin[d,,,],c(2,3),sum)
  GINI.day.contin[d,,] <- apply(GINI.contin[d,,,],c(2,3),sum)
  JSd.day.contin[d,,] <- apply(JSd.contin[d,,,],c(2,3),sum)
  print(paste("Day",d,(toc[3]-tic[3])))
  
}

# Time to build a function for creating a contingency table for these days

build.contin <- function(matched.df, p1.doctopic, p2.doctopic, cutoff=FALSE, cutoff.num=10, cutoff.prop=FALSE, cutoff.prop.val=0.5, cutoff.value.bool=FALSE, cutoff.value=0.5,threshold=0.1, plurality=FALSE, plurality.tol=0.2) {
  # Purpose: 
  #   This function takes in a matched data frame of 2 topic distributions and outputs a contingency table of correctly matched topics to incorrect matches
  # Inputs:
  #   matched.df: A matched (ordered) data frame (K by 3) with the first two columns representing matches and the third column representing the probability
  #   p1.doctopic, p2.doctopic: A document-topic data frame (N by K) with the rows representing the topic distribution for a certain document
  #   cutoff, cutoff.num: Cutoff is a boolean representing of a number top matched from the K topics (less than K of course)
  #   cutoff.value: a value based on the third column of the data frame where it is higher than a certain value
  #   cutoff.prop, cutoff.prop.val: Cutoff.prop is a boolean representing the proportion of top matched K topics to be accounted for
  #   threshold: A value between 0 and 1 that represents the candidate set of topics to be considered when calculating numbers
  #   plurality, plurality.tol: Plurality allows for the maximum topic to be chosen in addition to a small tolerance around the maximum to be considered
  contin <- matrix(NA,nrow=2,ncol=2)
  
  
  
  K <- nrow(matched.df)
  N <- nrow(p1.doctopic)
  
  # Do I want a cutoff based on the distribution of the values in the third column? Like greater than the mean?
  
  if(cutoff) { topmatch <- cutoff.num }
  else if(cutoff.prop){ topmatch <- round(K*cutoff.prop.val,0) }
  else if(cutoff.value.bool) {topmatch <- length(which(matched.df[,3] >= cutoff.value))}
  else{ topmatch <- K  }
  
  contin <- replicate(topmatch,contin)
  contin.percent <- contin
  # Instead I'll loop through the top topic matches
  # I'm probably going to have to go by row through row
  max1 <- apply(p1.doctopic,1,which.max)
  max2 <- apply(p2.doctopic,1,which.max)
  matches <- matrix(rep(0,N*topmatch),nrow=N, ncol=topmatch)
  if(plurality==T) {
    # Make a new matrix of N by k indicating larger number
    # So find numbers that are close
    ll1 <- apply(p1.doctopic,1,max)
    ll2 <- apply(p2.doctopic,1,max)
    plu.max1 <- which(p1.doctopic >= ll1*(1-plurality.tol),arr.ind = T)
    plu.max2 <- which(p2.doctopic >= ll2*(1-plurality.tol),arr.ind = T)
    
    for(k in 1:topmatch) {
      # Filter out rows that are less than the threshold for this combination of topics
      # Need to threshold later
      # Find all the rows in each set that respect the threshold
      # The ones that are combined to hold the threshold
      # Can be double counting here
      # Need to determine a new max list
      thres <- which(p1.doctopic[,matched.df[k,1]] >= threshold | p2.doctopic[,matched.df[k,2]] >= threshold)
      matchvec <- intersect(plu.max1[plu.max1[,2]==matched.df[k,1],1],plu.max2[plu.max2[,2]==matched.df[k,2],1])
      contin[1,1,k] <- length(matchvec)
      #print(contin[1,1,k])
      contin[1,2,k] <- length(plu.max2[plu.max2[,2]==matched.df[k,2],1])-contin[1,1,k]
      #print(contin[1,2,k])
      contin[2,1,k] <- length(plu.max1[plu.max1[,2]==matched.df[k,1],1])-contin[1,1,k]
      #print(contin[2,1,k])
      contin[2,2,k] <- length(thres)-contin[1,1,k]-contin[1,2,k]-contin[2,1,k]
      #print(contin[2,2,k])
      matches[matchvec,k] <- rep(1,contin[1,1,k])
      contin.percent[,,k] <- contin[,,k]/sum(contin[,,k])
    }
  }
  else {
    for(k in 1:topmatch) {
      thres <- which(p1.doctopic[,matched.df[k,1]] > threshold | p2.doctopic[,matched.df[k,2]] > threshold )
      contin[1,1,k] <- length(which(matched.df[k,1] == max1[thres] & matched.df[k,2] == max2[thres]))
      contin[1,2,k] <- length(which(matched.df[k,1] == max1[thres] & matched.df[k,2] != max2[thres]))
      contin[2,1,k] <- length(which(matched.df[k,1] != max1[thres] & matched.df[k,2] == max2[thres]))
      contin[2,2,k] <- length(which(matched.df[k,1] != max1[thres] & matched.df[k,2] != max2[thres]))
      
      contin.percent[,,k] <- contin[,,k]/sum(contin[,,k])
    }
  }
  
  
  return (list(matches=matches,table=contin,table.percent=contin.percent))
  
  
}

hi<- build.contin(matched.df=allmatches$taub[1,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=10, cutoff.prop=FALSE, cutoff.prop.val=0.5, cutoff.value.bool=FALSE, cutoff.value=0.5,threshold=0.1, plurality=TRUE, plurality.tol=0.2)

which(rowSums(hi$matches)>1) # no dupes

# Lets go through all the days and see how much there is

save.cutoff10.plural0.2 <- array(dim=c(6,30,2,2,10))
numdoublecount <- matrix(NA,nrow=30,ncol=6)

for(i in 1:30) {
  doctopic1 <- posterior(LDAlist.ss[[i]]$LDAmodel)$topics
  doctopic2 <- posterior(LDAlist.sl[[i]]$LDAmodel)$topics
  f <- build.contin(matched.df = allmatches$taub[i,,],p1.doctopic = doctopic1, p2.doctopic = doctopic2, cutoff=TRUE, cutoff.num = 10, plurality=TRUE)
  save.cutoff10.plural0.2[1,i,,,] <- f$table
  numdoublecount[i,1] <- length(f$matches[rowSums(f$matches) > 1,])
  # Taub w/o plurality
  f2 <- build.contin(matched.df = allmatches$taub[i,,],p1.doctopic = doctopic1, p2.doctopic = doctopic2, cutoff=TRUE, cutoff.num = 10, plurality = FALSE)
  save.cutoff10.plural0.2[2,i,,,] <- f2$table
  numdoublecount[i,2] <- length(f2$matches[rowSums(f2$matches) > 1,])
  # Gini w/ plurality
  f3 <- build.contin(matched.df = allmatches$gini[i,,],p1.doctopic = doctopic1, p2.doctopic = doctopic2, cutoff=TRUE, cutoff.num = 10, plurality = TRUE)
  save.cutoff10.plural0.2[3,i,,,] <- f3$table
  numdoublecount[i,3] <- length(f3$matches[rowSums(f3$matches) > 1,])
  # Gini w/o plurality
  f4 <- build.contin(matched.df = allmatches$gini[i,,],p1.doctopic = doctopic1, p2.doctopic = doctopic2, cutoff=TRUE, cutoff.num = 10, plurality = FALSE)
  save.cutoff10.plural0.2[4,i,,,] <- f4$table
  numdoublecount[i,4] <- length(f4$matches[rowSums(f4$matches) > 1,])
  # Jsd w/ plurality
  f5 <- build.contin(matched.df = allmatches$jsd[i,,],p1.doctopic = doctopic1, p2.doctopic = doctopic2, cutoff=TRUE, cutoff.num = 10, plurality = TRUE)
  save.cutoff10.plural0.2[5,i,,,] <- f5$table
  numdoublecount[i,5] <- length(f5$matches[rowSums(f5$matches) > 1,])
  # Jsd w/o plurality
  f6 <- build.contin(matched.df = allmatches$jsd[i,,],p1.doctopic = doctopic1, p2.doctopic = doctopic2, cutoff=TRUE, cutoff.num = 10, plurality = FALSE)
  save.cutoff10.plural0.2[6,i,,,] <- f6$table
  numdoublecount[i,6] <- length(f6$matches[rowSums(f6$matches) > 1,])
  
}

apply(save.cutoff10.plural0.2[1,,,,],c(2,3),sum)
apply(save.cutoff10.plural0.2[2,,,,],c(2,3),sum)
apply(save.cutoff10.plural0.2[3,,,,],c(2,3),sum)
apply(save.cutoff10.plural0.2[4,,,,],c(2,3),sum)
apply(save.cutoff10.plural0.2[5,,,,],c(2,3),sum)
apply(save.cutoff10.plural0.2[6,,,,],c(2,3),sum)