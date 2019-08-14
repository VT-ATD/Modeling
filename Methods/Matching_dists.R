#fulldata <- read.csv("AllAprilRedoProcessedData.csv",header=T)
fulldata <- read.csv("~/Documents/ATD group/ALLFILES/AllAprilRedoProcessedData.csv",header=T)
uniquedf1 <- distinct(fulldata)
uniquedf1$date <- substr(uniquedf1$publishedAt,1,10)
uniquedf1$time <- substr(uniquedf1$publishedAt,12,19)
uniquedf1$chrontime <- chron(dates=uniquedf1$date,times=uniquedf1$time, format = c(dates = "y-m-d",times="h:m:s"))
complete.fulldf <- uniquedf1[which(uniquedf1$content!="" & uniquedf1$full.content!=""),]


# Now I want to calculate all of the Kendall's Taus for every day for up to 5000n 

sepvecshyam <- c(seq(10,100,by=5),500,1000,1500)

savemat <- matrix(NA,nrow=30,ncol=length(sepvecshyam))

# Now the loop to collect all the info


# daydf <- complete.fulldf %>%
#   filter(date == "2019-04-05")
# 
# newlongfreq <- list()
# # Need to change a's 
# # Remove quotations
# for(i in 1:nrow(daydf)) {
#   # 20 words
#   long <- unlist(strsplit(tform(as.character(daydf$fullprocessed.content[i])),split=" "))
#   if(length(long) < 20) next
#   newlongfreq[[i]] <- sample(long,size=20)
#   #print(i)
# }
# sf <- unlist(strsplit(tform(as.character(daydf$shortprocessed.content)), split = " "))
# #ktaub(topwords=10,sf,unlist(newlongfreq))
# 
# 
# 
# #sepvecshyam = seq(10,150,by=5)
# saveshyamtau <- c()
# for(j in 1:length(sepvecshyam)) {
#   saveshyamtau[j] <- ktaub(sepvecshyam[j],sf,unlist(newlongfreq))$taub
#   print(j)
# }


# This will be the loop that calculates all 30 days Taubs

#sepvec1 <- seq(10,100,by=5)
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
    if(length(long) < 20) next
    newlongfreq[[i]] <- sample(long,size=20)
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

# Times take between 500-700 seconds for weekends
# Times take between 800-1000 seconds for weekdays

# Plot the Savemat

plot(sepvecshyam,savemat[1,],type="l",main="Random sample of 20 words over 30 days",ylab="Kendall's Tau-b",xlab="N top words",ylim=c(-.7,.7))
for(i in 2:30) {
  lines(sepvecshyam,savemat[i,])
}

####################################################
# Topic Matching #
####################################################

# Need to do short on short
# Short on long
# Long on long
# Long on short (but this hasn't come yet)

####################
# LDA formations #
####################
set.seed(5)
swordlist <- c("said","will","chars","-","â€™","â€œ","\"â€\u009d\"","â€”","â€”","â€¦","sunday","monday","tuesday","wednesday","thursday","friday","saturday")

april15th <- complete.fulldf %>%
  filter(date == "2019-04-15")

shortshort <- LDAgraphic(tform(april15th$shortprocessed.content),numtopics=20,
                         stopwordlist = swordlist)

shortlong <- LDAgraphic(tform(april15th$shortonlong.content),numtopics = 20,
                        stopwordlist = swordlist) # For some reason the encoding is messing up on this computer

longlong <- LDAgraphic(tform(april15th$fullprocessed.content),numtopics=20,
                       stopwordlist = swordlist)



# Part 1:
# Matching just based on ranked Gini

maxtomin.shortshortgini <- order(shortshort$gtunif.gini,decreasing=T) 
maxtomin.shortlonggini <- order(shortlong$gtunif.gini, decreasing=T)
newdists.short <- comp.Dists(shortshort,shortlong)
JSdiv.short <- c()

for(j in 1:20) {
  JSdiv.short[j] <- JSdiverge(newdists.short$d1$beta[newdists.short$d1$topic==maxtomin.shortshortgini[j]],newdists.short$d2$beta[newdists.short$d2$topic==maxtomin.shortlonggini[j]])
}

# Part 2:
# Calculate JS divergence between all combinations of topics
# Match based on JS between the distributions

saveJSshort <- matrix(NA,nrow=400,ncol=3)

# For example let's do short on short vs short on long
shortcomp <- comp.Dists(shortshort,shortlong)

######################################
# The encoding did not come through, changed stopwordlist
####################################



iter <- 1
for(i in 1:20) {
  for(j in (1:20)) {
    saveJSshort[iter,] <- c(i,j,JSdiverge(shortcomp$d1$beta[shortcomp$d1$topic==i],
                                          shortcomp$d2$beta[shortcomp$d2$topic==j]))
    iter <- iter + 1 
  }
}

saveJSshort <- saveJSshort[1:(min(which(is.na(saveJSshort)))-1),]

# Now I need to match them

# So when I match the JS diverge, look at mins for the pairs
# Then for whatever the minimum is between the two,
# Remove that match and continue down the list

# numselects <- c()
# 
# JSchange <- saveJSshort
# 
# saveJSchange <- matrix(NA,nrow=10,ncol=3)
# 
# for(k in 1:10) {
#   # Find the smallest JSdist
#   # Remove those i and j all rows from the JS change and redo
#   # Save that pair
#   minrow <- which.min(JSchange[,3])
#   saveJSchange[k,] <- JSchange[minrow,]
#   
#   # Now I need to remove the i and j rows from the JS change
#   
#   newi <- which(JSchange[,1] == JSchange[minrow,1])
#   newj <- which(JSchange[,2] == JSchange[minrow,2])
#   
#   JSchange <- JSchange[-c(newi,newj),]
# 
#   
# }


# function to match up based on JSdivergences

matchJS <- function(jsmat,min=T) {
  uniques <- c(length(unique(jsmat[,1])),length(unique(jsmat[,2])))
  nummatches <- uniques[which.min(uniques)]
  if(min==F) {
    nummatches <- uniques[which.max(uniques)]
  }
  JSchange <- jsmat
  
  saveJSchange <- matrix(NA,nrow=nummatches,ncol=3)
  
  iter <- 1
  while(nummatches > 0) {
    if(is.vector(JSchange)) {
      minrow <- 1
      saveJSchange[iter,] <- JSchange
      #print(saveJSchange[iter,])
      nummatches <- nummatches - 1
    }
    else {
      if(min==F) {
        minrow <- which.max(JSchange[,3])
      }
      else {minrow <- which.min(JSchange[,3])}
      saveJSchange[iter,] <- JSchange[minrow,]
      #print(minrow)
      
      # Now I need to remove the i and j rows from the JS change
      
      newi <- which(JSchange[,1] == JSchange[minrow,1])
      newj <- which(JSchange[,2] == JSchange[minrow,2])
      
      JSchange <- JSchange[-c(newi,newj),]
      iter <- iter + 1
      nummatches <- nummatches - 1
      #print(is.vector(JSchange))
    }
  }
  
  return (saveJSchange)
  
}


# Part 3:
# Calculate Gini for all topic dists
# For each topic distribution, find the optimal number of topics
# Utilize only those
# Compare only those selcted topics
# Match distributions based on JS divergence

# Since the Gini has been already calculated

# For shortshort I will use 0.45 to filter
# For shortlong I will use 0.52
# Then go through and do part 2
indexs <- seq(1,20,by=1)
set.seed(5)
shortshortindex <- indexs[which(shortshort$gtunif.gini > 0.40)] 
shortlongindex <- indexs[which(shortlong$gtunif.gini > 0.52)] 

saveJSJim <- matrix(NA,nrow=400,ncol=3)



iter <- 1
for(i in 1:20) {
  # Pass if something doesn't exist
  if(is.element(i,shortshortindex)) {
    for(j in (1:20)) {
      if(is.element(j,shortlongindex)) {
        saveJSJim[iter,] <- c(i,j,JSdiverge(shortcomp$d1$beta[shortcomp$d1$topic==i],
                                            shortcomp$d2$beta[shortcomp$d2$topic==j]))
        iter <- iter + 1 
      }
    }
  }
}

saveJSJim <- saveJSJim[1:(min(which(is.na(saveJSJim)))-1),]

par(mfcol=c(1,2))
hist(saveJSshort[,3],freq=F,xlim=c(0.56,0.7))
abline(v=log(2),col="red")
hist(saveJSJim[,3],freq=F,xlim=c(0.56,0.7))
abline(v=log(2),col="red")

# Show the means and variances for each 

# Then I need to actually match them based on the JSdiverge

matchedshort <- matchJS(saveJSshort)
matchedJim <- matchJS(saveJSJim)

# Make a table 

colnames(matchedshort) <- c("Truncated Topic","Full Topic","Jensen Shannon")
colnames(matchedJim) <- colnames(matchedshort)

# Make a table for matches

Ginimat <- cbind(maxtomin.shortshortgini,maxtomin.shortlonggini,JSdiv.short)
colnames(Ginimat) <- colnames(matchedJim)
library(data.table)
foo.dt <- data.table(Ginimat, key="Jensen Shannon")



kable(matchedshort)
kable(matchedJim)

par(mfcol=c(1,3))
hist(JSdiv.short,xlim=c(0.58,0.7),freq=F,main="Matching on Gini",xlab="Jensen-Shannon divergence",ylim=c(0,35),cex.lab=1.3,cex.axis=1.5,cex.main=2)
abline(v=log(2),col="red")
hist(matchedshort[,3],xlim=c(0.58,0.7),freq=F,main="No filtering topics",xlab="Jensen-Shannon divergence",ylim=c(0,35),cex.lab=1.3,cex.axis=1.5,cex.main=2)
abline(v=log(2),col="red")
hist(matchedJim[,3],xlim=c(0.58,0.7),freq=F,main="Filtering topics",xlab="Jensen-Shannon divergence",ylim=c(0,35),cex.lab=1.3,cex.axis=1.5,cex.main=2)
abline(v=log(2),col="red")
# Then make a histogram comparing both



######################################################
# Short versus long long
######################################################

#############################################
# This was all previous analysis 
# Down below is for doing 30 topic models
##############################################
set.seed(5)
swordlist <- c("said","will","chars","-","â€™","â€œ","\"â€\u009d\"","â€”","â€”","â€¦","sunday","monday","tuesday","wednesday","thursday","friday","saturday","can","new","april")

LDAlist.ss <- list()
LDAlist.sl <- list()
LDAlist.ll <- list()
for(i in 1:length(datevec)) {
  tic <- proc.time()
  thedate <- complete.fulldf %>%
    filter(date == datevec[i])
  
  LDAlist.ss[[i]] <- LDAgraphic(tform(thedate$shortprocessed.content),numtopics=20,
                                stopwordlist = swordlist)
  
  LDAlist.sl[[i]] <- LDAgraphic(tform(thedate$shortonlong.content),numtopics = 20,
                                stopwordlist = swordlist) # For some reason the encoding is messing up on this computer
  
  LDAlist.ll[[i]] <- LDAgraphic(tform(thedate$fullprocessed.content),numtopics=20,
                                stopwordlist = swordlist)
  toc <- proc.time()
  print(i)
  print(toc[3]-tic[3])
}