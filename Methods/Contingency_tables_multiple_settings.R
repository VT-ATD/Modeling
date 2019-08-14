pluralityvec <- c(0,0.1,0.2,0.3,0.4,0.5)
thresholdvec <- c(0,0.05,0.1,0.15,0.2,0.25)



# what elements stand for in this array
# First: Contingency table and the percent contingency table
# Second: Which method, ktau, gini, or jsd
# Third: Number of days
# Fourth: Amount of plurality conditions
# Fifth: Amount of threshold conditions
# Sixth/Seventh: The 2x2 contingency table
# Eighth: The amount of topics matched

continarray2topic <- array(dim=c(4,2,30,length(pluralityvec),length(thresholdvec),2,2,2))
continarray4topic <- array(dim=c(4,2,30,length(pluralityvec),length(thresholdvec),2,2,4))
continarray6topic <- array(dim=c(4,2,30,length(pluralityvec),length(thresholdvec),2,2,6))
continarray10topic <- array(dim=c(4,2,30,length(pluralityvec),length(thresholdvec),2,2,10))
continarray20topic <- array(dim=c(4,2,30,length(pluralityvec),length(thresholdvec),2,2,20))

doublecount <- array(dim=c(4,30,length(pluralityvec),length(thresholdvec),5))

# The 20 topic [4,,,,,] will most certainly be NA's

iteration <- 1
for(d in 1:30) {
  doctopic1 <- posterior(LDAlist.ss[[d]]$LDAmodel)$topics
  doctopic2 <- posterior(LDAlist.sl[[d]]$LDAmodel)$topics
  for(p in 1:length(pluralityvec)) {
    for(t in 1:length(thresholdvec)) {
      tic <- proc.time()
      a.tau <-  build.contin(matched.df=allmatches$taub[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray2topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray2topic[1,1,d,p,t,,,] <- a.tau$table
      continarray2topic[1,2,d,p,t,,,] <- a.tau$table.percent
      doublecount[1,d,p,t,1] <- length(a.tau$matches[rowSums(a.tau$matches) > 1,])
      
      a.gini <- build.contin(matched.df=allmatches$gini[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray2topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray2topic[2,1,d,p,t,,,] <- a.gini$table
      continarray2topic[2,2,d,p,t,,,] <- a.gini$table.percent
      doublecount[2,d,p,t,1] <- length(a.gini$matches[rowSums(a.gini$matches) > 1,])
      
      a.jsd <-  build.contin(matched.df=allmatches$jsd[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray2topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray2topic[3,1,d,p,t,,,] <- a.jsd$table
      continarray2topic[3,2,d,p,t,,,] <- a.jsd$table.percent
      doublecount[3,d,p,t,1] <- length(a.jsd$matches[rowSums(a.jsd$matches) > 1,])
      
      if((dim(continarray2topic)[8]) < nrow(allmatches$jsdr[[d]])) {
        a.jsdr <- build.contin(matched.df=allmatches$jsdr[[d]], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray2topic)[8],
                               threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
        continarray2topic[4,1,d,p,t,,,] <- a.jsdr$table
        continarray2topic[4,2,d,p,t,,,] <- a.jsdr$table.percent
        doublecount[4,d,p,t,1] <- length(a.jsdr$matches[rowSums(a.jsdr$matches) > 1,])
      }
      ###############################
      # Now for 4 topic cutoff
      ###############################
      a.tau <-  build.contin(matched.df=allmatches$taub[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray4topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray4topic[1,1,d,p,t,,,] <- a.tau$table
      continarray4topic[1,2,d,p,t,,,] <- a.tau$table.percent
      doublecount[1,d,p,t,2] <- length(a.tau$matches[rowSums(a.tau$matches) > 1,])
      
      a.gini <- build.contin(matched.df=allmatches$gini[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray4topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray4topic[2,1,d,p,t,,,] <- a.gini$table
      continarray4topic[2,2,d,p,t,,,] <- a.gini$table.percent
      doublecount[2,d,p,t,2] <- length(a.gini$matches[rowSums(a.gini$matches) > 1,])
      
      a.jsd <-  build.contin(matched.df=allmatches$jsd[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray4topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray4topic[3,1,d,p,t,,,] <- a.jsd$table
      continarray4topic[3,2,d,p,t,,,] <- a.jsd$table.percent
      doublecount[3,d,p,t,2] <- length(a.jsd$matches[rowSums(a.jsd$matches) > 1,])
      
      if((dim(continarray4topic)[8]) < nrow(allmatches$jsdr[[d]])) {
        a.jsdr <- build.contin(matched.df=allmatches$jsdr[[d]], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray4topic)[8],
                               threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
        continarray4topic[4,1,d,p,t,,,] <- a.jsdr$table
        continarray4topic[4,2,d,p,t,,,] <- a.jsdr$table.percent
        doublecount[4,d,p,t,2] <- length(a.jsdr$matches[rowSums(a.jsdr$matches) > 1,])
      }
      ###############################
      # Now for 6 topic cutoff
      ###############################
      a.tau <-  build.contin(matched.df=allmatches$taub[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray6topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray6topic[1,1,d,p,t,,,] <- a.tau$table
      continarray6topic[1,2,d,p,t,,,] <- a.tau$table.percent
      doublecount[1,d,p,t,3] <- length(a.tau$matches[rowSums(a.tau$matches) > 1,])
      
      a.gini <- build.contin(matched.df=allmatches$gini[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray6topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray6topic[2,1,d,p,t,,,] <- a.gini$table
      continarray6topic[2,2,d,p,t,,,] <- a.gini$table.percent
      doublecount[2,d,p,t,3] <- length(a.gini$matches[rowSums(a.gini$matches) > 1,])
      
      a.jsd <-  build.contin(matched.df=allmatches$jsd[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray6topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray6topic[3,1,d,p,t,,,] <- a.jsd$table
      continarray6topic[3,2,d,p,t,,,] <- a.jsd$table.percent
      doublecount[3,d,p,t,3] <- length(a.jsd$matches[rowSums(a.jsd$matches) > 1,])
      
      if((dim(continarray6topic)[8]) < nrow(allmatches$jsdr[[d]])) {
        a.jsdr <- build.contin(matched.df=allmatches$jsdr[[d]], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray6topic)[8],
                               threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
        continarray6topic[4,1,d,p,t,,,] <- a.jsdr$table
        continarray6topic[4,2,d,p,t,,,] <- a.jsdr$table.percent
        doublecount[4,d,p,t,3] <- length(a.jsdr$matches[rowSums(a.jsdr$matches) > 1,])
      }
      ###############################
      # Now for 10 topic cutoff
      ###############################
      a.tau <-  build.contin(matched.df=allmatches$taub[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray10topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray10topic[1,1,d,p,t,,,] <- a.tau$table
      continarray10topic[1,2,d,p,t,,,] <- a.tau$table.percent
      doublecount[1,d,p,t,4] <- length(a.tau$matches[rowSums(a.tau$matches) > 1,])
      
      a.gini <- build.contin(matched.df=allmatches$gini[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray10topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray10topic[2,1,d,p,t,,,] <- a.gini$table
      continarray10topic[2,2,d,p,t,,,] <- a.gini$table.percent
      doublecount[2,d,p,t,4] <- length(a.gini$matches[rowSums(a.gini$matches) > 1,])
      
      a.jsd <-  build.contin(matched.df=allmatches$jsd[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray10topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray10topic[3,1,d,p,t,,,] <- a.jsd$table
      continarray10topic[3,2,d,p,t,,,] <- a.jsd$table.percent
      doublecount[3,d,p,t,4] <- length(a.jsd$matches[rowSums(a.jsd$matches) > 1,])
      
      if((dim(continarray10topic)[8]) < nrow(allmatches$jsdr[[d]])) {
        a.jsdr <- build.contin(matched.df=allmatches$jsdr[[d]], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray10topic)[8],
                               threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
        continarray10topic[4,1,d,p,t,,,] <- a.jsdr$table
        continarray10topic[4,2,d,p,t,,,] <- a.jsdr$table.percent
        doublecount[4,d,p,t,4] <- length(a.jsdr$matches[rowSums(a.jsdr$matches) > 1,])
      }
      ###############################
      # Now for 20 topic cutoff
      ###############################
      a.tau <-  build.contin(matched.df=allmatches$taub[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray20topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray20topic[1,1,d,p,t,,,] <- a.tau$table
      continarray20topic[1,2,d,p,t,,,] <- a.tau$table.percent
      doublecount[1,d,p,t,5] <- length(a.tau$matches[rowSums(a.tau$matches) > 1,])
      
      a.gini <- build.contin(matched.df=allmatches$gini[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray20topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray20topic[2,1,d,p,t,,,] <- a.gini$table
      continarray20topic[2,2,d,p,t,,,] <- a.gini$table.percent
      doublecount[2,d,p,t,5] <- length(a.gini$matches[rowSums(a.gini$matches) > 1,])
      
      a.jsd <-  build.contin(matched.df=allmatches$jsd[d,,], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray20topic)[8],
                             threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
      continarray20topic[3,1,d,p,t,,,] <- a.jsd$table
      continarray20topic[3,2,d,p,t,,,] <- a.jsd$table.percent
      doublecount[3,d,p,t,5] <- length(a.jsd$matches[rowSums(a.jsd$matches) > 1,])
      
      if((dim(continarray20topic)[8]) < nrow(allmatches$jsdr[[d]])) {
        a.jsdr <- build.contin(matched.df=allmatches$jsdr[[d]], p1.doctopic=doctopic1, p2.doctopic=doctopic2, cutoff=TRUE, cutoff.num=dim(continarray20topic)[8],
                               threshold=thresholdvec[t], plurality=TRUE, plurality.tol=pluralityvec[p]) 
        continarray20topic[4,1,d,p,t,,,] <- a.jsdr$table
        continarray20topic[4,2,d,p,t,,,] <- a.jsdr$table.percent
        doublecount[4,d,p,t,5] <- length(a.jsdr$matches[rowSums(a.jsdr$matches) > 1,])
      }
      toc <- proc.time()
      print(paste0("Current loop: ", iteration, " out of ",  30*length(pluralityvec)*length(thresholdvec)," Day ",d," Plurality: ",pluralityvec[p], " Threshold: ",thresholdvec[t], " "))
      print(paste0("Time elapsed ", round(toc[3]-tic[3],3)))
      iteration <- iteration + 1
    }
  }
}
# This took around 15 minutes
# Around 8.1% of the time we are double counts, at most by 140 (6 times)


# 2   4   6   10   12  18  20    24  30  40    50  60  80  100 120  140 
# 18  102 210 414  42  12  540   6   12  126   6   84  12  12  24   6

# Summary without 0's (1626 of them)
# Min.   1st Qu.  Median    Mean    3rd Qu.    Max. 
# 2.00   10.00    20.00     20.99   20.00      140.00

# Summary with 0's (21600-1656)
# Min.    1st Qu.  Median    Mean   3rd Qu.    Max.       NA's 
# 0.000   0.000    0.000    1.711   0.000      140.000    1656 


# 0      2     4     6     10     12    18    20      24    30    40      50    60    80    100   120    140 
# 18318  18    102   210   414    42    12    540     6     12    126     6     84    12    12    24     6

allcontin20topictables <- list(t2=continarray2topic,t4=continarray4topic,t6=continarray6topic,t10=continarray10topic,t20=continarray20topic,doublecount=doublecount)
save(allcontin20topictables,file="~/Documents/ATD group/LDAmodelsApril/20topic30daymatchingcontingencytables.Rdata")
