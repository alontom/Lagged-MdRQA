laggedMdrqa <- function(maxlag, ts1, ts2, delay, embed, rescale,
                        radius, normalize, mindiagline, minvertline, tw, method) {
  # Note:
  # This function wraps the crqa()-function from 'crqa' package.
  # In order to run the function, the packages 'crqa' and 'tcpl',
  # as well as their dependencies, need to be installed and loaded.

# load libraries
library(tcpl)
library(crqa)

# infer parameters to create empty matrix
noTs <- dim(ts1)[2]
lagList <- matrix(0, ncol = noTs+1, nrow = (maxlag+1)^noTs)

# compute list of lag combinations
for(i in 1:noTs) {
    lagList[,i] <-rep(0:maxlag, each = (maxlag+1)^(i-1), (maxlag+1)^(noTs-i))
}

# equate lags on lag0
for(i in 1:(maxlag+1)^noTs) {
  lagList[i,1:noTs] <- lagList[i,1:noTs]-min(lagList[i,1:noTs])
}

# compute lag identifier
for(i in 1:noTs) {
  lagList[,noTs+1] <- lagList[,noTs+1] + lagList[,i]*100^(noTs-i)
}

# discard all lags that are not unique
lagList <- subset(lagList, duplicated(lagList[,noTs+1]) == FALSE)

# create results matrix and add lag parameters
results <- matrix(nrow = dim(lagList)[1], ncol = (9+dim(lagList)[2]-1))
results[,10:(10+(dim(lagList)[2]-2))] <- lagList[,1:(dim(lagList)[2]-1)]

# run lagged mdrqa
for(i in 1:dim(lagList)[1]) {
  
  # create temporary data matrix
  temp_ts1 <- matrix(ncol = dim(ts1)[2], nrow = length((1+lagList[i,1]):(dim(ts1)[1]-maxlag+lagList[i,1])))
  temp_ts2 <- matrix(ncol = dim(ts2)[2], nrow = length((1+lagList[i,1]):(dim(ts2)[1]-maxlag+lagList[i,1])))
  
  # construct lagged time series
  for(j in 1:dim(ts1)[2]) {
    temp_ts1[,j] <- ts1[(1+lagList[i,j]):(dim(ts1)[1]-maxlag+lagList[i,j]),j]
  }
  for(j in 1:dim(ts2)[2]) {
    temp_ts2[,j] <- ts2[(1+lagList[j,1]):(dim(ts2)[1]-maxlag+lagList[j,1]),j]
  }
  
  # rund mdrqa
        temp_res <- crqa(ts1=temp_ts1, ts2=temp_ts2, delay = delay, embed = embed,
                         rescale = rescale, radius = radius, normalize = normalize,
                         mindiagline = mindiagline, minvertline = minvertline,
                         tw = tw, method = "mdcrqa")
        
        # store recurrence measures on each iteration
        results[i,1:9] <- unlist(temp_res[1:9])
} 

# convert results to data frame
results <- as.data.frame(results)

# generate and add labels
newLabels <- c("RR","DER","NRLINE","maxL","L","ENTR","rENTR","LAM","TT")
j <- 0
for(i in 10:(10+(dim(lagList)[2]-2))) {
  j <- j+1
  newLabels[i] <- paste("ts",as.character(j),sep="")
}
colnames(results) <- newLabels

# sort data frame by RR
results <- results[order(results$RR, decreasing = TRUE),]

# return results
return(results)
}



