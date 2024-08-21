laggedJrqa <- function(maxlag, ts1, ts2, delay, embed, rescale,
                        radius, normalize, mindiagline, minvertline, tw, method) {
  # Note:
  # This function wraps the crqa()-function from 'crqa' package.
  # In order to run the function, the packages 'crqa'
  # as well as its dependencies, need to be installed and loaded.
  # This function suits only 3-variate time series and computes only RR.
  # The authors give no warranty for the correct functioning of the software and cannot be held legally accountable.

# load libraries
require(crqa)

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
results <- matrix(nrow = dim(lagList)[1], ncol = (1+dim(lagList)[2]-1))
results[,2:(2+(dim(lagList)[2]-2))] <- lagList[,1:(dim(lagList)[2]-1)]
# run lagged jrqa
for(i in 1:dim(lagList)[1]) {
  # Generate 3 unidimensional RP
  temp_rqa1 <- crqa(ts1=ts1[,1], ts2=ts1[,1], delay = delay, embed = embed,
                   rescale = rescale, radius = radius, normalize = normalize,
                   mindiagline = mindiagline, minvertline = minvertline,
                   tw = tw, method = "rqa")
  temp_rqa2 <- crqa(ts1=ts1[,2], ts2=ts1[,2], delay = delay, embed = embed,
                    rescale = rescale, radius = radius, normalize = normalize,
                    mindiagline = mindiagline, minvertline = minvertline,
                    tw = tw, method = "rqa")
  temp_rqa3 <- crqa(ts1=ts1[,3], ts2=ts1[,3], delay = delay, embed = embed,
                    rescale = rescale, radius = radius, normalize = normalize,
                    mindiagline = mindiagline, minvertline = minvertline,
                    tw = tw, method = "rqa")
  RP1=as.matrix(temp_rqa1$RP)
  RP2=as.matrix(temp_rqa2$RP)
  RP3=as.matrix(temp_rqa3$RP)
  # Apply the lags
  tempRP1=RP1[(1+lagList[i,1]):(dim(RP1)[1]-maxlag+lagList[i,1]),(1+lagList[i,1]):(dim(RP1)[1]-maxlag+lagList[i,1])]
  tempRP2=RP2[(1+lagList[i,2]):(dim(RP2)[1]-maxlag+lagList[i,2]),(1+lagList[i,2]):(dim(RP2)[1]-maxlag+lagList[i,2])]
  tempRP3=RP3[(1+lagList[i,3]):(dim(RP3)[1]-maxlag+lagList[i,3]),(1+lagList[i,3]):(dim(RP3)[1]-maxlag+lagList[i,3])]
  # Multiply to create JRQA
  tempRP=tempRP1*tempRP2*tempRP3
  # Compute RR
  results[i,1]=sum(tempRP)/(dim(tempRP)[1]*dim(tempRP)[2])
} 

# convert results to data frame
results <- as.data.frame(results)

# generate and add labels
newLabels <- c("RR")
j <- 0
for(i in 2:(2+(dim(lagList)[2]-2))) {
  j <- j+1
  newLabels[i] <- paste("ts",as.character(j),sep="")
}
colnames(results) <- newLabels

# sort data frame by RR
results <- results[order(results$RR, decreasing = TRUE),]

# return results
return(results)
}



