# The authors give no warranty for the correct functioning of the software and cannot be held legally accountable.

### Correlated Noise

# Clean workspace
rm(list=ls())

# Load libraries
require(crqa)

# Set wd
setwd("...")

# Noise - MdRQA

RECBest_cor_md=c()
noise=rnorm(510)
for(lagpar in 0:5) {
  for(it in 1:100) {
    ts1 <- noise[1:(500+lagpar)]
    
    # Create a lagged version of ts1 with added noise
    ts2 <- scale(ts1[(lagpar+1):length(ts1)] + rnorm(500)*sd(ts1)/10)
    ts3 <- scale(ts1[(lagpar+1):length(ts1)] + rnorm(500)*sd(ts1)/10)
    ts1 <- scale(ts1 + rnorm(length(ts1))*sd(ts1)/10)
    
    # Merge the three time series
    data <- as.matrix(cbind(ts1[1:500],ts2,ts3))
    res=laggedMdrqa(maxlag=lagpar+2, ts1=data, ts2=data, delay=1, embed=1, rescale=0,
                    radius=0.6, normalize=0, mindiagline=2, minvertline=2, tw=1, method='mdcrqa')
    
    # Find the best lag combination based on RR
    RECBest_cor_md=rbind(RECBest1,cbind(lagpar,it,res$RR[1],res$ts1[1],res$ts2[1],res$ts3[1]))
    print(cbind(lagpar,'x',it))
  }
  saveRDS(RECBest_cor_md,'MdRQA_noEmbnoise_RR.RDS')}


# Noise - JRQA

RECBest_cor_j=c()
noise=rnorm(510)
for(lagpar in 0:5) {
  for(it in 1:100) {
    ts1 <- noise[1:(500+lagpar)]
    
    # Create a lagged version of ts1 with added noise
    ts2 <- scale(ts1[(lagpar+1):length(ts1)] + rnorm(500)*sd(ts1)/10)
    ts3 <- scale(ts1[(lagpar+1):length(ts1)] + rnorm(500)*sd(ts1)/10)
    ts1 <- scale(ts1 + rnorm(length(ts1))*sd(ts1)/10)
    
    # Merge the three time series
    data <- as.matrix(cbind(ts1[1:500],ts2,ts3))
    res=laggedJrqa(maxlag=lagpar+2, ts1=data, ts2=data, delay=1, embed=1, rescale=0,
                   radius=0.6, normalize=0, mindiagline=2, minvertline=2, tw=1, method='mdcrqa')
    
    # Find the best lag combination based on RR
    RECBest_cor_j=rbind(RECBest3,cbind(lagpar,it,res$RR[1],res$ts1[1],res$ts2[1],res$ts3[1]))
    print(cbind(lagpar,'x',it))
  }
  saveRDS(RECBest_cor_j,'JRQA_noEmbnoise_RR_2.RDS')}



### Uncorrelated Noise

# Clean workspace
rm(list=ls())

# Load libraries
require(crqa)

# Set wd
setwd("...")

#Noise - MdRQA

RECBest_un_md=c()
for(lagpar in 0:5) {
  for(it in 1:100) {
    ts1 <- scale(rnorm(500))
    ts2 <- scale(rnorm(500))
    ts3 <- scale(rnorm(500))
    
    # Merge the three time series
    data <- as.matrix(cbind(ts1,ts2,ts3))
    res=laggedMdrqa(maxlag=lagpar+2, ts1=data, ts2=data, delay=1, embed=1, rescale=0,
                    radius=0.6, normalize=0, mindiagline=2, minvertline=2, tw=1, method='mdcrqa')
    
    # Find the best lag combination based on RR
    RECBest_un_md=rbind(RECBest1,cbind(lagpar,it,res$RR[1],res$ts1[1],res$ts2[1],res$ts3[1]))
    print(cbind(lagpar,'x',it))
  }
  saveRDS(RECBest_un_md,'MdRQA_noEmbnoise_RR_no_relationship_2.RDS')
}

# Noise - JRQA

RECBest_un_j=c()
for(lagpar in 0:5) {
  for(it in 1:100) {
    ts1 <- scale(rnorm(500))
    ts2 <- scale(rnorm(500))
    ts3 <- scale(rnorm(500))
    
    # Merge the three time series
    data <- as.matrix(cbind(ts1,ts2,ts3))
    res=laggedJrqa(maxlag=lagpar+2, ts1=data, ts2=data, delay=1, embed=1, rescale=0,
                   radius=0.6, normalize=0, mindiagline=2, minvertline=2, tw=1, method='mdcrqa')
    
    # Find the best lag combination based on RR
    RECBest_un_j=rbind(RECBest3,cbind(lagpar,it,res$RR[1],res$ts1[1],res$ts2[1],res$ts3[1]))
    print(cbind(lagpar,'x',it))
  }
  saveRDS(RECBest_un_j,'JRQA_noEmbnoise_RR_no_relationship_2.RDS')
}


### Sine Wave - MdRQA

# Unembedded

RECBest_emb_md=c()
t=seq(0,4+36/500,4/500)
sine=sin(t*2*pi)
for(lagpar in 0:5) {
  for(it in 1:500) {
    ts1 <- sine[1:(500+lagpar)]
    
    # Create a lagged version of ts1 with added noise
    ts2 <- scale(ts1[(lagpar+1):length(ts1)] + rnorm(500)*sd(ts1)/5)
    ts3 <- scale(ts1[(lagpar+1):length(ts1)] + rnorm(500)*sd(ts1)/5)
    
    # Merge the three time series
    data <- as.matrix(cbind(ts1[1:500],ts2,ts3))
    res=laggedMdrqa(maxlag=lagpar+2, ts1=data, ts2=data, delay=1, embed=1, rescale=0,
                    radius=0.5, normalize=0, mindiagline=1, minvertline=1, tw=1, method='mdcrqa')
    
    # Find the best lag combination based on RR
    RECBest_emb_md=rbind(RECBest,cbind(lagpar,it,res$RR[1],res$ts1[1],res$ts2[1],res$ts3[1]))
    print(cbind(lagpar,'x',it))
  }
  saveRDS(RECBest_emb_md,'MdRQA_noEmb_RR.RDS')}

# Embedded

RECBest_unemb_md=c()
t=seq(0,4+36/500,4/500)
sine=sin(t*2*pi)
for(lagpar in 0:5) {
  for(it in 1:500) {
    ts1 <- sine[1:(500+lagpar)]
    
    # Create a lagged version of ts1 with added noise
    ts2 <- scale(ts1[(lagpar+1):length(ts1)] + rnorm(500)*sd(ts1)/5)
    ts3 <- scale(ts1[(lagpar+1):length(ts1)] + rnorm(500)*sd(ts1)/5)
    
    # Merge the three time series
    data <- as.matrix(cbind(ts1[1:500],ts2,ts3))
    res=laggedMdrqa(maxlag=lagpar+2, ts1=data, ts2=data, delay=30, embed=2, rescale=0,
                    radius=1.2, normalize=0, mindiagline=1, minvertline=1, tw=1, method='mdcrqa')
    
    # Find the best lag combination based on RR
    RECBest_unemb_md=rbind(RECBest,cbind(lagpar,it,res$RR[1],res$ts1[1],res$ts2[1],res$ts3[1]))
    print(cbind(lagpar,'x',it))
  }
  saveRDS(RECBest_unemb_md,'MdRQA_Emb_RR.RDS')}

# Sine wave - JRQA

RECBest_emb_j=c()
t=seq(0,4+36/500,4/500)
sine=sin(t*2*pi)
for(lagpar in 0:5) {
  for(it in 1:100) {
    ts1 <- sine[1:(500+lagpar)]
    
    # Create a lagged version of ts1 with added noise
    ts2 <- scale(ts1[(lagpar+1):length(ts1)] + rnorm(500)*sd(ts1)/5)
    ts3 <- scale(ts1[(lagpar+1):length(ts1)] + rnorm(500)*sd(ts1)/5)
    
    # Merge the three time series
    data <- as.matrix(cbind(ts1[1:500],ts2,ts3))
    res=laggedJrqa(maxlag=lagpar+2, ts1=data, ts2=data, delay=1, embed=2, rescale=0,
                   radius=0.6, normalize=0, mindiagline=1, minvertline=1, tw=1, method='mdcrqa')
    
    # Find the best lag combination based on RR
    RECBest_emb_j=rbind(RECBest,cbind(lagpar,it,res$RR[1],res$ts1[1],res$ts2[1],res$ts3[1]))
    print(cbind(lagpar,'x',it))
  }
  saveRDS(RECBest_emb_j,'JRQA_RR.RDS')}
