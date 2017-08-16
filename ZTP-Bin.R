
# Title: zero-truncated Poisson and Binomial (ZTP-Bin) modelling for length and conservation bias correction
# Date: 2017/07/31
# Description: The function estimates, for a metagenomic sample, the number of DNA coding sequences (CDS)
#              associated to a Cluster of Orthologous Groups of proteins (COG) after the sequence alignment 
#              output being corrected against the artificats and cross-annotations.
# The relative abundance of the CDSs given by the function is usually treated as a functional profile of the sample.

#rm(list=ls(all=TRUE)) # remove all objects in the environment, if needed

ZTPBin <- function (estParam, crtCount, nonCDS='Non-CDS') {
  
  # Estimates the number of CDS
  #
  # Args:
  #   estParm: a dataframe of estimated lambda and p for COG, containing three columns: COG, lambda, p
  #   crtCount: a dataframe of corrected counts to COGs from a query sample, containing two columns: COG, count
  #   nonCDS: the name for the non-coding sequence (Non-CDS) in the COG column from crtCount; if it does not exist, set nonCDS=NA 
  
  #
  # Returns:
  #   a dataframe that contains:
  #      COG: the COG group after artifacts being removed
  #      count: the sequence count to a COG after cross-annotation being corrected
  #      lambda: the estimated parameter representing the average count of reads fragmented from a COG
  #      p: the estimated parameter representing the proportion of the reads retained after alignment and the score cutoff
  #      NCDS: the predicted number of CDS associated to a COG
  #      relAbun: the relative abundance of a COG in the query sample
  
  if (!is.na(nonCDS))
    crtCount <- crtCount[!crtCount$COG==nonCDS, ] # if the query dataset has the Non-CDS row then remove it
  
  mlambda <- round(mean(estParam$lambda)) # mean of lambda among all the recorded COGs in estParam
  mp <- mean(estParam$p) # mean of p among all the recorded COGs in estParam

  out <- merge(crtCount, estParam, by.x="COG", by.y="COG", all.x=T) # query dataset merged with estimated lambda and p 
  
  if (length(which(is.na(out$lambda))) > 0)
    out$lambda[is.na(out$lambda)] <- mlambda # a missing lambda is replaced by the mean of lambda 

  if (length(which(is.na(out$p))) > 0)
    out$p[is.na(out$p)] <- mp # a missing p is replaced by the mean of p
  
  if (length(which(out$p==0)) > 0)
    out$p[out$p==0] <- mp # if a p equals zero, its' set to the mean of p

  out$NCDS<-round(out$count*(1-exp(-out$lambda))/(out$lambda*out$p)) # calculate the predicated number of CDS to each COG
  
  out$relAbun <- as.numeric(out$NCDS/sum(out$NCDS)) # calculate the realative abundance of a COG
  
  out <- out[order(-out$NCDS),] # set the order from most abundance to lease abundance
  
  return(out)
}
