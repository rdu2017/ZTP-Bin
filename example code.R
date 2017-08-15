
######## set up environment

#rm(list=ls(all=TRUE)) # remove all objects in the environment, if needed

setwd('.../') # set path to the working directory 

source('Param.est.R') # get the mle.lambda function into current environment

source('ZTP-Bin.R') # get the ZTPBin function into current environment


######## estimate lambda 

Simu_Ui <- read.table('./data/Simu_Ui.txt', header=T, stringsAsFactors = F) # read in the file containing Ui(s) to each COG
Simu_Ui<-Simu_Ui[-which(Simu_Ui$Ui==0),] # remove those COGs to which no reads being aligned

Simu_lambda_mle <- NULL
for (c in unique(as.character(Simu_Ui$COG))){
  loc <- which(Simu_Ui$COG==c)
  Ui <- Simu_Ui[loc,]$Ui
  Lambda <- round(mle.lambda(Ui,10)) # use mle.lambda function to esimtate lambda for each COG
  Simu_lambda_mle <- rbind(Simu_lambda_mle,
                           data.frame(COG=as.character(c),
                                    lambda=as.numeric(Lambda),
                                    stringsAsFactors = F))} 

Simu_lambda_mle <- Simu_lambda_mle[!Simu_lambda_mle$lambda==0,] # remove those COGs with lambda estimated 0


######## estimate p

Simu_blaBH <- read.table("./data/Simu_blastBH.txt", header=T) # read in the file containing alignment scores, and to which COG the read is trully associated 

Simu_blaBHc <- Simu_blaBH[which(Simu_blaBH$score > 66),] # retain reads with alignment score greater than the cutoff 

Zsum <- table(Simu_blaBHc$COG) # Zsum: the total number of reads assigned to a COG after score cutoff
Simu_Zsum <- data.frame(COG = as.character(names(Zsum)), Zsum = as.numeric(Zsum)) 

Simu_Usum <- NULL
for (c in unique(as.character(Simu_Ui$COG)))
{
  loc<-which(Simu_Ui$COG==c)
  Usum<-sum(Simu_Ui$Ui[loc])
  Simu_Usum<-rbind(Simu_Usum,data.frame(COG=as.character(c),Usum=as.numeric(Usum)))
} # Usum: the total number of reads fragmented from a COG

Simu_Usum_Zsum <- merge(Simu_Usum, Simu_Zsum, by.x="COG", by.y="COG")
Simu_Usum_Zsum <- Simu_Usum_Zsum[!Simu_Usum_Zsum$Zsum==0,]

Simu_p_mle <- data.frame(COG = Simu_Usum_Zsum$COG,
                         p = as.numeric(Simu_Usum_Zsum$Zsum/Simu_Usum_Zsum$Usum))


Simu_lambda_p <- merge(Simu_lambda_mle, Simu_p_mle, by='COG') # this dataset contains the MLEs of lambda and p


######## ZTP-Bin model fitting

HOT25_crtCount <- read.table('./data/HOT25_crtCount.txt', header=T) # read in a real dataset which has gone through Step 1&2

HOT25_relAbun <- ZTPBin(estParam = Simu_lambda_p, 
                        crtCount = HOT25_crtCount, 
                        nonCDS='Non-CDS') # output after ZTP-Bin model fitting 
