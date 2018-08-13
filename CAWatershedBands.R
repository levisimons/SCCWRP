library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)

#This script looks at co-occurrence networks generated on a HUC-8 watershed scale from the SCCWRP
#California data set.  Each watershed is divided into land use bands, and then the core networks
#which persist across all land use bands by watershed are compared to infer the relationship between
#co-occurrence network structure and stress.
setwd("~/Desktop/SCCWRP")

#Read in site data containing biological counts, water chemistry, and land usage values.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#How many samples per group for generating a co-occurrence network?
groupNum=10
#Run through analysis on SCCWRP archive on a watershed-level scale.
#Select watersheds with a large enough set of samples for analysis.
watersheds <- subset(as.data.frame(table(SCCWRP$Watershed)),Freq>=3*groupNum)
colnames(watersheds) <- c("Watershed","Samples")

#Watersheds with a large enough set of samples to generate co-occurrence networks
#for each land use band.
LargeWS <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(LargeWS) <- c("Watershed","Samples")

#Set land use thresholds.
PD <- 5 #Partially degraded land use threshold.
HD <- 15 #Highly degraded land use threshold.
#Subset data by watershed, and then by land use bands.
for(WS in watersheds$Watershed){
  #Get samples per watershed.
  GISBioDataSubset <- subset(GISBioData,Watershed==WS)
  j=0
  for(i in 1:3){
    if(i==1){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K<=PD)
      if(length(unique(GISBioDataSubsetLU$UniqueID))>=groupNum){
        j=j+1
      }
    }
    if(i==2){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K > PD & LU_2000_5K<=HD)
      if(length(unique(GISBioDataSubsetLU$UniqueID))>=groupNum){
        j=j+1
      }
    }
    if(i==3){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K > HD)
      if(length(unique(GISBioDataSubsetLU$UniqueID))>=groupNum){
        j=j+1
      }
    }
    #Get watersheds with a minimum number of samples per watershed and land use band.
    if(j==3){
      LargeWS <- rbind(LargeWS,subset(watersheds,watersheds$Watershed==WS))
      }
  }
}

for(WS in LargeWS$Watershed){
  GISBioDataSubset <- subset(GISBioData,Watershed==WS)
  for(i in 1:3){
    if(i==1){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K<=PD)
    }
    if(i==2){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K > PD & LU_2000_5K<=HD)
    }
    if(i==3){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K > HD)
    }
    selected <- GISBioDataSubsetLU
    eLSAInput <- as.data.frame(unique(selected$FinalID))
    colnames(eLSAInput)<-c("FinalID")
    eLSAInput <- as.data.frame(eLSAInput[order(as.character(eLSAInput$FinalID)),])
    colnames(eLSAInput)<-c("FinalID")
    
    #Add the relative taxa abundances by column to a new dataframe.
    #The rows are the unique taxa in a given subset of data.
    selected <- selected[order(selected$Year,selected$UniqueID,selected$FinalID),]
    for(ID in unique(selected$UniqueID)){
      tmp <- filter(selected, UniqueID == ID)[,c("FinalID","Measurement","UniqueID")]
      tmp <- as.data.frame(tmp[order(tmp$FinalID),])
      tmp <- tmp[-c(3)]
      colnames(tmp)<-c("FinalID",paste("Measurement",ID,sep=" "))
      eLSAInput <- join(eLSAInput,tmp,by="FinalID")
      eLSAInput$FinalID=as.character(eLSAInput$FinalID)
      eLSAInput <- eLSAInput %>% group_by(FinalID) %>% summarise_if(is.numeric,mean,na.rm=TRUE)
    }
    eLSAInput[is.na(eLSAInput)] <- "NA"
    
    #Determine the number of time points in the eLSA input file.
    spotNum = length(unique(selected$Year))
    #Determine the number of replicates per time point in the eLSA input file.
    #In order to ensure a uniform number of replicates per year this needs to
    #be the maximum number of replicates for all of the years available.
    repMax = 0
    for(year in unique(selected$Year)){
      tmp <- filter(selected, Year == year)[,c("UniqueID","Year")]
      repNum = length(unique(tmp$UniqueID))
      if(repNum >= repMax){repMax = repNum}
      #print (paste(repMax,repNum,year,sep=" "))
    }
    repNum = repMax
    
    #Now insert the replicates with actual data in between the "NA" dummy columns
    #which ensure that the final eLSA input file has an even number of replicates
    #per year regardless of the variations in the actual number of sites (replicates)
    #sampled per year.
    eLSAtmp <- eLSAInput[,1]
    j=1
    k=1
    nulCol <- data.frame(matrix(ncol=repNum*spotNum,nrow=length(unique(selected$FinalID))))
    nulCol[,1] <- eLSAInput[,1]
    for(year in unique(selected$Year)){
      tmp <- filter(selected, Year == year)
      rep = length(unique(tmp$UniqueID))
      for(i in 1:repNum){
        if(i <= rep){
          repLabel = paste(year,"DoneRep",i,sep="")
          j=j+1
          k=k+1
          eLSAtmp[,k] <- eLSAInput[,j]
        }
        else{
          repLabel = as.character(paste(year,"Rep",i,sep=""))
          k=k+1
          eLSAtmp[,k] <- "NA"
          #print(paste(k,repLabel,sep=" "))
        }
      }
    }
    
    eLSAInput <- eLSAtmp
    #Designate a unique filename.
    #S is the number of spots, or years represented in the subsample group.
    #R is the number of replicates per year.  Many of the years will have null replicates, but a uniform number is needed for eLSA.
    #M is the mean LU_2000_5K score per subsample group.
    sampleNum <- length(unique(selected$UniqueID))
    meanLU = mean(na.omit(selected$LU_2000_5K))
    filename <- paste("Watershed",WS,"NSamples",sampleNum,"S",spotNum,"R",repNum,"M",meanLU,sep="")
    #Output file for use in eLSA.
    write.table(eLSAInput,paste(filename,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    eLSACommand = paste("lsa_compute ",filename,".txt ","-r ",repNum," -s ",spotNum," ",filename,"Network.txt;",sep="")
    print(eLSACommand)
  }
}
