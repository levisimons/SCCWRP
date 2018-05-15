library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)

setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage
#values.  If this file is not yet generated then proceed with the following commands
#to generate it in the first place.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Ensure that all sites have a LU_2000_5K value.
GISBioData <- subset(GISBioData, LU_2000_5K != "NA")
#Order data by LU_2000_5K.
GISBioData <- arrange(GISBioData,LU_2000_5K)

#Get number of unique LU_2000_5K values.
sitesNum <- length(unique(GISBioData$UniqueID))
#Enter number of divisions for subsampling.
divisionNum = 60
#Obtain subsampling number.
sampleNum <- as.integer(sitesNum/divisionNum)
uniqueSamples <- as.data.frame(unique(GISBioData$UniqueID))
colnames(uniqueSamples) <- c("UniqueID")

for(i in 1:divisionNum){
  lowNum=(i-1)*sampleNum+1
  highNum=i*sampleNum
  GISBioData <- arrange(GISBioData,LU_2000_5K)
  uniqueSampleSubset <- as.data.frame(uniqueSamples[lowNum:highNum,1])
  colnames(uniqueSampleSubset) <- c("UniqueID")
  GISBioDataSubset <- GISBioData[GISBioData$UniqueID %in% as.vector(uniqueSampleSubset$UniqueID),]
  #Determine the average LU_2000_5K per subsample of sites.
  meanLU_2000_5K = mean(na.omit(GISBioDataSubset$LU_2000_5K))
  #Initialize a data frame where the rows are all of the unique measurements for a given
  #subset of the data.
  #Order the data frame by measurement name.
  selected <- arrange(GISBioDataSubset,Year,UniqueID)
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
    #print(ID)
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
  #N is the number of samples in the subsample group.
  #S is the number of spots, or years represented in the subsample group.
  #R is the number of replicates per year.  Many of the years will have null replicates, but a uniform number is needed for eLSA.
  #M is the mean LU_2000_5K score per subsample group.
  filename = paste("LUSweepCA",sampleNum,"Samples",lowNum,"to",highNum,"S",spotNum,"R",repNum,"M",meanLU_2000_5K,sep="")
  
  #Output file for use in eLSA.
  #write.table(eLSAInput,paste(filename,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  eLSACommand = paste("lsa_compute ",filename,".txt ","-r ",repNum," -s ",spotNum," ",filename,"Network.txt;",sep="")
  print(eLSACommand)
}
