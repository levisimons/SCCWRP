library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)
library(hierDiversity)
library(zetadiv)


#This script focuses on zeta diversity patterns across the SCCWRP California streams data set.
setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage
#values.  If this file is not yet generated then proceed with the following commands
#to generate it in the first place.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Remove duplicate measures.
GISBioData <- GISBioData[!duplicated(GISBioData[,c("UniqueID","FinalID","Count")]),]
#Order data by LU_2000_5K.
GISBioData <- arrange(GISBioData,LU_2000_5K)
#Add taxa counts by sample.
tmp <- data.frame(table(GISBioData$UniqueID))
colnames(tmp) <- c("UniqueID","nTaxa")
GISBioData <- join(GISBioData,tmp,by=c("UniqueID"))
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Get latitude and longitude bounds of data set.
latMin <- min(SCCWRP$Latitude)
latMax <- max(SCCWRP$Latitude)
lonMin <- min(SCCWRP$Longitude)
lonMax <- max(SCCWRP$Longitude)
#Merge in sample altitude.
GISBioData <- join(GISBioData,SCCWRP[,c("UniqueID","altitude")],by=c("UniqueID"))
#Get samples per watershed.
watersheds <- as.data.frame((table(SCCWRP$Watershed)))
colnames(watersheds) <- c("Watershed","Samples")
#Get the samples per watershed for watersheds with at least a certain number of samples.
LargeWatersheds <- subset(watersheds,Samples>=40)
#Taxa frequency table.
taxaFreq <- as.data.frame(table(GISBioData$FinalID))
colnames(taxaFreq) <- c("FinalID","Freq")
#Find the total number of taxa in the full data set.
taxaMax <- length(unique(GISBioData$FinalID))

#Get number of unique LU_2000_5K values.
sitesNum <- length(unique(GISBioData$UniqueID))
#Enter number of divisions for subsampling.
divisionNum = 60
#Obtain subsampling number.
sampleNum <- as.integer(sitesNum/divisionNum)
uniqueSamples <- as.data.frame(unique(GISBioData$UniqueID))
colnames(uniqueSamples) <- c("UniqueID")

zetaAnalysis <- data.frame()
for(i in 1:divisionNum){
  lowNum=(i-1)*sampleNum+1
  highNum=i*sampleNum
  GISBioData <- arrange(GISBioData,LU_2000_5K)
  uniqueSampleSubset <- as.data.frame(uniqueSamples[lowNum:highNum,1])
  colnames(uniqueSampleSubset) <- c("UniqueID")
  GISBioDataSubset <- GISBioData[GISBioData$UniqueID %in% as.vector(uniqueSampleSubset$UniqueID),]
  #Determine the average LU_2000_5K per subsample of sites.
  meanLU_2000_5K = mean(na.omit(GISBioDataSubset$LU_2000_5K))
  print(paste(lowNum,highNum,meanLU_2000_5K))
  #Initialize a data frame where the rows are all of the unique measurements for a given
  #subset of the data.
  #Order the data frame by measurement name.
  selected <- arrange(GISBioDataSubset,Year,UniqueID)
  
  #Generating a presence/absence matrix for California SCCWRP data.
  eLSAInput <- as.data.frame(unique(selected$FinalID))
  colnames(eLSAInput)<-c("FinalID")
  eLSAInput <- as.data.frame(eLSAInput[order(as.character(eLSAInput$FinalID)),])
  colnames(eLSAInput)<-c("FinalID")
  taxa <- eLSAInput
  #Add the relative taxa abundances by column to a new dataframe.
  #The rows are the unique taxa in a given subset of data.
  selected <- selected[order(selected$Year,selected$UniqueID,selected$FinalID),]
  for(ID in unique(selected$UniqueID)){
    tmp <- filter(selected, UniqueID == ID)[,c("FinalID","Measurement","UniqueID")]
    tmp <- as.data.frame(tmp[order(tmp$FinalID),])
    tmp <- tmp[-c(3)]
    colnames(tmp)<-c("FinalID",ID)
    tmp <- tmp %>% group_by(FinalID) %>% summarise_if(is.numeric,mean,na.rm=TRUE)
    tmp <- join(tmp,taxa,type="full",by=c("FinalID"))
    tmp <- as.data.frame(tmp[order(tmp$FinalID),])
    eLSAInput <- cbind(eLSAInput,tmp)
    eLSAInput <- eLSAInput[,!duplicated(colnames(eLSAInput))]
  }
  
  #Generate a presence/absence dataframe for zeta diversity analysis.
  #Rows for samples, columns for taxa IDs.
  eLSAInput[is.na(eLSAInput)] <- 0
  eLSANames <- eLSAInput$FinalID
  data.SCCWRP <- as.data.frame(t(eLSAInput[,-c(1)]))
  colnames(data.SCCWRP) <- eLSANames
  data.SCCWRP[data.SCCWRP > 0] <- 1
  
  #Computes zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
  #using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
  zetaDecay <- Zeta.decline.mc(data.SCCWRP,xy=NULL,orders=1:10,sam=1000)
  
  dat <- data.frame()
  dat[1,1] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
  dat[1,2] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
  dat[1,3] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
  dat[1,4] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
  dat[1,5] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
  dat[1,6] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
  
  zetaAnalysis <- rbind(zetaAnalysis,dat)
  print(dat)
}
colnames(zetaAnalysis) <- c("ZetaExponentialIntercept","ZetaExponentialExponent","ZetaExponentialAIC","ZetaPLIntercept","ZetaPLExponent","ZetaPLAIC")
write.table(zetaAnalysis,"LU_2000_5KSiteSweepCAZeta.txt",quote=FALSE,sep="\t",row.names = FALSE)

networkAnalysis <- read.table("LU_2000_5KSiteSweepCA.txt",header=TRUE)
networkAnalysis <- cbind(networkAnalysis,zetaAnalysis)
