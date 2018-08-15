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
#Read in site data containing biological counts, water chemistry, and land usage values.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Remove duplicate measures.
GISBioData <- GISBioData[!duplicated(GISBioData[,c("UniqueID","FinalID","Count")]),]
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#How many samples per watershed?
groupNum=20
#Run through analysis on SCCWRP archive on a watershed-level scale.
#Select watersheds with a large enough set of samples for analysis.
watersheds <- subset(as.data.frame(table(SCCWRP$Watershed)),Freq>=groupNum)
colnames(watersheds) <- c("Watershed","Samples")
#Get samples only found in more heavily sampled watersheds.
GISBioDataLargeWS <- subset(GISBioData,Watershed %in% watersheds$Watershed)

#Get the frequency of pairs of taxa showing up by watershed.
#For example, a pair showing up three times in one watershed, and eight times in a second watershed,
#will be counted as having shown up in two unique watersheds.
CAMatch <- data.frame()
for(WS in watersheds$Watershed){
  WSSubset <- subset(GISBioData,Watershed==WS)
  WSMatch <- data.frame()
  for(ID in unique(WSSubset$UniqueID)){
    WSSample <- subset(WSSubset,UniqueID==ID)
    if(length(unique(WSSample$FinalID))>2){
      taxaMatch <- as.data.frame(t(combn(unique(WSSample$FinalID),2)))
      WSMatch <- rbind(WSMatch,taxaMatch)
      WSMatch <- WSMatch[!duplicated(WSMatch[,c("V1","V2")]),]
      print(paste(WS,ID,length(unique(WSSample$FinalID))))
    }
  }
  WSName <- data.frame(matrix(nrow=nrow(WSMatch),ncol=1))
  colnames(WSName) <- c("Watershed")
  WSName$Watershed <- WS
  WSMatch <- cbind(WSMatch,WSName)
  CAMatch <- rbind(CAMatch,WSMatch)
}
CAMatch <- ddply(CAMatch, .(CAMatch$V1,CAMatch$V2),nrow)
colnames(CAMatch) <- c("V1","V2","NumWS")
write.table(CAMatch,"AssociationsPerWaterhsed.txt",quote=FALSE,sep="\t",row.names = FALSE)
