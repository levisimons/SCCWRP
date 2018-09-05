library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)
library(MASS)
library(zetadiv)
library(magrittr)
library(stats)
library(CINNA)
library(fitdistrplus)

#How does endemism change with land use given the sampling year and watershed?
setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage values.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Remove duplicate measures.
GISBioData <- GISBioData[!duplicated(GISBioData[,c("UniqueID","FinalID","Count")]),]
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Select watersheds with a large enough set of samples for analysis.
watersheds <- as.data.frame(table(SCCWRP$Watershed))
colnames(watersheds) <- c("Watershed","Samples")
GISBioData <- join(GISBioData,watersheds,by=c("Watershed"))
#Add the number of genera per sample
taxaBySample <- count(GISBioData,UniqueID)
colnames(taxaBySample) <- c("UniqueID","nTaxa")
GISBioData <- join(GISBioData,taxaBySample,by=c("UniqueID"))

#Randomize the the location of taxa present by watershed to see if trends in endemism are non-random
#GISBioDataRand <- GISBioData %>% group_by(Year) %>% mutate(FinalID=sample(FinalID))
GISBioDataRand <- GISBioData %>% mutate(Watershed=sample(Watershed))

selected <- GISBioDataRand

LUquantile <- quantile(GISBioData$LU_2000_5K,probs=seq(0,1,0.1))#To get land use quantiles.

#Find genera unique to a particular watershed for each year. 
endemism <-  data.frame()
for(year in unique(selected$Year)){
  annualSubset <- subset(selected, Year==year)
  for(watershed in unique(annualSubset$Watershed)){
    for(i in 1:length(LUquantile)){
      if(i>1){
        LULow <- LUquantile[i-1]
        LUHigh <- LUquantile[i]
        MidLU <- 0.5*(LULow+LUHigh)
        localSubset <- subset(annualSubset,Watershed==watershed & LU_2000_5K >= LULow & LU_2000_5K <= LUHigh)
        localTaxa <- unique(localSubset$FinalID) #Taxa found in a set of samples defined by sampling year, target watershed, and land use.
        externalSubset <- subset(annualSubset,Watershed!=watershed & LU_2000_5K >= LULow & LU_2000_5K <= LUHigh)
        externalTaxa <- unique(externalSubset$FinalID) #Taxa found in a set of samples defined by sampling year and land use, and which are external to a target watershed.
        endemicTaxa <- setdiff(localTaxa,externalTaxa) #Taxa unique to a set of samples defined by sampling year, target watershed, and land use.
        if(length(localTaxa) >0 & length(endemicTaxa) > 0){
          endemicRatio <- length(endemicTaxa)/length(localTaxa)
        } else{
          endemicRatio <- NA
        }
        row <- t(as.data.frame(c(watershed,year,LULow,LUHigh,MidLU,length(endemicTaxa),length(unique(localSubset$UniqueID)),length(unique(externalSubset$UniqueID)),endemicRatio)))
        print(paste(watershed,year,LULow,LUHigh,MidLU,length(endemicTaxa),length(unique(localSubset$UniqueID)),length(unique(externalSubset$UniqueID)),endemicRatio))
        endemism <- rbind(row,endemism)
      }
    }
  }
}
colnames(endemism) <- c("Watershed","Year","LULow","LUHigh","MidLU","EndemicTaxa","LocalSamples","ExternalSamples","endemicRatio")
rownames(endemism) <- 1:nrow(endemism)
endemism$Year <- as.numeric(as.character(endemism$Year))
endemism$LULow <- as.numeric(as.character(endemism$LULow))
endemism$LUHigh <- as.numeric(as.character(endemism$LUHigh))
endemism$MidLU <- as.numeric(as.character(endemism$MidLU))
endemism$EndemicTaxa <- as.numeric(as.character(endemism$EndemicTaxa))
endemism$LocalSamples <- as.numeric(as.character(endemism$LocalSamples))
endemism$ExternalSamples <- as.numeric(as.character(endemism$ExternalSamples))
endemism$endemicRatio <- as.numeric(as.character(endemism$endemicRatio))

endemismFull <- endemism

endemism <- subset(endemism, LocalSamples!=0 & ExternalSamples!=0) #Remove rows without local and/or external samples.
endemism <- na.omit(endemism) #Remove trivial endemic ratio entries

#Plotting trends in endemism
plot(endemism$MidLU,endemism$endemicRatio)
lines(smooth.spline(endemism$MidLU,endemism$endemicRatio,df=3),col="red")
cor.test(endemism$MidLU,endemism$endemicRatio,method="spearman")
