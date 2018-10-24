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

#Create co-occurrence networks by watershed and land use decile in the SCCWRP data set.
setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage values.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Remove duplicate measures.
GISBioData <- GISBioData[!duplicated(GISBioData[,c("UniqueID","FinalID","Count")]),]
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Read in functional feeding group for each taxon.
#Abbreviations used in denoting functional feeding groups are as follows ( http://www.safit.org/Docs/CABW_std_taxonomic_effort.pdf ):
#P= predator MH= macrophyte herbivore OM= omnivore
#PA= parasite PH= piercer herbivore XY= xylophage (wood eater)
#CG= collector-gatherer SC= scraper
#CF= collector filterer SH= shredder 
FFG <- read.table("metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
# Filter data so only known functional feeding groups are kept.
FFG <- subset(FFG, FunctionalFeedingGroup != "")
# Generate functional feeding group data frame.
FFG <- FFG[,c("FinalID","LifeStageCode","FunctionalFeedingGroup")]
FFG <- subset(FFG,LifeStageCode=="L" | LifeStageCode=="X" | FinalID=="Hydrophilidae" | FinalID=="Hydraenidae")
#Merge in functional feeding groups into sample data.
GISBioData <- join(GISBioData,FFG[,c("FinalID","FunctionalFeedingGroup")],by=c("FinalID"))
FFGCounts <- na.omit(as.data.frame(unique(GISBioData$FunctionalFeedingGroup)))
colnames(FFGCounts) <- c("FunctionalFeedingGroups")
#How many samples per watershed?
groupNum=25
#Select watersheds with a large enough set of samples for analysis.
watersheds <- as.data.frame(table(SCCWRP$Watershed))
colnames(watersheds) <- c("Watershed","Samples")
GISBioData <- join(GISBioData,watersheds,by=c("Watershed"))
#Get samples only found in more heavily sampled watersheds.
GISBioDataLargeWS <- subset(GISBioData,Samples>=groupNum)

LUquantile <- quantile(GISBioDataLargeWS$LU_2000_5K,probs=seq(0,1,0.2))#To get land use quantiles.

for(watershed in unique(GISBioDataLargeWS$Watershed)){
  for(i in 1:length(LUquantile)){
    if(i>1){
      LULow <- LUquantile[i-1]
      LUHigh <- LUquantile[i]
      MidLU <- 0.5*(LULow+LUHigh)
      localSubset <- subset(GISBioDataLargeWS,Watershed==watershed & LU_2000_5K >= LULow & LU_2000_5K <= LUHigh)
      if(length(unique(localSubset$UniqueID)) >= groupNum){
        localTaxa <- unique(localSubset$FinalID) #Taxa found in a set of samples defined by target watershed and land use.
        print(paste(watershed,LULow,LUHigh,length(localTaxa)))
      }
    }
  }
}
