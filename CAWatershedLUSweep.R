library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)

#This script focuses on generating co-occurrence networks on a HUC-8 watershed scale within the SCCWRP archive.
setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage
#values.  If this file is not yet generated then proceed with the following commands
#to generate it in the first place.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#How many samples per group for generating a co-occurrence network?
groupNum=30
#Run through analysis on SCCWRP archive on a watershed-level scale.
#Select watersheds with a large enough set of samples for analysis.
watersheds <- subset(as.data.frame(table(SCCWRP$Watershed)),Freq>=2*groupNum)
colnames(watersheds) <- c("Watershed","Samples")
for(WS in watersheds$Watershed){
  #Get samples per watershed.
  GISBioDataSubset <- subset(GISBioData,Watershed==WS)
  GISBioDataSubset <- arrange(GISBioDataSubset,LU_2000_5K)
  sampleNum <- length(unique(GISBioDataSubset$UniqueID))
  GISBioDataSubsetLow <- GISBioDataSubset[GISBioDataSubset$UniqueID %in% as.vector(unique(GISBioDataSubset$UniqueID)[1:groupNum]),]
  GISBioDataSubsetHigh <- GISBioDataSubset[GISBioDataSubset$UniqueID %in% as.vector(unique(GISBioDataSubset$UniqueID)[as.integer((sampleNum-groupNum)+1):as.integer(sampleNum)]),]
  for(i in 1:2){
    if(i==1){selected <- GISBioDataSubsetLow}
    if(i==2){selected <- GISBioDataSubsetHigh}
  }
  }
