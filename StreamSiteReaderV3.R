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
#Read in site data containing biological counts, CSCI, and land usage values.

#Read in aggregated biological data organized by taxonomic counts.
bioDataRaw <- read.table("CombinedTaxonomy.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
bioData <- filter(bioDataRaw, Replicate==1)
#Filter out entries without a FinalID
bioData <- filter(bioData,FinalID!="")

#Read in older merged data set to get unique algal and invertebrate taxa IDs.
oldSet <- read.table("GISBioDataSmallSet.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)

#Separate algae counts.
algae <- filter(oldSet,MeasurementType=="Soft-bodied algal relative abundance" | MeasurementType=="Benthic algal relative abundance")
algaeID <- unique(algae$FinalID)
algaeData <- subset(bioData,bioData$FinalID %in% algaeID)
#Separate out benthic algae counts.
benthicAlgaeData <- algaeData[!is.na(algaeData$BAResult),]
benthicAlgaeData <- benthicAlgaeData[,c("StationCode","SampleDate","FinalID","BAResult")]
#Separate out soft-bodied algae counts.
softAlgaeData <- algaeData[!is.na(algaeData$Result),]
softAlgaeData <- softAlgaeData[,c("StationCode","SampleDate","FinalID","Result")]

#Separate benthic macroinvertebrate counts.
BMI <- filter(oldSet,MeasurementType=="Invertebrate relative abundance" | MeasurementType=="Invertebrate relative abundances")
BMIID <- unique (BMI$FinalID)
BMIData <- subset(bioData,bioData$FinalID %in% BMIID)
BMIData <- BMIData[,c("StationCode","SampleDate","FinalID","BAResult")]
tmp1<-data.table(BMIData)
#Sum insect counts for matching IDs, but with different life stages.
BMIData<-tmp1[,.(BAResult=sum(BAResult)),by=.(StationCode,SampleDate,FinalID)]

#Determine the benthic algae totals count column and make it a temporary dataframe.
tmp1<-data.table(benthicAlgaeData)
tmp1<-tmp1[,.(BAResult=sum(BAResult)),by=.(StationCode,SampleDate)]
colnames(tmp1) <- c("StationCode","SampleDate","ActualOrganismCount")
#Join the total benthic algal count by sample to each unique sample.
benthicAlgaeData <- join(benthicAlgaeData,tmp1,by=c("StationCode","SampleDate"))
#Calculate the relative abundance of benthic algae data.
benthicAlgaeData$Measurement <- with(benthicAlgaeData,BAResult/ActualOrganismCount)
#Add measurement type label.
benthicAlgaeData$MeasurementType <- "benthic alage relative abundance"
#Add unique ID per sample.
benthicAlgaeData$UniqueID <- paste(benthicAlgaeData$StationCode,benthicAlgaeData$SampleDate)
#Subset key columns.
benthicAlgaeData <- benthicAlgaeData[,c("StationCode","SampleDate","UniqueID","FinalID","MeasurementType","Measurement")]

#Determine the soft-bodied algae totals count column and make it a temporary dataframe.
tmp1<-data.table(softAlgaeData)
tmp1<-tmp1[,.(Result=sum(Result)),by=.(StationCode,SampleDate)]
colnames(tmp1) <- c("StationCode","SampleDate","ActualOrganismCount")
#Join the total benthic algal count by sample to each unique sample.
softAlgaeData <- join(softAlgaeData,tmp1,by=c("StationCode","SampleDate"))
#Calculate the relative abundance of benthic algae data.
softAlgaeData$Measurement <- with(softAlgaeData,Result/ActualOrganismCount)
#Add measurement type label.
softAlgaeData$MeasurementType <- "soft-bodied alage relative abundance"
#Add unique ID per sample.
softAlgaeData$UniqueID <- paste(softAlgaeData$StationCode,softAlgaeData$SampleDate)
#Subset key columns.
softAlgaeData <- softAlgaeData[,c("StationCode","SampleDate","UniqueID","FinalID","MeasurementType","Measurement")]

#Determine the benthic algae totals count column and make it a temporary dataframe.
tmp1<-data.table(BMIData)
tmp1<-tmp1[,.(BAResult=sum(BAResult)),by=.(StationCode,SampleDate)]
colnames(tmp1) <- c("StationCode","SampleDate","ActualOrganismCount")
#Join the total benthic macroinvertebrate count by sample to each unique sample.
BMIData <- join(BMIData,tmp1,by=c("StationCode","SampleDate"))
#Calculate the relative abundance of benthic algae data.
BMIData$Measurement <- with(BMIData,BAResult/ActualOrganismCount)
#Add measurement type label.
BMIData$MeasurementType <- "benthic macroinvertebrate relative abundance"
#Add unique ID per sample.
BMIData$UniqueID <- paste(BMIData$StationCode,BMIData$SampleDate)
#Subset key columns.
BMIData <- BMIData[,c("StationCode","SampleDate","UniqueID","FinalID","MeasurementType","Measurement")]

#Bind the data frames back together as a unified biological set.
bioData <- rbind(benthicAlgaeData,softAlgaeData,BMIData)
bioData <- bioData[order(bioData$UniqueID),]

#Read in CSCI and land usage data by sample.
GISData <- read.csv("CSCI.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)

#Join CSCI data set with the biological data set.
GISBioData <- join(bioData,GISData,by=c("StationCode","SampleDate"))

#Add in aggregated land coverage at a 5km watershed buffer.
GISBioData$LU_2000_5K <- with(GISBioData,Ag_2000_5K+CODE_21_2000_5K+URBAN_2000_5K)

#Write out merged data set to read back in the future as opposed to 
#generating it each time from component data files.
write.csv(GISBioData,file="CAGISBioData.csv",row.names=FALSE)
