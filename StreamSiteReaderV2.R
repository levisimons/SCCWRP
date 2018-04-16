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
#values.  Generate a merged data set.

#GISBiochemData <- read.table("GISBiochemData.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)

#If you need to aggregate site data please proceed here.
#Read in algae data from SMC sites.
algaeDataSMCRaw <- read.table("AlgaeTax_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
algaeDataSMC <- filter(algaeDataSMCRaw, Replicate==1)
#Subset columns of interest for the SMC sites.
algaeDataSMC <- algaeDataSMCRaw[,c(2,3,43,44,40)]
#Determine the algal totals count column and make it a temporary dataframe.
tmp1 <- as.data.frame(xtabs(BAResult ~ StationCode,algaeDataSMC))
colnames(tmp1) <- c("StationCode","ActualOrganismCount")
#Determine the algal volumes column and make it a temporary dataframe.
tmp2 <- as.data.frame(xtabs(Result ~ StationCode,algaeDataSMC))
colnames(tmp2) <- c("StationCode","ActualOrganismVolume")
#Add algal totals count column to algae dataframe.
algaeDataSMC <- merge(algaeDataSMC,tmp1,"StationCode")
algaeDataSMC <- merge(algaeDataSMC,tmp2,"StationCode")
#Temporarily split data frame into soft-bodied and benthic algae sets.
#Determine each algal sets relative abundances and then merge data frames back.
tmp1 <- filter(algaeDataSMC,BAResult!="NA")
#Calculate the relative abundance of benthic algae data.
tmp1$Measurement <- with(tmp1,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
tmp1$MeasurementType <- with(tmp1,"Benthic algal relative abundance")
tmp2 <- filter(algaeDataSMC, Result!="NA")
#Calculate the relative abundance of soft-bodied algae data.
tmp2$Measurement <- with(tmp2,Result/ActualOrganismVolume)
#Add organism type for later use in merged data sets.
tmp2$MeasurementType <- with(tmp2,"Soft-bodied algal relative abundance")
algaeDataSMC <- rbind(tmp1,tmp2)
#Force a uniform date format
algaeDataSMC$SampleDate <- mdy(algaeDataSMC$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
algaeDataSMC$UniqueID <- with(algaeDataSMC,paste(algaeDataSMC$StationCode,"SMC",algaeDataSMC$SampleDate))
#Find sampling year.
algaeDataSMC$Year <- year(algaeDataSMC$SampleDate)

#Read in algae data from CEDEN sites.
algaeDataCEDENRaw <- read.table("AlgaeTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
algaeDataCEDEN <- filter(algaeDataCEDENRaw, CollectionReplicate==1)
#Subset columns of interest for the CEDEN sites.
algaeDataCEDEN <- algaeDataCEDENRaw[,c(6,11,36,35,26)]
names(algaeDataCEDEN)[names(algaeDataCEDEN)=="Counts"]<-"Result"
#Determine the algal totals count column and make it a temporary dataframe.
tmp1 <- as.data.frame(xtabs(BAResult ~ StationCode,algaeDataCEDEN))
colnames(tmp1) <- c("StationCode","ActualOrganismCount")
#Determine the algal volumes column and make it a temporary dataframe.
tmp2 <- as.data.frame(xtabs(Result ~ StationCode,algaeDataCEDEN))
colnames(tmp2) <- c("StationCode","ActualOrganismVolume")
#Add algal totals count column to algae dataframe.
algaeDataCEDEN <- merge(algaeDataCEDEN,tmp1,"StationCode")
algaeDataCEDEN <- merge(algaeDataCEDEN,tmp2,"StationCode")
#Temporarily split data frame into soft-bodied and benthic algae sets.
#Determine each algal sets relative abundances and then merge data frames back.
tmp1 <- filter(algaeDataCEDEN,BAResult!="NA")
#Calculate the relative abundance of benthic algae data.
tmp1$Measurement <- with(tmp1,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
tmp1$MeasurementType <- with(tmp1,"Benthic algal relative abundance")
tmp2 <- filter(algaeDataCEDEN, Result!="NA")
#Calculate the relative abundance of soft-bodied algae data.
tmp2$Measurement <- with(tmp2,Result/ActualOrganismVolume)
#Add organism type for later use in merged data sets.
tmp2$MeasurementType <- with(tmp2,"Soft-bodied algal relative abundance")
algaeDataCEDEN <- rbind(tmp1,tmp2)
#Force a uniform date format
algaeDataCEDEN$SampleDate <- ymd(algaeDataCEDEN$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
algaeDataCEDEN$UniqueID <- with(algaeDataCEDEN,paste(algaeDataCEDEN$StationCode,"CEDEN",algaeDataCEDEN$SampleDate))
#Find sampling year.
algaeDataCEDEN$Year <- year(algaeDataCEDEN$SampleDate)

#The SWAMP data file is in a somewhat irregular format and this is accounted for
#when being read in.
algaeDataSWAMPRaw <- read.table("AlgaeTaxonomy_dnaSamples_SWAMP.csv", fill=TRUE,header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
algaeDataSWAMP <- filter(algaeDataSWAMPRaw, Replicate==1)
algaeDataSWAMP <- algaeDataSWAMPRaw[,c(6,8,97,99,90)]
#Determine the algal totals count column and make it a temporary dataframe.
tmp1 <- as.data.frame(xtabs(BAResult ~ StationCode,algaeDataSWAMP))
colnames(tmp1) <- c("StationCode","ActualOrganismCount")
#Determine the algal volumes column and make it a temporary dataframe.
tmp2 <- as.data.frame(xtabs(Result ~ StationCode,algaeDataSWAMP))
colnames(tmp2) <- c("StationCode","ActualOrganismVolume")
#Add algal totals count column to algae dataframe.
algaeDataSWAMP <- merge(algaeDataSWAMP,tmp1,"StationCode")
algaeDataSWAMP <- merge(algaeDataSWAMP,tmp2,"StationCode")
#Temporarily split data frame into soft-bodied and benthic algae sets.
#Determine each algal sets relative abundances and then merge data frames back.
tmp1 <- filter(algaeDataSWAMP,BAResult!="NA")
#Calculate the relative abundance of benthic algae data.
tmp1$Measurement <- with(tmp1,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
tmp1$MeasurementType <- with(tmp1,"Benthic algal relative abundance")
tmp2 <- filter(algaeDataSWAMP, Result!="NA")
#Calculate the relative abundance of soft-bodied algae data.
tmp2$Measurement <- with(tmp2,Result/ActualOrganismVolume)
#Add organism type for later use in merged data sets.
tmp2$MeasurementType <- with(tmp2,"Soft-bodied algal relative abundance")
algaeDataSWAMP <- rbind(tmp1,tmp2)
#Force a uniform date format
algaeDataSWAMP$SampleDate <- mdy(algaeDataSWAMP$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
algaeDataSWAMP$UniqueID <- with(algaeDataSWAMP,paste(algaeDataSWAMP$StationCode,"SWAMP",algaeDataSWAMP$SampleDate))
#Find sampling year.
algaeDataSWAMP$Year <- year(algaeDataSWAMP$SampleDate)

#Create merged algae data set.
algaeData <- do.call("rbind",list(algaeDataSMC,algaeDataSWAMP,algaeDataCEDEN))
#Remove raw count data.
algaeData <- within(algaeData,rm("BAResult","Result","ActualOrganismCount","ActualOrganismVolume"))
#Reorder columns prior to merger.
algaeData <- algaeData[c("StationCode","SampleDate","FinalID","Measurement","MeasurementType","UniqueID","Year")]

#Read in insect data from SMC sites.
insectDataSMCRAW <- read.csv("BugTax_dnaSites_SMC.csv")
#Subset only replicate 1
insectDataSMC <- filter(insectDataSMCRAW, FieldReplicate==1)
#Subset columns of interest.
insectDataSMC <- insectDataSMCRAW[,c(2,3,9,6)]
#Sum insect counts for matching IDs, but with different life stages.
insectDataSMC <- insectDataSMC %>% group_by(StationCode,SampleDate,FinalID) %>% summarize(BAResult=sum(BAResult))
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ StationCode,insectDataSMC))
colnames(tmp) <- c("StationCode","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataSMC <- merge(insectDataSMC,tmp,"StationCode")
#Calculate the relative abundance of insect data.
insectDataSMC$Measurement <- with(insectDataSMC,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataSMC$MeasurementType <- with(insectDataSMC,"Invertebrate relative abundances")
#Force a uniform date format
insectDataSMC$SampleDate <- mdy(insectDataSMC$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
insectDataSMC$UniqueID <- with(insectDataSMC,paste(insectDataSMC$StationCode,"SMC",insectDataSMC$SampleDate))
#Find sampling year.
insectDataSMC$Year <- year(insectDataSMC$SampleDate)

#Read in insect data from CEDEN sites.
insectDataCEDENRAW <- read.table("BugTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
insectDataCEDEN <- filter(insectDataCEDENRAW, CollectionReplicate==1)
#Subset columns of interest.
insectDataCEDEN <- insectDataCEDENRAW[,c(6,11,36,26)]
#Sum insect counts for matching IDs, but with different life stages.
insectDataCEDEN <- insectDataCEDEN %>% group_by(StationCode,SampleDate,FinalID) %>% summarize(BAResult=sum(BAResult))
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ StationCode,insectDataCEDEN))
colnames(tmp) <- c("StationCode","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataCEDEN <- merge(insectDataCEDEN,tmp,"StationCode")
#Calculate the relative abundance of insect data.
insectDataCEDEN$Measurement <- with(insectDataCEDEN,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataCEDEN$MeasurementType <- with(insectDataCEDEN,"Invertebrate relative abundance")
#Force a uniform date format
insectDataCEDEN$SampleDate <- ymd(insectDataCEDEN$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
insectDataCEDEN$UniqueID <- with(insectDataCEDEN,paste(insectDataCEDEN$StationCode,"CEDEN",insectDataCEDEN$SampleDate))
#Find sampling year.
insectDataCEDEN$Year <- year(insectDataCEDEN$SampleDate)

#Read in insect data from SWAMP sites.
insectDataSWAMPRAW <- read.table("BugTaxonomy_dnaSamples_SWAMP.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Subset only replicate 1
insectDataSWAMP <- filter(insectDataSWAMPRAW, Replicate==1)
#Subset columns of interest.
insectDataSWAMP <- insectDataSWAMPRAW[,c(6,8,97,90)]
#Sum insect counts for matching IDs, but with different life stages.
insectDataSWAMP <- insectDataSWAMP %>% group_by(StationCode,SampleDate,FinalID) %>% summarize(BAResult=sum(BAResult))
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ StationCode,insectDataSWAMP))
colnames(tmp) <- c("StationCode","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataSWAMP <- merge(insectDataSWAMP,tmp,"StationCode")
#Calculate the relative abundance of insect data.
insectDataSWAMP$Measurement <- with(insectDataSWAMP,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataSWAMP$MeasurementType <- with(insectDataSWAMP,"Invertebrate relative abundance")
#Force a uniform date format
insectDataSWAMP$SampleDate <- mdy(insectDataSWAMP$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
insectDataSWAMP$UniqueID <- with(insectDataSWAMP,paste(insectDataSWAMP$StationCode,"SWAMP",insectDataSWAMP$SampleDate))
#Find sampling year.
insectDataSWAMP$Year <- year(insectDataSWAMP$SampleDate)

#Create merged insect data set.
insectData <- do.call("rbind",list(insectDataSMC,insectDataSWAMP,insectDataCEDEN))
#Remove raw count data.
insectData <- within(insectData,rm("BAResult","ActualOrganismCount"))
#Reorder columns prior to merger.
insectData <- insectData[c("StationCode","SampleDate","FinalID","Measurement","MeasurementType","UniqueID","Year")]

#Merge insect and algae data.
bioData <- do.call("rbind",list(insectData,algaeData))
bioData <- bioData[!duplicated(bioData),]
#Reorder columns post merger.
bioData <- bioData[c("StationCode","SampleDate","FinalID","Measurement","MeasurementType","UniqueID","Year")]
