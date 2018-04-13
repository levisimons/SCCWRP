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
algaeDataSMC <- algaeDataSMCRaw[,c(1,3,43,40)]
#Change the header name for station ID.
names(algaeDataSMC)[names(algaeDataSMC)=="Sample Station ID"]<-"SampleStationID"
#Determine the algal totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,algaeDataSMC))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add algal totals count column to algae dataframe.
algaeDataSMC <- merge(algaeDataSMC,tmp,"SampleStationID")
#Calculate the relative abundance of algae data.
algaeDataSMC$Measurement <- with(algaeDataSMC,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
algaeDataSMC$MeasurementType <- with(algaeDataSMC,"Algal relative abundance")
#Force a uniform date format
algaeDataSMC$SampleDate <- mdy(algaeDataSMC$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
algaeDataSMC$UniqueID <- with(algaeDataSMC,paste(algaeDataSMC$SampleStationID,"SMC",algaeDataSMC$SampleDate))
#Find sampling year.
algaeDataSMC$Year <- year(algaeDataSMC$SampleDate)

#Read in algae data from CEDEN sites.
algaeDataCEDENRaw <- read.table("AlgaeTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
algaeDataCEDEN <- filter(algaeDataCEDENRaw, CollectionReplicate==1)
#Subset columns of interest for the CEDEN sites.
algaeDataCEDEN <- algaeDataCEDENRaw[,c(6,11,36,26)]
#Change names to uniforma schema.
names(algaeDataCEDEN)[names(algaeDataCEDEN)=="StationCode"]<-"SampleStationID"
#names(algaeDataCEDEN)[names(algaeDataCEDEN)=="Species"]<-"FinalID"
#Determine the algal totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,algaeDataCEDEN))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add algal totals count column to algae dataframe.
algaeDataCEDEN <- merge(algaeDataCEDEN,tmp,"SampleStationID")
#Calculate the relative abundance of algae data.
algaeDataCEDEN$Measurement <- with(algaeDataCEDEN,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
algaeDataCEDEN$MeasurementType <- with(algaeDataCEDEN,"Algal relative abundance")
#Force a uniform date format
algaeDataCEDEN$SampleDate <- ymd(algaeDataCEDEN$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
algaeDataCEDEN$UniqueID <- with(algaeDataCEDEN,paste(algaeDataCEDEN$SampleStationID,"CEDEN",algaeDataCEDEN$SampleDate))
#Find sampling year.
algaeDataCEDEN$Year <- year(algaeDataCEDEN$SampleDate)

#The SWAMP data file is in a somewhat irregular format and this is accounted for
#when being read in.
algaeDataSWAMPRaw <- read.table("AlgaeTaxonomy_dnaSamples_SWAMP.csv", fill=TRUE,header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
algaeDataSWAMP <- filter(algaeDataSWAMPRaw, Replicate==1)
algaeDataSWAMP <- algaeDataSWAMPRaw[,c(6,8,97,90)]
#Change names to uniforma schema.
names(algaeDataSWAMP)[names(algaeDataSWAMP)=="StationCode"]<-"SampleStationID"
#Determine the algal totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,algaeDataSWAMP))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add algal totals count column to algae dataframe.
algaeDataSWAMP <- merge(algaeDataSWAMP,tmp,"SampleStationID")
#Calculate the relative abundance of algae data.
algaeDataSWAMP$Measurement <- with(algaeDataSWAMP,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
algaeDataSWAMP$MeasurementType <- with(algaeDataSWAMP,"Algal relative abundance")
#Force a uniform date format
algaeDataSWAMP$SampleDate <- mdy(algaeDataSWAMP$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
algaeDataSWAMP$UniqueID <- with(algaeDataSWAMP,paste(algaeDataSWAMP$SampleStationID,"SWAMP",algaeDataSWAMP$SampleDate))
#Find sampling year.
algaeDataSWAMP$Year <- year(algaeDataSWAMP$SampleDate)

#Create merged algae data set.
algaeData <- do.call("rbind",list(algaeDataSMC,algaeDataSWAMP,algaeDataCEDEN))
algaeData <- na.omit(algaeData)

#Read in insect data from SMC sites.
insectDataSMCRAW <- read.csv("BugTax_dnaSites_SMC.csv")
#Subset only replicate 1
insectDataSMC <- filter(insectDataSMCRAW, FieldReplicate==1)
#Subset columns of interest.
insectDataSMC <- insectDataSMCRAW[,c(1,3,9,6)]
#Change names to uniforma schema.
names(insectDataSMC)[names(insectDataSMC)=="Sample.Station.ID"]<-"SampleStationID"
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,insectDataSMC))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataSMC <- merge(insectDataSMC,tmp,"SampleStationID")
#Calculate the relative abundance of insect data.
insectDataSMC$Measurement <- with(insectDataSMC,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataSMC$MeasurementType <- with(insectDataSMC,"Invertebrate relative abundances")
#Force a uniform date format
insectDataSMC$SampleDate <- mdy(insectDataSMC$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
insectDataSMC$UniqueID <- with(insectDataSMC,paste(insectDataSMC$SampleStationID,"SMC",insectDataSMC$SampleDate))
#Find sampling year.
insectDataSMC$Year <- year(insectDataSMC$SampleDate)

#Read in insect data from CEDEN sites.
insectDataCEDENRAW <- read.table("BugTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
insectDataCEDEN <- filter(insectDataCEDENRAW, CollectionReplicate==1)
#Subset columns of interest.
insectDataCEDEN <- insectDataCEDENRAW[,c(6,11,36,26)]
#Change names to uniforma schema.
names(insectDataCEDEN)[names(insectDataCEDEN)=="StationCode"]<-"SampleStationID"
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,insectDataCEDEN))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataCEDEN <- merge(insectDataCEDEN,tmp,"SampleStationID")
#Calculate the relative abundance of insect data.
insectDataCEDEN$Measurement <- with(insectDataCEDEN,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataCEDEN$MeasurementType <- with(insectDataCEDEN,"Invertebrate relative abundance")
#Force a uniform date format
insectDataCEDEN$SampleDate <- ymd(insectDataCEDEN$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
insectDataCEDEN$UniqueID <- with(insectDataCEDEN,paste(insectDataCEDEN$SampleStationID,"CEDEN",insectDataCEDEN$SampleDate))
#Find sampling year.
insectDataCEDEN$Year <- year(insectDataCEDEN$SampleDate)

#Read in insect data from SWAMP sites.
insectDataSWAMPRAW <- read.table("BugTaxonomy_dnaSamples_SWAMP.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Subset only replicate 1
insectDataSWAMP <- filter(insectDataSWAMPRAW, Replicate==1)
#Subset columns of interest.
insectDataSWAMP <- insectDataSWAMPRAW[,c(1,8,97,90)]
#Change names to uniforma schema.
names(insectDataSWAMP)[names(insectDataSWAMP)=="Sample Station ID"]<-"SampleStationID"
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,insectDataSWAMP))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataSWAMP <- merge(insectDataSWAMP,tmp,"SampleStationID")
#Calculate the relative abundance of insect data.
insectDataSWAMP$Measurement <- with(insectDataSWAMP,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataSWAMP$MeasurementType <- with(insectDataSWAMP,"Invertebrate relative abundance")
#Force a uniform date format
insectDataSWAMP$SampleDate <- mdy(insectDataSWAMP$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
insectDataSWAMP$UniqueID <- with(insectDataSWAMP,paste(insectDataSWAMP$SampleStationID,"SWAMP",insectDataSWAMP$SampleDate))
#Find sampling year.
insectDataSWAMP$Year <- year(insectDataSWAMP$SampleDate)

#Create merged insect data set.
insectData <- do.call("rbind",list(insectDataSMC,insectDataSWAMP,insectDataCEDEN))
insectData <- na.omit(insectData)
insectData$UniqueTaxa <- paste(insectData$UniqueID,insectData$FinalID)
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ UniqueTaxa,insectData))
colnames(tmp) <- c("UniqueTaxa","BAResult")
#Add insect totals count column to insect dataframe.
insectData <- merge(subset(insectData,select=-c(3)),tmp,"UniqueTaxa")
#Calculate the relative abundance of insect data.
insectData$Measurement <- with(insectData,BAResult/ActualOrganismCount)
insectData <- insectData[!duplicated(insectData),]
#Remove UniqueTaxa column
insectData <- subset(insectData,select=-c(1))
#Reorder columns prior to merger.
insectData <- insectData[c("SampleStationID","SampleDate","BAResult","FinalID","ActualOrganismCount","Measurement","MeasurementType","UniqueID","Year")]
#Determine the insect Shannon diversity and make it a temporary dataframe.

#Merge insect and algae data.
bioData <- do.call("rbind",list(insectData,algaeData))
#bioData <- bioData[,-c(3,5)]
bioData <- bioData[!duplicated(bioData),]
#Reorder columns post merger.
bioData <- bioData[c("SampleStationID","SampleDate","FinalID","Measurement","MeasurementType","UniqueID","Year","BAResult","ActualOrganismCount")]
