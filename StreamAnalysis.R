library("plyr")
library(dplyr)
library("ggplot2")
library(phyloseq)
library("ape")
library("vegan")
library("microbiome")
library(data.table)

setwd("~/Desktop/SCCWRP")
#Read in algae data from SMC sites.
algaeDataSMCRaw <- read.table("AlgaeTax_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset columns of interest for the SMC sites.
algaeDataSMC <- algaeDataSMCRaw[,c(1,3,43,40,34)]
#Change the header name for station ID.
names(algaeDataSMC)[names(algaeDataSMC)=="Sample Station ID"]<-"SampleStationID"
#Calculate the relative abundance of algae data.
algaeDataSMC$RAbund <- with(algaeDataSMC,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
algaeDataSMC$OrganismType <- with(algaeDataSMC,"algae")

#Read in algae data from CEDEN sites.
algaeDataCEDENRaw <- read.table("AlgaeTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset columns of interest for the CEDEN sites.
algaeDataCEDEN <- algaeDataCEDENRaw[,c(6,11,36,32)]
#Fix the date format.
algaeDataCEDEN$SampleDate <- factor(gsub("-","/",algaeDataCEDEN$SampleDate))
#Remove rows with blank data.
algaeDataCEDEN <- na.omit(algaeDataCEDEN)
#Change names to uniforma schema.
names(algaeDataCEDEN)[names(algaeDataCEDEN)=="StationCode"]<-"SampleStationID"
names(algaeDataCEDEN)[names(algaeDataCEDEN)=="Species"]<-"FinalID"
#Determine the algal totals count column and make it a temporary dataframe.
tmp <- ddply(algaeDataCEDEN,"SampleStationID",numcolwise(sum))[c(1,2)]
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add algal totals count column to insect dataframe.
algaeDataCEDEN <- merge(algaeDataCEDEN,tmp,"SampleStationID")
#Calculate the relative abundance of algae data.
algaeDataCEDEN$RAbund <- with(algaeDataCEDEN,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
algaeDataCEDEN$OrganismType <- with(algaeDataCEDEN,"algae")

#The SWAMP data file is in a somewhat irregular format and this is accounted for
#when being read in.
algaeDataSWAMPRaw <- read.table("AlgaeTaxonomy_dnaSamples_SWAMP.csv", fill=TRUE,header=TRUE, sep=",",as.is=T,check.names=FALSE)
algaeDataSWAMP <- algaeDataSWAMPRaw[,c(6,8,97,90,73)]
algaeDataSWAMP <- na.omit(algaeDataSWAMP)
#Change names to uniforma schema.
names(algaeDataSWAMP)[names(algaeDataSWAMP)=="StationCode"]<-"SampleStationID"
#Calculate the relative abundance of algae data.
algaeDataSWAMP$RAbund <- with(algaeDataSWAMP,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
algaeDataSWAMP$OrganismType <- with(algaeDataSWAMP,"algae")

#Create merged algae data set.
algaeData <- do.call("rbind",list(algaeDataSMC,algaeDataSWAMP,algaeDataCEDEN))
algaeData <- na.omit(algaeData)
#Make the date format uniform
algaeData$SampleDate <- sub("-","/",algaeData$SampleDate)

#Read in insect data from SMC sites.
insectDataSMCRAW <- read.table("BugTax_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset columns of interest.
insectDataSMC <- insectDataSMCRAW[,c(1,3,9,6)]
#Change names to uniforma schema.
names(insectDataSMC)[names(insectDataSMC)=="Sample Station ID"]<-"SampleStationID"
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- ddply(insectDataSMC,"SampleStationID",numcolwise(sum))[c(1,2)]
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataSMC <- merge(insectDataSMC,tmp,"SampleStationID")
#Calculate the relative abundance of insect data.
insectDataSMC$RAbund <- with(insectDataSMC,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataSMC$OrganismType <- with(insectDataSMC,"invertebrate")

#Read in insect data from CEDEN sites.
insectDataCEDENRAW <- read.table("BugTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset columns of interest.
insectDataCEDEN <- insectDataCEDENRAW[,c(6,11,36,26)]
#Change names to uniforma schema.
names(insectDataCEDEN)[names(insectDataCEDEN)=="StationCode"]<-"SampleStationID"
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- ddply(insectDataCEDEN,"SampleStationID",numcolwise(sum))[c(1,2)]
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataCEDEN <- merge(insectDataCEDEN,tmp,"SampleStationID")
#Calculate the relative abundance of insect data.
insectDataCEDEN$RAbund <- with(insectDataCEDEN,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataCEDEN$OrganismType <- with(insectDataCEDEN,"invertebrate")

#Read in insect data from SWAMP sites.
insectDataSWAMPRAW <- read.table("BugTaxonomy_dnaSamples_SWAMP.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Subset columns of interest.
insectDataSWAMP <- insectDataSWAMPRAW[,c(1,8,97,90,73)]
#Change names to uniforma schema.
names(insectDataSWAMP)[names(insectDataSWAMP)=="Sample Station ID"]<-"SampleStationID"
#Calculate the relative abundance of insect data.
insectDataSWAMP$RAbund <- with(insectDataSWAMP,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataSWAMP$OrganismType <- with(insectDataSWAMP,"invertebrate")

#Create merged insect data set.
insectData <- do.call("rbind",list(insectDataSMC,insectDataSWAMP,insectDataCEDEN))
insectData <- na.omit(insectData)
#Make the date format uniform
insectData$SampleDate <- sub("-","/",insectData$SampleDate)

#Merge insect and algae data.
bioData <- do.call("rbind",list(insectData,algaeData))

#Read in chemical data for the test sites.
chemDataRAW <- read.table("Chem_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T)
#Subset columns.
chemData <- chemDataRAW[,-c(2,4:12,14,16,18,21:27)]
names(chemData)[names(chemData)=="Sample.Station.ID"]<-"SampleStationID"
chemData <- na.omit(chemData)

#Merge chemical and biological data.
biochemData <- join(bioData,chemData,by=c("SampleStationID","SampleDate"))

#Read in geospatial data.
GISDataRAW <- read.table("GIS_dnaSites.csv", header=TRUE, sep=",",as.is=T)
#Subset columns of interest.
GISData <- GISDataRAW[,-c(2:5,8:10,15)]
names(GISData)[names(GISData)=="Sample.Station.ID"]<-"SampleStationID"
names(GISData)[names(GISData)=="New_Lat"]<-"Latitude"
names(GISData)[names(GISData)=="New_Long"]<-"Longitude"
GISData <- na.omit(GISData)
#Merge geospatial data with biological observations.
GISBiochemData <- join(biochemData,GISData,by="SampleStationID")
GISBiochemData <- GISBiochemData[,-c(55:59,87:89,108)]
GISBiochemData <- na.omit(GISBiochemData)

#Calculate land usage index based on 1K, 5K, and catchment zone values.
GISBiochemData$LU_2011_1K <- with(GISBiochemData,Ag_2011_1K+CODE_21_2011_1K+URBAN_2011_1K)
