library("plyr")
library(dplyr)
library("ggplot2")
library(phyloseq)
library("ape")
library("vegan")
library("microbiome")
library(data.table)

setwd("~/Desktop/SCCWRP/DNA_Sites/SMC")
#Read in algae from SMC sites.
#Merge data into a single frame for biological observations.
algaeDataRAW <- read.table("AlgaeTax_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T)
#Subset columns of interest.
algaeData <- algaeDataRAW[,c(1,3,34,40,43)]
#Remove rows with missing data.
algaeData <- na.omit(algaeData)
#Calculate the relative abundance of algae data.
algaeData$AlgaeRAbund <- with(algaeData,BAResult/ActualOrganismCount)
colnames(algaeData) <- c("SampleStationID","SampleDate","AlgaeCountTotal","AlgaeID","AlgaeCount","AlgaeRAbund")
#Sort dataframe by Sample station ID and then by date.
algaeData <- algaeData[order(algaeData$SampleStationID),]
algaeData <- algaeData[order(algaeData$SampleDate),]
#Create unique ID combining the sample site ID with the sample date
algaeData$UniqueID <- with(algaeData,paste(SampleStationID,SampleDate))
#Initialize a diversity data frame in order to calculate and store alpha diversity
#measures for each unique ID in a data set.
algaeDiversity <- as.data.frame(matrix(ncol=2,nrow=length(unique(algaeData$UniqueID))))
colnames(algaeDiversity) <- c("UniqueID","ShannonIndex")
for(ID in unique(algaeData$UniqueID)){
  Shannon=diversity(algaeData[algaeData$UniqueID==ID,]$AlgaeRAbund)
  algaeDiversity[ID,] <- c(ID,Shannon)
}
algaeDiversity <- na.omit(algaeDiversity)
#Merge alpha diversity measures onto the algal data set.
algaeData <- join(algaeData,algaeDiversity,by="UniqueID")

#Read in insect data from SMC sites.
insectDataRAW <- read.table("BugTax_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T)
#Subset columns of interest.
insectData <- insectDataRAW[,c(1,3,6,9)]
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- ddply(insectDataRAW,"Sample.Station.ID",numcolwise(sum))[c(1,4)]
colnames(tmp) <- c("Sample.Station.ID","InsectCountTotal")
#Add insect totals count column to insect dataframe.
insectData <- merge(insectData,tmp,"Sample.Station.ID")
#Calculate the relative abundance of insect data.
insectData$InsectRAbund <- with(insectData,BAResult/InsectCountTotal)
colnames(insectData) <- c("SampleStationID","SampleDate","InsectID","InsectCount","InsectCountTotal","InsectRAbund")
insectData <- insectData[order(insectData$SampleStationID),]
insectData <- insectData[order(insectData$SampleDate),]

#Merge insect and algae data.
bioData <- merge(algaeData,insectData,by=c("SampleStationID","SampleDate"),all.x=TRUE, all.y=FALSE)
bioData <- na.omit(bioData)

#Read in chemical data for the test sites.
chemDataRAW <- read.table("Chem_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T)
#Subset columns.
chemData <- chemDataRAW[,-c(2,4:12,14,16,18,21:27)]
names(chemData)[names(chemData)=="Sample.Station.ID"]<-"SampleStationID"
chemData <- na.omit(chemData)

#Merge chemical and biological data.
biochemData <- join(bioData,chemData,by=c("SampleStationID","SampleDate"))

setwd("~/Desktop/SCCWRP/DNA_Sites/")
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
