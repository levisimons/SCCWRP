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
#Create unique ID combining the sample station ID and sampling date.
algaeDataSMC$UniqueID <- with(algaeDataSMC,paste(algaeDataSMC$SampleStationID,algaeDataSMC$SampleDate))

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
#Create unique ID combining the sample station ID and sampling date.
algaeDataCEDEN$UniqueID <- with(algaeDataCEDEN,paste(algaeDataCEDEN$SampleStationID,algaeDataCEDEN$SampleDate))

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
#Create unique ID combining the sample station ID and sampling date.
algaeDataSWAMP$UniqueID <- with(algaeDataSWAMP,paste(algaeDataSWAMP$SampleStationID,algaeDataSWAMP$SampleDate))

#Create merged algae data set.
algaeData <- do.call("rbind",list(algaeDataSMC,algaeDataSWAMP,algaeDataCEDEN))
algaeData <- na.omit(algaeData)
#Make the date format uniform
algaeData$SampleDate <- sub("-","/",algaeData$SampleDate)
algaeData$UniqueID <- sub("-","/",algaeData$UniqueID)

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
#Create unique ID combining the sample station ID and sampling date.
insectDataSMC$UniqueID <- with(insectDataSMC,paste(insectDataSMC$SampleStationID,insectDataSMC$SampleDate))

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
#Create unique ID combining the sample station ID and sampling date.
insectDataCEDEN$UniqueID <- with(insectDataCEDEN,paste(insectDataCEDEN$SampleStationID,insectDataCEDEN$SampleDate))

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
#Create unique ID combining the sample station ID and sampling date.
insectDataSWAMP$UniqueID <- with(insectDataSWAMP,paste(insectDataSWAMP$SampleStationID,insectDataSWAMP$SampleDate))

#Create merged insect data set.
insectData <- do.call("rbind",list(insectDataSMC,insectDataSWAMP,insectDataCEDEN))
insectData <- na.omit(insectData)
#Make the date format uniform
insectData$SampleDate <- sub("-","/",insectData$SampleDate)
insectData$UniqueID <- sub("-","/",insectData$UniqueID)

#Merge insect and algae data.
bioData <- do.call("rbind",list(insectData,algaeData))

#Read in chemical data for the test sites.
chemDataRAW <- read.table("Chem_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset columns.
chemData <- chemDataRAW[,-c(2,4:12,14,16,18,21:27)]
names(chemData)[names(chemData)=="Sample Station ID"]<-"SampleStationID"
#Create unique ID combining the sample station ID and sampling date.
chemData$UniqueID <- with(chemData,paste(chemData$SampleStationID,chemData$SampleDate))
chemData <- na.omit(chemData)

#Merge chemical and biological data.
biochemData <- join(bioData,chemData,by="UniqueID")

#Read in geospatial data.
GISDataRAW <- read.table("GIS_dnaSites.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset columns of interest.
GISData <- GISDataRAW[,-c(2:5,8:10,15)]
names(GISData)[names(GISData)=="Sample Station ID"]<-"SampleStationID"
names(GISData)[names(GISData)=="New_Lat"]<-"Latitude"
names(GISData)[names(GISData)=="New_Long"]<-"Longitude"
GISData <- na.omit(GISData)
#Merge geospatial data with biological observations.
GISBiochemData <- join(biochemData,GISData,by="SampleStationID")
#Create unique ID combining the sample station ID and sampling date.
GISBiochemData$UniqueID <- with(GISBiochemData,paste(GISBiochemData$SampleStationID,GISBiochemData$SampleDate))
GISBiochemData <- na.omit(GISBiochemData)
GISBiochemData <- GISBiochemData[,-c(9:10,22:30,56:59,87:98,108)]

#Calculate land usage index based on 1K, 5K, and catchment zone values.
GISBiochemData$LU_2011_1K <- with(GISBiochemData,Ag_2011_1K+CODE_21_2011_1K+URBAN_2011_1K)
GISBiochemData$LU_2011_5K <- with(GISBiochemData,Ag_2011_5K+CODE_21_2011_5K+URBAN_2011_5K)
GISBiochemData$LU_2011_WS <- with(GISBiochemData,Ag_2011_WS+CODE_21_2011_WS+URBAN_2011_WS)

#Subset site data based on land usage index within 1K catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISBiochemDataLD1K <- GISBiochemData[which(GISBiochemData$LU_2011_1K < 5),]
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISBiochemDataMD1K <- GISBiochemData[which(GISBiochemData$LU_2011_1K < 15 & GISBiochemData$LU_2011_1K >= 5),]
#HD = low disturbance.  Land usage index is greater than 15%.
GISBiochemDataHD1K <- GISBiochemData[which(GISBiochemData$LU_2011_1K >= 15),]

#Subset site data based on land usage index within 5K catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISBiochemDataLD5K <- GISBiochemData[which(GISBiochemData$LU_2011_5K < 5),]
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISBiochemDataMD5K <- GISBiochemData[which(GISBiochemData$LU_2011_5K < 15 & GISBiochemData$LU_2011_5K >= 5),]
#HD = low disturbance.  Land usage index is greater than 15%.
GISBiochemDataHD5K <- GISBiochemData[which(GISBiochemData$LU_2011_5K >= 15),]

#Subset site data based on land usage index within the full water drainage catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISBiochemDataLDWS <- GISBiochemData[which(GISBiochemData$LU_2011_WS < 5),]
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISBiochemDataMDWS <- GISBiochemData[which(GISBiochemData$LU_2011_WS < 15 & GISBiochemData$LU_2011_WS >= 5),]
#HD = low disturbance.  Land usage index is greater than 15%.
GISBiochemDataHDWS <- GISBiochemData[which(GISBiochemData$LU_2011_WS >= 15),]

#Initialize a data frame where the rows are all of the unique taxa for a given
#subset of the data.
numrow = length(unique(GISBiochemDataHD1K$FinalID))
numcol = length(unique(GISBiochemDataHD1K$UniqueID))
HD1K <- as.data.frame(unique(GISBiochemDataHD1K$FinalID))
names(HD1K)[names(HD1K)=="unique(GISBiochemDataHD1K$FinalID)"]<-"FinalID"
HD1K[order(HD1K$FinalID),]

#Add the relative taxa abundances by column to a new dataframe.
#The rows are the unique taxa in a given subset of data.
i=1
for(ID in unique(GISBiochemDataHD1K$UniqueID)){
  print(ID)
  i=i+1
  print(i)
  tmp <- GISBiochemDataHD1K[which(GISBiochemDataHD1K$UniqueID==ID),][,c(4,6)]
  tmp[order(tmp$FinalID),]
  tmp <- tmp[!duplicated(tmp),]
  print(length(tmp$FinalID))
  names(tmp)[names(tmp)=="RAbund"]<-paste("RAbund",ID)
  HD1K <- join(HD1K,tmp,by="FinalID")
}
#Output dataframe for use in eLSA.
#Note that the the data needs to have at least two location replicates per time point
#and that the number of replicates per time point needs to be uniform.
#This may involve subsampling data depending on the variation in the number of replicates per time point.
names(HD1K)[names(HD1K)=="FinalID"]<-"#FinalID"
write.table(HD1K,"HD1K.tsv",quote=FALSE,sep="\t",row.names = FALSE)

#Read in eLSA output.
#Compute network statistics of the likeliest association networks between taxa.
library(igraph)
library(network)
networkdata <- read.table("HD1KNetwork.txt",header=TRUE, sep="\t",as.is=T)
names(networkdata)[names(networkdata)=="PCC"]<-"weight"
networkgraph=graph.data.frame(networkdata,directed=FALSE)
plot(networkgraph,layout=layout.circle(networkgraph),edge.width=E(networkgraph)$weight*10,edge.color=ifelse(E(networkgraph)$weight > 0, "blue","red"))
# Calculate the average network path length
mean_distance(networkgraph)
# Calculate the clustering coefficient
transitivity(networkgraph)
# Generate adjacency matrix of relative taxa abundance correlations
adj= as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
# Get the number of unique network edges
0.5*network.edgecount(adj)
# Get the network density.
network.density(adj)

library(vcd)
library(MASS)
# Get degree distribution of network.
DDN <- degree(networkgraph)
# Fit a poisson distribution to the link distribution of the network
poissonFit <- fitdistr(DDN,"Poisson")
# Get the value of lambda for the Poisson distribution
coef(poissonFit)
# Get the log-likelihood for this fit
logLik(poissonFit)
