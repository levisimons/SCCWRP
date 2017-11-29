library("plyr")
library(dplyr)
library("ggplot2")
library(phyloseq)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)

setwd("~/Desktop/SCCWRP")
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
#Create unique ID combining the sample station ID and sampling date.
algaeDataSMC$UniqueID <- with(algaeDataSMC,paste(algaeDataSMC$SampleStationID,"SMC",algaeDataSMC$SampleDate))
#Find sampling year.
algaeDataSMC$Year <- with(algaeDataSMC,lapply(strsplit(algaeDataSMC$SampleDate,split="/"),"[",3))

#Read in algae data from CEDEN sites.
algaeDataCEDENRaw <- read.table("AlgaeTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
algaeDataCEDEN <- filter(algaeDataCEDENRaw, CollectionReplicate==1)
#Subset columns of interest for the CEDEN sites.
algaeDataCEDEN <- algaeDataCEDENRaw[,c(6,11,36,32)]
#Fix the date format.
algaeDataCEDEN$SampleDate <- factor(gsub("-","/",algaeDataCEDEN$SampleDate))
#Change names to uniforma schema.
names(algaeDataCEDEN)[names(algaeDataCEDEN)=="StationCode"]<-"SampleStationID"
names(algaeDataCEDEN)[names(algaeDataCEDEN)=="Species"]<-"FinalID"
#Determine the algal totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,algaeDataCEDEN))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add algal totals count column to algae dataframe.
algaeDataCEDEN <- merge(algaeDataCEDEN,tmp,"SampleStationID")
#Calculate the relative abundance of algae data.
algaeDataCEDEN$Measurement <- with(algaeDataCEDEN,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
algaeDataCEDEN$MeasurementType <- with(algaeDataCEDEN,"Algal relative abundance")
#Create unique ID combining the sample station ID and sampling date.
algaeDataCEDEN$UniqueID <- with(algaeDataCEDEN,paste(algaeDataCEDEN$SampleStationID,"CEDEN",algaeDataCEDEN$SampleDate))
#Find sampling year.
algaeDataCEDEN$Year <- with(algaeDataCEDEN,lapply(strsplit(as.character(algaeDataCEDEN$SampleDate),split="/"),"[",1))

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
#Create unique ID combining the sample station ID and sampling date.
algaeDataSWAMP$UniqueID <- with(algaeDataSWAMP,paste(algaeDataSWAMP$SampleStationID,"SWAMP",algaeDataSWAMP$SampleDate))
#Find sampling year.
algaeDataSWAMP$Year <- with(algaeDataSWAMP,lapply(strsplit(algaeDataSWAMP$SampleDate,split="/"),"[",3))

#Create merged algae data set.
algaeData <- do.call("rbind",list(algaeDataSMC,algaeDataSWAMP,algaeDataCEDEN))
algaeData <- na.omit(algaeData)
#Make the date format uniform
algaeData$SampleDate <- sub("-","/",algaeData$SampleDate)
algaeData$UniqueID <- sub("-","/",algaeData$UniqueID)

#Read in insect data from SMC sites.
insectDataSMCRAW <- read.table("BugTax_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
insectDataSMC <- filter(insectDataSMCRAW, FieldReplicate==1)
#Subset columns of interest.
insectDataSMC <- insectDataSMCRAW[,c(1,3,9,6)]
#Change names to uniforma schema.
names(insectDataSMC)[names(insectDataSMC)=="Sample Station ID"]<-"SampleStationID"
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,insectDataSMC))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataSMC <- merge(insectDataSMC,tmp,"SampleStationID")
#Calculate the relative abundance of insect data.
insectDataSMC$Measurement <- with(insectDataSMC,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataSMC$MeasurementType <- with(insectDataSMC,"Invertebrate relative abundances")
#Create unique ID combining the sample station ID and sampling date.
insectDataSMC$UniqueID <- with(insectDataSMC,paste(insectDataSMC$SampleStationID,"SMC",insectDataSMC$SampleDate))
#Find sampling year.
insectDataSMC$Year <- with(insectDataSMC,lapply(strsplit(insectDataSMC$SampleDate,split="/"),"[",3))

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
#Create unique ID combining the sample station ID and sampling date.
insectDataCEDEN$UniqueID <- with(insectDataCEDEN,paste(insectDataCEDEN$SampleStationID,"CEDEN",insectDataCEDEN$SampleDate))
#Find sampling year.
insectDataCEDEN$Year <- with(insectDataCEDEN,lapply(strsplit(as.character(insectDataCEDEN$SampleDate),split="-"),"[",1))

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
#Create unique ID combining the sample station ID and sampling date.
insectDataSWAMP$UniqueID <- with(insectDataSWAMP,paste(insectDataSWAMP$SampleStationID,"SWAMP",insectDataSWAMP$SampleDate))
#Find sampling year.
insectDataSWAMP$Year <- with(insectDataSWAMP,lapply(strsplit(as.character(insectDataSWAMP$SampleDate),split="/"),"[",3))

#Create merged insect data set.
insectData <- do.call("rbind",list(insectDataSMC,insectDataSWAMP,insectDataCEDEN))
insectData <- na.omit(insectData)
#Make the date format uniform
insectData$SampleDate <- sub("-","/",insectData$SampleDate)
insectData$UniqueID <- sub("-","/",insectData$UniqueID)

#Merge insect and algae data.
bioData <- do.call("rbind",list(insectData,algaeData))
bioData <- bioData[,-c(3,5)]

#Read in chemical data for the SMC test sites.
chemDataSMCRAW <- read.table("Chem_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
chemDataSMC <- filter(chemDataSMCRAW, FieldReplicate==1)
#Subset columns.
chemDataSMC <- chemDataSMCRAW[,c(1,3,13,17,15)]
#Introduce common naming schema.
names(chemDataSMC)[names(chemDataSMC)=="Sample Station ID"]<-"SampleStationID"
names(chemDataSMC)[names(chemDataSMC)=="AnalyteName"]<-"FinalID"
names(chemDataSMC)[names(chemDataSMC)=="Result"]<-"Measurement"
names(chemDataSMC)[names(chemDataSMC)=="Unit"]<-"MeasurementType"
#Create unique ID combining the sample station ID and sampling date.
chemDataSMC$UniqueID <- with(chemDataSMC,paste(chemDataSMC$SampleStationID,"SMC",chemDataSMC$SampleDate))
#Find sampling year.
chemDataSMC$Year <- with(chemDataSMC,lapply(strsplit(chemDataSMC$SampleDate,split="/"),"[",3))

#Read in chemical data for the CEDEN test sites.
chemDataCEDENRAW <- read.table("Chem_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
chemDataCEDEN <- filter(chemDataCEDENRAW, CollectionReplicate==1)
#Subset the data.
chemDataCEDEN <- chemDataCEDENRAW[,c(5,6,18,20,19)]
#Introduce common naming schema.
names(chemDataCEDEN)[names(chemDataCEDEN)=="StationCode"]<-"SampleStationID"
names(chemDataCEDEN)[names(chemDataCEDEN)=="Analyte"]<-"FinalID"
names(chemDataCEDEN)[names(chemDataCEDEN)=="Result"]<-"Measurement"
names(chemDataCEDEN)[names(chemDataCEDEN)=="Unit"]<-"MeasurementType"
#Create unique ID combining the sample station ID and sampling date.
chemDataCEDEN$UniqueID <- with(chemDataCEDEN,paste(chemDataCEDEN$SampleStationID,"CEDEN",chemDataCEDEN$SampleDate))
#Find sampling year.
chemDataCEDEN$Year <- with(chemDataCEDEN,lapply(strsplit(chemDataCEDEN$SampleDate,split="-"),"[",1))

#Read in chemical data for the SWAMP test sites.
chemDataSWAMPRAW <- read.table("Chem_dnaSamples_SWAMP.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Subset only replicate 1
chemDataSWAMP <- filter(chemDataSWAMPRAW, Replicate==1)
#Subset the data.
chemDataSWAMP <- chemDataSWAMPRAW[,c(1,8,70,89,74)]
#Introduce common naming schema.
names(chemDataSWAMP)[names(chemDataSWAMP)=="Sample Station ID"]<-"SampleStationID"
names(chemDataSWAMP)[names(chemDataSWAMP)=="AnalyteName"]<-"FinalID"
names(chemDataSWAMP)[names(chemDataSWAMP)=="Result"]<-"Measurement"
names(chemDataSWAMP)[names(chemDataSWAMP)=="UnitName"]<-"MeasurementType"
#Create unique ID combining the sample station ID and sampling date.
chemDataSWAMP$UniqueID <- with(chemDataSWAMP,paste(chemDataSWAMP$SampleStationID,"SWAMP",chemDataSWAMP$SampleDate))
#Find sampling year.
chemDataSWAMP$Year <- with(chemDataSWAMP,lapply(strsplit(chemDataSWAMP$SampleDate,split="/"),"[",3))
#Remove last erroneous row from data frame.
chemDataSWAMP <- head(chemDataSWAMP,-1)
#Create merged chemical data frame.
chemData <- do.call("rbind",list(chemDataSMC,chemDataSWAMP,chemDataCEDEN))
chemData <- na.omit(chemData)
#Make the date format uniform
chemData$SampleDate <- sub("-","/",chemData$SampleDate)
chemData$UniqueID <- sub("-","/",chemData$UniqueID)

#Merge chemical and biological data.
biochemData <- do.call("rbind",list(bioData,chemData))
#Sort bio-chem data by year.
biochemData <- biochemData[order(as.numeric(biochemData$Year)),]

#Read in geospatial data.
GISDataRAW <- read.table("GIS_dnaSites.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset columns of interest.
GISData <- GISDataRAW[,-c(2:5,8:10,15)]
names(GISData)[names(GISData)=="Sample Station ID"]<-"SampleStationID"
names(GISData)[names(GISData)=="New_Lat"]<-"Latitude"
names(GISData)[names(GISData)=="New_Long"]<-"Longitude"

#Merge geospatial data with biological observations.
GISBiochemData <- join(biochemData,GISData,by="SampleStationID")
GISBiochemData <- GISBiochemData[,-c(10:11,14:22,47:51,82:90,100)]

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

#Initialize a data frame where the rows are all of the unique measurements for a given
#subset of the data.
#Order the data frame by measurement name.
eLSAInput <- as.data.frame(unique(GISBiochemDataHD1K$FinalID))
colnames(eLSAInput)<-c("FinalID")
eLSAInput <- as.data.frame(eLSAInput[order(eLSAInput$FinalID),])
colnames(eLSAInput)<-c("FinalID")
#eLSAInput <- eLSAInput %>% group_by(FinalID,UniqueID) %>% summarise_all((funs(mean)))

#Add the relative taxa abundances by column to a new dataframe.
#The rows are the unique taxa in a given subset of data.
for(ID in unique(GISBiochemDataHD1K$UniqueID)){
  tmp <- filter(GISBiochemDataHD1K, UniqueID == ID)[,c(3,4,6)]
  tmp <- as.data.frame(tmp[order(tmp$FinalID),])
  tmp <- tmp[-c(3)]
  colnames(tmp)<-c("FinalID",paste("Measurement",ID,sep=" "))
  eLSAInput <- join(eLSAInput,tmp,by="FinalID")
  #eLSAInput <- join(eLSAInput,tmp,by="FinalID")
  #for(ID in unique(tmp$UniqueID)){
    #label=paste("Measurement",year,"Replicate",i,sep=" ")
    #print(label)
    #tmp = filter(tmp, UniqueID==ID)
    #colnames(tmp)<-c("FinalID",label,"UniqueID")
    #eLSAInput <- join(eLSAInput,tmp,by="FinalID")
    #eLSAInput$UniqueID <- NULL
  #}
}
eLSAInput$FinalID=as.character(eLSAInput$FinalID)
eLSAInput <- eLSAInput %>% group_by(FinalID) %>% summarise_if(is.numeric,mean,na.rm=TRUE)

eLSAInput[is.na(eLSAInput)] <- "NA"

#Output dataframe for use in eLSA.
#Note that the the data needs to have at least two location replicates per time point
#and that the number of replicates per time point needs to be uniform.
#This may involve subsampling data depending on the variation in the number of replicates per time point.
names(eLSAInput)[names(eLSAInput)=="FinalID"]<-"#FinalID"
write.table(eLSAInput,"eLSAInput.txt",quote=FALSE,sep="\t",row.names = FALSE)

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
