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
#Force a uniform date format
chemDataSMC$SampleDate <- mdy(chemDataSMC$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
chemDataSMC$UniqueID <- with(chemDataSMC,paste(chemDataSMC$SampleStationID,"SMC",chemDataSMC$SampleDate))
#Find sampling year.
chemDataSMC$Year <- year(chemDataSMC$SampleDate)

#Read in chemical data for the CEDEN test sites.
chemDataCEDENRAW <- read.csv("Chem_dnaSites_CEDEN.csv")
#Subset only replicate 1
chemDataCEDEN <- filter(chemDataCEDENRAW, CollectionReplicate==1)
#Subset the data.
chemDataCEDEN <- chemDataCEDEN[,c(5,6,18,20,19)]
#Introduce common naming schema.
names(chemDataCEDEN)[names(chemDataCEDEN)=="StationCode"]<-"SampleStationID"
names(chemDataCEDEN)[names(chemDataCEDEN)=="Analyte"]<-"FinalID"
names(chemDataCEDEN)[names(chemDataCEDEN)=="Result"]<-"Measurement"
names(chemDataCEDEN)[names(chemDataCEDEN)=="Unit"]<-"MeasurementType"
#Force a uniform date format
chemDataCEDEN$SampleDate <- as.Date(chemDataCEDEN$SampleDate,format="%m/%d/%y")
#Create unique ID combining the sample station ID and sampling date.
chemDataCEDEN$UniqueID <- with(chemDataCEDEN,paste(chemDataCEDEN$SampleStationID,"CEDEN",chemDataCEDEN$SampleDate))
#Find sampling year.
chemDataCEDEN$Year <- year(chemDataCEDEN$SampleDate)

#Read in chemical data for the SWAMP test sites.
chemDataSWAMPRAW <- read.csv("Chem_dnaSamples_SWAMP.csv")
#Subset only replicate 1
chemDataSWAMP <- filter(chemDataSWAMPRAW,Replicate==1)
#Subset the data.
chemDataSWAMP <- chemDataSWAMP[,c(1,8,70,89,74)]
#Introduce common naming schema.
names(chemDataSWAMP)[names(chemDataSWAMP)=="Sample.Station.ID"]<-"SampleStationID"
names(chemDataSWAMP)[names(chemDataSWAMP)=="AnalyteName"]<-"FinalID"
names(chemDataSWAMP)[names(chemDataSWAMP)=="Result"]<-"Measurement"
names(chemDataSWAMP)[names(chemDataSWAMP)=="UnitName"]<-"MeasurementType"
#Force a uniform date format
chemDataSWAMP$SampleDate <- mdy(chemDataSWAMP$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
chemDataSWAMP$UniqueID <- with(chemDataSWAMP,paste(chemDataSWAMP$SampleStationID,"SWAMP",chemDataSWAMP$SampleDate))
#Find sampling year.
chemDataSWAMP$Year <- year(chemDataSWAMP$SampleDate)

#Create merged chemical data frame.
chemData <- do.call("rbind",list(chemDataSMC,chemDataSWAMP,chemDataCEDEN))
chemData <- na.omit(chemData)

#Merge chemical and biological data.
bioChemData <- do.call("rbind",list(bioData,chemData))
#Sort bio-chem data by year.
bioChemData <- bioChemData[order(as.numeric(bioChemData$Year)),]

#Read in geospatial data.
GISDataRAW <- read.table("GIS_dnaSites.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset columns of interest.
GISData <- GISDataRAW[,-c(1,3:5,8:10,15)]
names(GISData)[names(GISData)=="StationCode"]<-"SampleStationID"
names(GISData)[names(GISData)=="New_Lat"]<-"Latitude"
names(GISData)[names(GISData)=="New_Long"]<-"Longitude"

#Merge geospatial data with biological observations.
GISBiochemData <- join(bioChemData,GISData,by="SampleStationID")
#GISBiochemData <- GISBiochemData[,-c(10:11,14:22,47:51,82:90,100)]
#Sort merged data set by year then measurement name.
GISBiochemData <- as.data.frame(GISBiochemData[order(as.numeric(GISBiochemData$Year),as.character(GISBiochemData$FinalID)),])

#Filter out erroneous negative values for physical parameter data.
chemID <- unique(chemData$FinalID)
GISBiochemData <- filter(GISBiochemData,Measurement>=0)

#Filter out duplicate rows.
GISBiochemData <- GISBiochemData[!duplicated(GISBiochemData),]

#Read in the covariant component of the networks generated via land usage subsetting
covariantNetwork <- read.csv("SCCWRPNetworkAnalysisLandP01Covariant.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
covariantNetwork[covariantNetwork=="#DIV/0!" | covariantNetwork=="#NUM!"] <- "NA"
covariantNetwork <- cbind(covariantNetwork[,1:3],data.frame(sapply(covariantNetwork[,4:ncol(covariantNetwork)],function(x) as.numeric(as.character(x)))))

#Read in the covariant component of the networks generated via land usage subsetting
contravariantNetwork <- read.csv("SCCWRPNetworkAnalysisLandP01Contravariant.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
contravariantNetwork[contravariantNetwork=="#DIV/0!" | contravariantNetwork=="#NUM!"] <- "NA"
contravariantNetwork <- cbind(contravariantNetwork[,1:3],data.frame(sapply(contravariantNetwork[,4:ncol(contravariantNetwork)],function(x) as.numeric(as.character(x)))))

#Read in California Stream Condition Index data
csciData <- read.csv("csci_scored_sites_tbl.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
csciData <- filter(csciData,CSCI!="NA")
names(csciData)[names(csciData)=="StationCode"]<-"SampleStationID"
#Add in qualifier columns based on California Streams Condition Index.  The cutoff is 0.79.
csciData$CSCIQualifier <- ifelse(csciData$CSCI >= 0.79, "Healthy","Disturbed")
csciData$CSCIQualNum <- ifelse(csciData$CSCI >= 0.79, 1,0)

#Calculate land usage index based on 1K, 5K, and catchment zone values.
#Use land usage data from 2011.
GISBiochemData$LU_2011_1K <- with(GISBiochemData,Ag_2011_1K+CODE_21_2011_1K+URBAN_2011_1K)
GISBiochemData$LU_2011_5K <- with(GISBiochemData,Ag_2011_5K+CODE_21_2011_5K+URBAN_2011_5K)
GISBiochemData$LU_2011_WS <- with(GISBiochemData,Ag_2011_WS+CODE_21_2011_WS+URBAN_2011_WS)

#Filter out data without a land usage index.
GISBiochemData <- filter(GISBiochemData,LU_2011_1K!="NA")
GISBiochemData <- filter(GISBiochemData,LU_2011_5K!="NA")
GISBiochemData <- filter(GISBiochemData,LU_2011_WS!="NA")

#Subset site data based on land usage index within 1K catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISBiochemDataLD1K <- GISBiochemData[which(GISBiochemData$LU_2011_1K < 5),]
GISBiochemDataLD1K$LUCategory <- "LD1K"
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISBiochemDataMD1K <- GISBiochemData[which(GISBiochemData$LU_2011_1K < 15 & GISBiochemData$LU_2011_1K >= 5),]
GISBiochemDataMD1K$LUCategory <- "MD1K"
#HD = low disturbance.  Land usage index is greater than 15%.
GISBiochemDataHD1K <- GISBiochemData[which(GISBiochemData$LU_2011_1K >= 15),]
GISBiochemDataHD1K$LUCategory <- "HD1K"

#Subset site data based on land usage index within 5K catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISBiochemDataLD5K <- GISBiochemData[which(GISBiochemData$LU_2011_5K < 5),]
GISBiochemDataLD5K$LUCategory <- "LD5K"
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISBiochemDataMD5K <- GISBiochemData[which(GISBiochemData$LU_2011_5K < 15 & GISBiochemData$LU_2011_5K >= 5),]
GISBiochemDataMD5K$LUCategory <- "MD5K"
#HD = low disturbance.  Land usage index is greater than 15%.
GISBiochemDataHD5K <- GISBiochemData[which(GISBiochemData$LU_2011_5K >= 15),]
GISBiochemDataHD5K$LUCategory <- "HD5K"

#Subset site data based on land usage index within the full water drainage catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISBiochemDataLDWS <- GISBiochemData[which(GISBiochemData$LU_2011_WS < 5),]
GISBiochemDataLDWS$LUCategory <- "LDWS"
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISBiochemDataMDWS <- GISBiochemData[which(GISBiochemData$LU_2011_WS < 15 & GISBiochemData$LU_2011_WS >= 5),]
GISBiochemDataMDWS$LUCategory <- "MDWS"
#HD = low disturbance.  Land usage index is greater than 15%.
GISBiochemDataHDWS <- GISBiochemData[which(GISBiochemData$LU_2011_WS >= 15),]
GISBiochemDataHDWS$LUCategory <- "HDWS"

#Merge land usage subsets back together for later analytical tools.
GISBiochemData <- do.call("rbind",list(GISBiochemDataLD1K,GISBiochemDataMD1K,GISBiochemDataHD1K,GISBiochemDataLD5K,GISBiochemDataMD5K,GISBiochemDataHD5K,GISBiochemDataLDWS,GISBiochemDataMDWS,GISBiochemDataHDWS))

#Merge in CSCI data and network statistics data
GISBiochemcsciData <- join(GISBiochemData,csciData,by="SampleStationID")
GISBiochemcsciDataCovariant <- join(GISBiochemcsciData,covariantNetwork,by="LUCategory")
GISBiochemcsciDataCovariant <- filter(GISBiochemcsciDataCovariant,GISBiochemcsciDataCovariant$CSCI!="NA" & GISBiochemcsciDataCovariant$l_rL!="NA")
GISBiochemcsciDataContravariant <- join(GISBiochemcsciData,contravariantNetwork,by="LUCategory")
GISBiochemcsciDataContravariant <- filter(GISBiochemcsciDataContravariant,GISBiochemcsciDataContravariant$CSCI!="NA" & GISBiochemcsciDataContravariant$l_rL!="NA")

#If you want to subselect data by geographic selection scales (1km, 5km, or watershed)
GISBiochemcsciDataContravariant <- filter(GISBiochemcsciDataContravariant,SelectionZone=="Watershed")
GISBiochemcsciDataCovariant <- filter(GISBiochemcsciDataCovariant,SelectionZone=="Watershed")

#Add in a scaled ratio in order to help plot the log ratio of network path lengths.
# -30*l_rL = rL for mapping.
GISBiochemcsciDataCovariant$rL <- -30*(GISBiochemcsciDataCovariant$l_rL)
GISBiochemcsciDataContravariant$rL <- -30*(GISBiochemcsciDataContravariant$l_rL)

#Generate map of data for a given chemical parameter in California.
library(ggmap)
library(maps)
library(mapdata)
dev.off()
MapCoordinates <- data.frame(GISBiochemcsciDataContravariant$rL,GISBiochemcsciDataContravariant$CSCIQualifier,GISBiochemcsciDataContravariant$Longitude,GISBiochemcsciDataContravariant$Latitude)
colnames(MapCoordinates) = c('SECND','CSCIQualifier','lon','lat')
MapCoordinates <- na.omit(MapCoordinates)
mapBoundaries <- make_bbox(lon=MapCoordinates$lon,lat=MapCoordinates$lat,f=0.1)
CalMap <- get_map(location=mapBoundaries,maptype="satellite",source="google")
CalMap <- ggmap(CalMap)+geom_point(data = MapCoordinates, mapping = aes(x = lon, y = lat, color = CSCIQualifier, size=SECND))+ggtitle("Stream disturbance and the Scaled Ecological Contravariant Network Distance (SECND)\n1 km selection scales",subtitle="SECND = -30*log(L_contravariant/L_random)")
CalMap

#Logistic regression between network parameters and CSCI
library(aod)
library(rcompanion)
logReg <- glm(formula = CSCIQualNum ~ rL, data = GISBiochemcsciDataContravariant, family=binomial(link="logit"))
#Determine pseudo-r^2 and p for logistic model.
nagelkerke(logReg)
#Summary statistics for logistic model.
#log(p/(1-p)) = k*x + b.
#The probability of the binomial model being in a state of 1, versus 0, is p.
#First decide a value of p, such as 0.95 for 95% probability.
#The estimate of variable, in this case l_rL, is k.
#The estimate of the intercept is b.
#Solve for the value of the target variable which gives a probability of p that the binomial model is in a state of 1 versus 0.
#x = -1*([log(p/(1-p))-b]/k) for a given value of p, b, and k.
summary(logReg)

#Check for correlations between network parameters and the CSCI
v <- cbind(as.numeric(GISBiochemcsciDataCovariant$CSCI),as.numeric(GISBiochemcsciDataCovariant$l_rL))
v <- na.omit(v)
landCor <- round(cor(v[,1],v[,2],method="pearson"),4)
landP <- round(cor.test(v[,1],v[,2],method="pearson")$p.value,4)
dev.off()
plot(v[,1],v[,2],main=paste("Covariant l_rL vs. CSCI, All selection scales","\n","Pearson r = ",landCor,"p = ",landP),ylab="Covariant l_rL",xlab="CSCI")
abline(lm(v[,2]~v[,1]),col="red")
v <- cbind(as.numeric(GISBiochemcsciDataContravariant$CSCI),as.numeric(GISBiochemcsciDataContravariant$rL))
v <- na.omit(v)
cor(v[,1],v[,2],method="pearson")
landCor <- round(cor(v[,1],v[,2],method="pearson"),4)
landP <- round(cor.test(v[,1],v[,2],method="pearson")$p.value,4)
dev.off()
plot(v[,1],v[,2],main=paste("SECND vs. CSCI, 1km selection scale","\n","Pearson r = ",landCor,"p = ",landP),ylab="SECND",xlab="CSCI")
abline(lm(v[,2]~v[,1]),col="red")

#Find the average value of the parameters most strongly correlated to the top
#two principal components which describe variations in chemical parameter space
#given changes in land usage intensity.  These averages will be used to split
#the site data according to sites above or below the average parameter values.
#parameter1 = strongly correlated parameter to axis composed of
#principal components 1 and 2.
i=0
parameter1 <- data.frame()
parameter1Name <- "Nitrate + Nitrite as N"
for(site in unique(GISBiochemData$UniqueID)){
  GISBiochemDataSite <- GISBiochemData[GISBiochemData$UniqueID == site,]
  if(parameter1Name %in% GISBiochemDataSite$FinalID){
    i=i+1
    tmp1 <- GISBiochemDataSite[which(GISBiochemDataSite$FinalID==parameter1Name),]
    if(tmp1$Measurement[1]>0){
      parameter1 <- rbind(parameter1,tmp1$Measurement[1])
    }
    print(paste(tmp1))
  }
}

parameter1Ave <- colMeans(as.data.frame(rowMeans(parameter1,na.rm=TRUE)),na.rm=TRUE)
print(paste("There are",i,"rows where this parameter is measured.",sep=" "))
print(paste("Average value",parameter1Name,":",parameter1Ave,sep=" "))

#Split site data based on the average value of the two parameters.
#Given the most significant chemical factors related to changes in land usage
#subset site data based on those values.
#H1 = high parameter 1 concentration.  L1 = low parameter 1 concentration.
GISBiochemDataH1 <- data.frame()
GISBiochemDataL1 <- data.frame()
for(site in unique(GISBiochemData$UniqueID)){
  GISBiochemDataSite <- GISBiochemData[GISBiochemData$UniqueID == site,]
  if(parameter1Name %in% GISBiochemDataSite$FinalID){
    tmp1 <- GISBiochemDataSite[which(GISBiochemDataSite$FinalID==parameter1Name),]
    if(tmp1$Measurement[1]>parameter1Ave){
      print(paste("H1: ",site,parameter1Name,tmp1$Measurement[1]))
      GISBiochemDataH1 <- rbind(GISBiochemDataH1,GISBiochemDataSite)
      }
    if(tmp1$Measurement[1]<parameter1Ave){
      print(paste("L1: ",site,parameter1Name,tmp1$Measurement[1]))
      GISBiochemDataL1 <- rbind(GISBiochemDataL1,GISBiochemDataSite)
      }
    }
}

#Select a geographic subset data frame from the total merged data set.
selected <- GISBiochemDataH1
suffix <- "HighNitrateNitrite"

#Initialize a data frame where the rows are all of the unique measurements for a given
#subset of the data.
#Order the data frame by measurement name.
eLSAInput <- as.data.frame(unique(selected$FinalID))
colnames(eLSAInput)<-c("FinalID")
eLSAInput <- as.data.frame(eLSAInput[order(as.character(eLSAInput$FinalID)),])
colnames(eLSAInput)<-c("FinalID")

#Add the relative taxa abundances by column to a new dataframe.
#The rows are the unique taxa in a given subset of data.
for(ID in unique(selected$UniqueID)){
  tmp <- filter(selected, UniqueID == ID)[,c(3,4,6)]
  tmp <- as.data.frame(tmp[order(tmp$FinalID),])
  tmp <- tmp[-c(3)]
  colnames(tmp)<-c("FinalID",paste("Measurement",ID,sep=" "))
  eLSAInput <- join(eLSAInput,tmp,by="FinalID")
  eLSAInput$FinalID=as.character(eLSAInput$FinalID)
  eLSAInput <- eLSAInput %>% group_by(FinalID) %>% summarise_if(is.numeric,mean,na.rm=TRUE)
  print(ID)
}

eLSAInput[is.na(eLSAInput)] <- "NA"

#Determine the number of time points in the eLSA input file.
spotNum = length(unique(selected$Year))
#Determine the number of replicates per time point in the eLSA input file.
#In order to ensure a uniform number of replicates per year this needs to
#be the maximum number of replicates for all of the years available.
repMax = 0
for(year in unique(selected$Year)){
  tmp <- filter(selected, Year == year)[,c(6,7)]
  repNum = length(unique(tmp$UniqueID))
  if(repNum >= repMax){repMax = repNum}
  print (paste(repMax,repNum,sep=" "))
}
repNum = repMax

#Now insert the replicates with actual data in between the "NA" dummy columns
#which ensure that the final eLSA input file has an even number of replicates
#per year regardless of the variations in the actual number of sites (replicates)
#sampled per year.
eLSAtmp <- eLSAInput[,1]
j=1
k=1
nulCol <- data.frame(matrix(ncol=repNum*spotNum,nrow=length(unique(selected$FinalID))))
nulCol[,1] <- eLSAInput[,1]
for(year in unique(selected$Year)){
  tmp <- filter(selected, Year == year)
  rep = length(unique(tmp$UniqueID))
  for(i in 1:repNum){
    if(i <= rep){
      repLabel = paste(year,"DoneRep",i,sep="")
      j=j+1
      k=k+1
      eLSAtmp[,k] <- eLSAInput[,j]
    }
    else{
      repLabel = as.character(paste(year,"Rep",i,sep=""))
      k=k+1
      eLSAtmp[,k] <- "NA"
      print(paste(k,repLabel,sep=" "))
    }
  }
}

eLSAInput <- eLSAtmp

#If you want to determine the average parameter values across the subsetted data set.
library(matrixStats)
eLSAAverage <- eLSAInput
eLSAAverage <- as.data.frame(sapply(eLSAAverage,as.numeric))
eLSAAverage$FinalID <- eLSAInput$FinalID
eLSAAverage$mean <- rowMeans(eLSAAverage[,2:ncol(eLSAAverage)],na.rm=TRUE)
eLSAAverage$median <- rowMedians(as.matrix(eLSAAverage[,2:ncol(eLSAAverage)]),na.rm=TRUE)
eLSAAverage$STDEV <- rowSds(as.matrix(eLSAAverage[,2:ncol(eLSAAverage)]),na.rm=TRUE)

#If you want to output the parameter averages for a given geographic subset of data.
chemID <- unique(chemData$FinalID)
eLSAAverage <- subset(eLSAAverage,eLSAAverage$FinalID %in% chemID)[,c(-2:-(ncol(eLSAAverage)-3))]
colnames(eLSAAverage) <- c("FinalID",paste("mean",suffix),"median","STDEV")
write.table(eLSAAverage,paste("eLSAAverage",suffix,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)

#Aggreate all of the mean physical parameter averages into a single data frame
#for analysis.  Make sure you've already generated these subsetted files first.
suffixList = c("HD1K","MD1K","LD1K","HD5K","MD5K","LD5K","HDWS","MDWS","LDWS")
means <- data.frame(chemID)
colnames(means) <- c("FinalID")
for(suffix in suffixList){
  print(paste("eLSAAverage",suffix,".txt",sep=""))
  parameter <- read.delim(paste("eLSAAverage",suffix,".txt",sep=""),header=TRUE, sep="\t",as.is=T,check.names=FALSE)
  print(parameter[,-c(3:4)])
  parameter <- data.frame(parameter[,-c(3:4)])
  means <- join(means,parameter,by="FinalID")
}
means <- t(means)
colnames(means) <- means[1,]
means <- means[-1,]
means <- as.data.frame(means)

#Scale your average physical parameter dataframe to perform PCA.
library("factoextra")
means <- data.frame(sapply(means, function(x) as.numeric(as.character(x))))
log.means <- log(means)
row.names(log.means) <- suffixList
log.means[is.na(log.means)] <- 0
log.means <- log.means[,which(apply(log.means,2,var)!=0)]
means.pca <- prcomp(log.means,center=TRUE,scale.=TRUE)
var <- get_pca_var(means.pca)
var$coord[,1:3]
# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
#Plot the contributions of each principal component to the overall
#variation in the data set.
fviz_screeplot(means.pca)

# Variable correlation/coordinates
loadings <- means.pca$rotation
sdev <- means.pca$sdev

# Determine which physical variables drive most of the variation
# in the principal components of the system.
var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
var.coord[, 1:3]
var.cos2 <- var.coord^2
var.cos2[,1:3]
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[,1:2])
#Plot the chemical parameters most strongly correlated to a designated
#individual or set of principal components.  Store the correlations in a dataframe.
fviz_cos2(means.pca, choice = "var", axes = 1:2,sort.val="des")
cos2Plot <- fviz_cos2(means.pca, choice = "var", axes = 1:2,sort.val="des",top=100)
cos2Values <- cos2Plot$data

#Output dataframe for use in eLSA.
#Note that the the data needs to have at least two location replicates per time point
#and that the number of replicates per time point needs to be uniform.
#This may involve subsampling data depending on the variation in the number of replicates per time point.
#The first character in an eLSA formatted file needs to be #.
names(eLSAInput)[names(eLSAInput)=="FinalID"]<-"#FinalID"
write.table(eLSAInput,paste("eLSAInput",suffix,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#At this point insert columns with all NA values for each year in order to even
#out the number of replicates per year.

#Read in eLSA output.
#Compute network statistics of the likeliest association networks between taxa.
library(igraph)
library(network)
suffix <- "LDWS"
networkdata <- read.delim(paste("eLSAOutput",suffix,".txt",sep=""),header=TRUE, sep="\t",as.is=T,check.names=FALSE)
#Filter out association network data based on P scores, for the local similarity
#between two factors, with values less than 0.05.
networkdata <- filter(networkdata, P <= 0.01)
names(networkdata)[names(networkdata)=="LS"]<-"weight"
#Filter network data based on local similarity scores.
#networkdata <- subset(networkdata,networkdata$weight<0)
algaeID <- unique(algaeData$FinalID)
insectID <- unique(insectData$FinalID)
chemID <- unique(chemData$FinalID)

#Define a 'not in' function.
'%!in%' <- function(x,y)!('%in%'(x,y))
#Remove some subset of chemical and biological factors as nodes from the network.
#networkdata1 <- subset(networkdata,networkdata$X %!in% chemID & networkdata$Y %in% chemID)
#networkdata2 <- subset(networkdata,networkdata$Y %!in% chemID & networkdata$X %in% chemID)
#networkdata <- rbind(networkdata1,networkdata2)
networkdata <- subset(networkdata,networkdata$X %!in% chemID)
networkdata <- subset(networkdata,networkdata$Y %!in% chemID)

#Generate network graph and begin calculating network parameters.
networkgraph=graph.data.frame(networkdata,directed=FALSE)
l <- layout_with_fr(networkgraph)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
par(mfrow=c(1,1), mar=c(0,0,0,0))
#Option 1: plot just the structure of the network.
V(networkgraph)$label <- ""
plot(networkgraph,vertex.size=3)
#Option 2: plot the network and weight links by LS scores and node size by number of links.
#plot(networkgraph,rescale=F,layout=l*1.0,vertex.size=5*degree(networkgraph),edge.width=abs(E(networkgraph)$weight*10),edge.color=ifelse(E(networkgraph)$weight > 0, "blue","red"))
# Calculate number of groups and the modularity of the network.
cluster_edge_betweenness(networkgraph, weights=NULL,directed=FALSE)
# Calculate the average network path length
mean_distance(networkgraph,directed=FALSE)
# Calculate the clustering coefficient
transitivity(networkgraph,type="globalundirected",isolate="zero")
# Generate adjacency matrix of relative taxa abundance correlations
adj= as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
# Get the number of unique network edges
network.edgecount(adj)
# Get the number of nodes
network.size(adj)
# Get the network density.
network.density(adj)

#Plot a corresponding random graph given a number of edges and nodes.
RandGraph <- erdos.renyi.game(network.size(adj),network.edgecount(adj),type="gnm")
plot(RandGraph, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 3, xlab = "Random Network: G(N,L) model")

library(vcd)
library(MASS)
# Get degree distribution of network.
DDN <- degree(networkgraph)
# Fit a poisson distribution to the link distribution of the network
curveFit <- fitdistr(DDN,"exponential")
# Get the fit parameters for the distribution.
# Scale-free networks are exponential and random networks are Poisson distributions.
coef(curveFit)
# Get the log-likelihood for this fit
logLik(curveFit)
# Histogram of node degree distribution.
hist(degree(networkgraph, mode="all"), breaks=1:vcount(networkgraph)-1, main="Histogram of node degree")
# Find the largest clique within the network.
largest_cliques(as.undirected(networkgraph, mode= "collapse", edge.attr.comb=list(weight="sum", "ignore")))
# List nodes by their number of links.
degree(networkgraph)
