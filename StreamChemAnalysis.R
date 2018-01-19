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

#Read in geospatial data.
GISDataRAW <- read.table("GIS_dnaSites.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset columns of interest.
GISData <- GISDataRAW[,-c(1,3:5,8:10,15)]
names(GISData)[names(GISData)=="StationCode"]<-"SampleStationID"
names(GISData)[names(GISData)=="New_Lat"]<-"Latitude"
names(GISData)[names(GISData)=="New_Long"]<-"Longitude"

#Merge geospatial data with biological observations.
GISChemData <- join(chemData,GISData,by="SampleStationID")
#Sort merged data set by year then measurement name.
GISChemData <- as.data.frame(GISChemData[order(as.numeric(GISChemData$Year),as.character(GISChemData$FinalID)),])
GISChemData <- as.data.frame(GISChemData[order(as.character(GISChemData$UniqueID)),])

#Filter out erroneous negative values for physical parameter data.
chemID <- as.data.frame(unique(chemData$FinalID))
colnames(chemID) <- c("FinalID")
GISChemData <- filter(GISChemData,Measurement>=0)

#Filter out duplicate rows.
GISChemData <- GISChemData[!duplicated(GISChemData),]

#Calculate land usage index based on 1K, 5K, and catchment zone values.
#Use land usage data from 2011.
GISChemData$LU_2011_1K <- with(GISChemData,Ag_2011_1K+CODE_21_2011_1K+URBAN_2011_1K)
GISChemData$LU_2011_5K <- with(GISChemData,Ag_2011_5K+CODE_21_2011_5K+URBAN_2011_5K)
GISChemData$LU_2011_WS <- with(GISChemData,Ag_2011_WS+CODE_21_2011_WS+URBAN_2011_WS)

#Filter out data without a land usage index.
GISChemData <- filter(GISChemData,LU_2011_1K!="NA")
GISChemData <- filter(GISChemData,LU_2011_5K!="NA")
GISChemData <- filter(GISChemData,LU_2011_WS!="NA")

#Subset site data based on land usage index within 1K catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISChemDataLD1K <- GISChemData[which(GISChemData$LU_2011_1K < 5),]
GISChemDataLD1K$LUCategory <- "LD1K"
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISChemDataMD1K <- GISChemData[which(GISChemData$LU_2011_1K < 15 & GISChemData$LU_2011_1K >= 5),]
GISChemDataMD1K$LUCategory <- "MD1K"
#HD = low disturbance.  Land usage index is greater than 15%.
GISChemDataHD1K <- GISChemData[which(GISChemData$LU_2011_1K >= 15),]
GISChemDataHD1K$LUCategory <- "HD1K"

#Subset site data based on land usage index within 5K catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISChemDataLD5K <- GISChemData[which(GISChemData$LU_2011_5K < 5),]
GISChemDataLD5K$LUCategory <- "LD5K"
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISChemDataMD5K <- GISChemData[which(GISChemData$LU_2011_5K < 15 & GISChemData$LU_2011_5K >= 5),]
GISChemDataMD5K$LUCategory <- "MD5K"
#HD = low disturbance.  Land usage index is greater than 15%.
GISChemDataHD5K <- GISChemData[which(GISChemData$LU_2011_5K >= 15),]
GISChemDataHD5K$LUCategory <- "HD5K"

#Subset site data based on land usage index within the full water drainage catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISChemDataLDWS <- GISChemData[which(GISChemData$LU_2011_WS < 5),]
GISChemDataLDWS$LUCategory <- "LDWS"
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISChemDataMDWS <- GISChemData[which(GISChemData$LU_2011_WS < 15 & GISChemData$LU_2011_WS >= 5),]
GISChemDataMDWS$LUCategory <- "MDWS"
#HD = low disturbance.  Land usage index is greater than 15%.
GISChemDataHDWS <- GISChemData[which(GISChemData$LU_2011_WS >= 15),]
GISChemDataHDWS$LUCategory <- "HDWS"

#Merge land usage subsets back together for later analytical tools.
GISChemData <- do.call("rbind",list(GISChemDataLD1K,GISChemDataMD1K,GISChemDataHD1K,GISChemDataLD5K,GISChemDataMD5K,GISChemDataHD5K,GISChemDataLDWS,GISChemDataMDWS,GISChemDataHDWS))

#Read in the covariant component of the networks generated via land usage subsetting
covariantNetwork <- read.csv("SCCWRPNetworkAnalysisLandP01Covariant.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
covariantNetwork[covariantNetwork=="#DIV/0!" | covariantNetwork=="#NUM!"] <- "NA"
covariantNetwork <- cbind(covariantNetwork[,1:3],data.frame(sapply(covariantNetwork[,4:ncol(covariantNetwork)],function(x) as.numeric(as.character(x)))))

#Read in the covariant component of the networks generated via land usage subsetting
contravariantNetwork <- read.csv("SCCWRPNetworkAnalysisLandP01Contravariant.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
contravariantNetwork[contravariantNetwork=="#DIV/0!" | contravariantNetwork=="#NUM!"] <- "NA"
contravariantNetwork <- cbind(contravariantNetwork[,1:3],data.frame(sapply(contravariantNetwork[,4:ncol(contravariantNetwork)],function(x) as.numeric(as.character(x)))))

#Read in California Streams Condition Index data
CSCIData <- read.csv("csci_scored_sites_tbl.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
names(CSCIData)[names(CSCIData)=="StationCode"]<-"SampleStationID"
CSCIData <- filter(CSCIData,CSCIData$CSCI!="NA")

#Merge in CSCI data and network statistics data
GISChemCSCIData <- join(GISChemData,CSCIData,by="SampleStationID")
GISChemCSCIDataCovariant <- join(GISChemCSCIData,covariantNetwork,by="LUCategory")
GISChemCSCIDataCovariant <- filter(GISChemCSCIDataCovariant,GISChemCSCIDataCovariant$CSCI!="NA" & GISChemCSCIDataCovariant$l_rL!="NA")
GISChemCSCIDataContravariant <- join(GISChemCSCIData,contravariantNetwork,by="LUCategory")
GISChemCSCIDataContravariant <- filter(GISChemCSCIDataContravariant,GISChemCSCIDataContravariant$CSCI!="NA" & GISChemCSCIDataContravariant$l_rL!="NA")

#Check for correlations between network parameters and the CSCI
v <- cbind(as.numeric(GISChemCSCIDataCovariant$CSCI),as.numeric(GISChemCSCIDataCovariant$l_rL))
v <- na.omit(v)
landCor <- round(cor(v[,1],v[,2],method="pearson"),4)
dev.off()
plot(v[,1],v[,2],main=paste("Covariant l_rL vs. CSCI","\n","Pearson r = ",landCor),ylab="Covariant l_rL",xlab="CSCI")
abline(lm(v[,2]~v[,1]),col="red")
v <- cbind(as.numeric(GISChemCSCIDataContravariant$CSCI),as.numeric(GISChemCSCIDataContravariant$l_rL))
v <- na.omit(v)
cor(v[,1],v[,2],method="pearson")
landCor <- round(cor(v[,1],v[,2],method="pearson"),4)
dev.off()
plot(v[,1],v[,2],main=paste("Contravariant l_rL vs. CSCI","\n","Pearson r = ",landCor),ylab="Contravariant l_rL",xlab="CSCI")
abline(lm(v[,2]~v[,1]),col="red")

#Create a unified dataframe of the averaged chemical parameters for each sample
#in order to perform principal component analysis on water chemistry data.
tmp <- as.data.frame(chemID)
colnames(tmp) <- c("FinalID")
for(site in unique(GISChemData$UniqueID)){
  GISChemDataSite <- GISChemData[GISChemData$UniqueID == site,][,3:4]
  GISChemDataSite <- aggregate(GISChemDataSite[,2],list(GISChemDataSite$FinalID),mean)
  colnames(GISChemDataSite) <- c("FinalID",paste(site,"Measurement"))
  tmp <- join(tmp,GISChemDataSite)
}

means <- as.data.frame(t(tmp))
chemNames <- t(tmp$FinalID)
means <- means[-1,]
colnames(means) <- chemNames

#Scale your average physical parameter dataframe to perform PCA.
library("factoextra")
library(ggfortify)
means <- data.frame(sapply(means, function(x) as.numeric(as.character(x))))
log.means <- log(means)
#row.names(log.means) <- suffixList
log.means[is.na(log.means)] <- 0
log.means <- log.means[,which(apply(log.means,2,var)!=0)]
means.pca <- prcomp(log.means,center=TRUE,scale.=TRUE)
autoplot(means.pca)
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
(var.contrib[,1:2])
fviz_cos2(means.pca, choice = "var", axes = 1)
chemCor<-fviz_cos2(means.pca, choice = "var", axes = 1:2)$data

#Determine how correlated each chemical parameter is with land usage indices.
corRow <- data.frame(matrix(nrow=1,ncol=4))
colnames(corRow) <- c("Parameter","Correlation with LU_2011_1K","Correlation with LU_2011_5K","Correlation with LU_2011_WS")
chemLand <- corRow
for(uniqueChem in unique(GISChemData$FinalID)){
  tmp2 <- filter(GISChemData,GISChemData$FinalID==uniqueChem & GISChemData$Measurement!="NA" & GISChemData$LU_2011_1K!="NA")
  corRow$Parameter <- uniqueChem
  corRow$`Correlation with LU_2011_1K` <- cor(tmp2$Measurement,tmp2$LU_2011_1K,method="spearman")
  corRow$`Correlation with LU_2011_5K` <- cor(tmp2$Measurement,tmp2$LU_2011_5K,method="spearman")
  corRow$`Correlation with LU_2011_WS` <- cor(tmp2$Measurement,tmp2$LU_2011_WS,method="spearman")
  chemLand <- rbind(chemLand,corRow)
  print(corRow[1,])
}
chemLand <- na.omit(chemLand)
