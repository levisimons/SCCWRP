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
#values.  If this file is not yet generated then proceed with the following commands
#to generate it in the first place.
GISBiochemData <- read.table("GISBiochemData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter the data set so only the parameter of interest remains.
parameter="Nitrogen, Total"
GISBiochemData <- filter(GISBiochemData,FinalID==parameter)

#Make a temporary data frame to aggregate the average measurement per site if
#there are duplicated sites with data taken concurrently.
tmp <- aggregate(GISBiochemData$Measurement,list(GISBiochemData$UniqueID), mean)
colnames(tmp) <- c("UniqueID","Measurement")

#Remove duplicate site rows.
GISBiochemData <- GISBiochemData[!duplicated(GISBiochemData$UniqueID),]

#Remove the original measurement column.
GISBiochemData <- subset(GISBiochemData,select=-c(Measurement))

#Merge back in the averaged measurements column.
GISBiochemData <- merge(GISBiochemData,tmp,by=("UniqueID"))

#Read in California Stream Condition Index data
csciData <- read.csv("csci_scored_sites_tbl.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
csciData <- filter(csciData,CSCI!="NA")
names(csciData)[names(csciData)=="StationCode"]<-"SampleStationID"
names(csciData)[names(csciData)=="SAMPLEDATE"]<-"SampleDate"
#Make the date format uniform.
csciData$SampleDate <- as.Date(csciData$SampleDate,format="%m/%d/%y")
#Add in qualifier columns based on California Streams Condition Index.  The cutoff is 0.79.
csciData$CSCIQualifier <- ifelse(csciData$CSCI >= 0.79, "Healthy","Disturbed")
csciData$CSCIQualNum <- ifelse(csciData$CSCI >= 0.79, 1,0)
#Add in qualifier column based on total nitrogen concentrations.  The cutoff is 0.42mg/L.
#csciData$TNQualifier <- ifelse(csciData$`Total_Nitrogen (mg/L)`<=0.42,"Healthy","Disturbed")
#Add in qualifier column based on total phosphate concentrations.  The cutoff is 0.03mg/L.
#csciData$TPQualifier <- ifelse(csciData$`Total_Phosphorous (mg/L)`<=0.03,"Healthy","Disturbed")
#Parse out the year the CSCI was calculated by site.
csciData$Year <- year(csciData$SampleDate)
#Keep only the first replicate.  Most sites only have one replicate.
csciData <- filter(csciData,REPLICATE==1)
#Filter by year so that only CSCI values are found within the same time window
#as the merged site data set.
csciData <- filter(csciData,Year>=min(GISBiochemData$Year) & Year<= max(GISBiochemData$Year))

#Merge CSCI data onto site data.
GISBiochemData <- merge(GISBiochemData,csciData,by=c("SampleStationID"))

#Read in network statistics generated under different chemical parameter selection windows.
chemNetworks <- read.csv("SCCWRPNetworkAnalysisChemicalThresholdsSummary.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)

#Assign network parameters based on a particular parameter concentration.
threshold=0.42
GISBiochemData$Con_l_rL <-ifelse(GISBiochemData$Measurement>threshold,chemNetworks[chemNetworks$Prefix=='HTN','Con_l_rL'],chemNetworks[chemNetworks$Prefix=='LTN','Con_l_rL'])  
GISBiochemData$Con_l_rCl <-ifelse(GISBiochemData$Measurement>threshold,chemNetworks[chemNetworks$Prefix=='HTN','Con_l_rCl'],chemNetworks[chemNetworks$Prefix=='LTN','Con_l_rCl'])
GISBiochemData$Con_l_rM <-ifelse(GISBiochemData$Measurement>threshold,chemNetworks[chemNetworks$Prefix=='HTN','Con_l_rM'],chemNetworks[chemNetworks$Prefix=='LTN','Con_l_rM'])  
GISBiochemData$Cov_l_rL <-ifelse(GISBiochemData$Measurement>threshold,chemNetworks[chemNetworks$Prefix=='HTN','Cov_l_rL'],chemNetworks[chemNetworks$Prefix=='LTN','Cov_l_rL'])  
GISBiochemData$Cov_l_rCl <-ifelse(GISBiochemData$Measurement>threshold,chemNetworks[chemNetworks$Prefix=='HTN','Cov_l_rCl'],chemNetworks[chemNetworks$Prefix=='LTN','Cov_l_rCl'])  
GISBiochemData$Cov_l_rM <-ifelse(GISBiochemData$Measurement>threshold,chemNetworks[chemNetworks$Prefix=='HTN','Cov_l_rM'],chemNetworks[chemNetworks$Prefix=='LTN','Cov_l_rM'])  

testData <- select(GISBiochemData,CSCI,Con_l_rL,Con_l_rCl,Con_l_rM,Cov_l_rL,Cov_l_rCl,Cov_l_rM)
testData <- na.omit(testData)

#Logistic regression between network parameters and CSCI
library(PerformanceAnalytics)
library(aod)
library(glmm)
library(rcompanion)
model.vars <- names(testData)[2:7]
model.list <- lapply(model.vars, function(x){
  lm(substitute(CSCI ~ i, list(i=as.name(x))),data=testData)
})
lapply(model.list,summary)
