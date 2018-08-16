library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)
library(hierDiversity)
library(zetadiv)


#This script focuses on zeta diversity patterns across the SCCWRP California streams data set.
setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage
#values.  If this file is not yet generated then proceed with the following commands
#to generate it in the first place.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Remove duplicate measures.
GISBioData <- GISBioData[!duplicated(GISBioData[,c("UniqueID","FinalID","Count")]),]
#Order data by LU_2000_5K.
GISBioData <- arrange(GISBioData,LU_2000_5K)
#Add taxa counts by sample.
tmp <- data.frame(table(GISBioData$UniqueID))
colnames(tmp) <- c("UniqueID","nTaxa")
GISBioData <- join(GISBioData,tmp,by=c("UniqueID"))
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Merge in sample altitude.
GISBioData <- join(GISBioData,SCCWRP[,c("UniqueID","altitude")],by=c("UniqueID"))
#Get samples per watershed.
watersheds <- as.data.frame((table(SCCWRP$Watershed)))
colnames(watersheds) <- c("Watershed","Samples")
GISBioData <- join(GISBioData,watersheds,by=c("Watershed"))
#Taxa frequency table.
taxaFreq <- as.data.frame(table(GISBioData$FinalID))
colnames(taxaFreq) <- c("FinalID","Freq")
#Find the total number of taxa in the full data set.
taxaMax <- length(unique(GISBioData$FinalID))

#Generating a presence/absence matrix for California SCCWRP data.
selected <- GISBioData
selected <- arrange(selected,Year,UniqueID)
eLSAInput <- as.data.frame(unique(selected$FinalID))
colnames(eLSAInput)<-c("FinalID")
eLSAInput <- as.data.frame(eLSAInput[order(as.character(eLSAInput$FinalID)),])
colnames(eLSAInput)<-c("FinalID")
taxa <- eLSAInput

#Add the relative taxa abundances by column to a new dataframe.
#The rows are the unique taxa in a given subset of data.
selected <- selected[order(selected$Year,selected$UniqueID,selected$FinalID),]
i=0
for(ID in unique(selected$UniqueID)){
  tmp <- filter(selected, UniqueID == ID)[,c("FinalID","Measurement","UniqueID")]
  tmp <- as.data.frame(tmp[order(tmp$FinalID),])
  tmp <- tmp[-c(3)]
  colnames(tmp)<-c("FinalID",ID)
  tmp <- tmp %>% group_by(FinalID) %>% summarise_if(is.numeric,mean,na.rm=TRUE)
  tmp <- join(tmp,taxa,type="full",by=c("FinalID"))
  tmp <- as.data.frame(tmp[order(tmp$FinalID),])
  eLSAInput <- cbind(eLSAInput,tmp)
  eLSAInput <- eLSAInput[,!duplicated(colnames(eLSAInput))]
  #eLSAInput <- eLSAInput[,-c(1)]
  #eLSAInput$FinalID=as.character(eLSAInput$FinalID)
  #eLSAInput <- eLSAInput %>% group_by(FinalID) %>% summarise_if(is.numeric,mean,na.rm=TRUE)
  i=i+1
  print(i)
}

#Generate a presence/absence dataframe for zeta diversity analysis.
#Rows for samples, columns for taxa IDs.
eLSAInput[is.na(eLSAInput)] <- 0
eLSANames <- eLSAInput$FinalID
data.SCCWRP <- as.data.frame(t(eLSAInput[,-c(1)]))
colnames(data.SCCWRP) <- eLSANames
data.SCCWRP[data.SCCWRP > 0] <- 1


#Get the samples per watershed for watersheds with at least a certain number of samples.
selected <- subset(GISBioData,Samples>=100)
LUhigh = 100
LULow = 15
data.SCCWRP.LUband <- data.frame()
for(WS in unique(selected$Watershed)){
  WSSubset <- subset(selected,Watershed==WS & LU_2000_5K <= LUhigh & LU_2000_5K > LULow)
  if(length(unique(WSSubset$UniqueID))>20){
    data.SCCWRP.subset <- subset(data.SCCWRP, rownames(data.SCCWRP) %in% sample(unique(WSSubset$UniqueID),20))
    data.SCCWRP.LUband <- rbind(data.SCCWRP.LUband,data.SCCWRP.subset)
    #Computes zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
    #using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
    zetaDecay <- Zeta.decline.mc(data.SCCWRP.subset,xy=NULL,orders=1:10,sam=1000)
    print(paste(WS,length(unique(WSSubset$UniqueID)),zetaDecay$zeta.pl$coefficients[2]))
  }
}
zetaDecay <- Zeta.decline.mc(data.SCCWRP.LUband,xy=NULL,orders=1:10,sam=1000)
print(paste(LULow,LUhigh,zetaDecay$zeta.pl$coefficients[2]))
#Store presence/absence matrix of California SCCWRP data.
#write.table(data.SCCWRP,"SCCWRPCAPresenceAbsence.txt",quote=FALSE,sep="\t",row.names = FALSE)

#If the presence/absence matrix already exists load it here.
#data.SCCWRP <- read.table("SCCWRPCAPresenceAbsence.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Subset location data.
xy.SCCWRP <- selected[,c("Latitude","Longitude","UniqueID")]
xy.SCCWRP <- xy.SCCWRP[!duplicated(xy.SCCWRP[c("UniqueID")]),]
xy.SCCWRP <- xy.SCCWRP[,c("Latitude","Longitude")]

#Subset environmental factor data.
env.SCCWRP <- selected[,c("UniqueID","Watershed","Ag_2000_5K","CODE_21_2000_5K","URBAN_2000_5K","LU_2000_5K","altitude","Year")]
env.SCCWRP <- env.SCCWRP[!duplicated(env.SCCWRP[c("UniqueID")]),]
env.SCCWRP <- env.SCCWRP[,c("Watershed","Ag_2000_5K","CODE_21_2000_5K","URBAN_2000_5K","LU_2000_5K","altitude","Year")]
env.SCCWRP$Year <- as.numeric(env.SCCWRP$Year)
env.SCCWRP$altitude <- as.numeric(env.SCCWRP$altitude)
#env.SCCWRP$Latitude <- as.numeric(env.SCCWRP$Latitude)
#env.SCCWRP$Longitude <- as.numeric(env.SCCWRP$Longitude)
#env.SCCWRP$nTaxa <- as.numeric(env.SCCWRP$nTaxa)
env.SCCWRP$Watershed <- as.factor(env.SCCWRP$Watershed)

#Plotting zeta diversity decay with distance.
Zeta.ddecay(xy.SCCWRP,data.SCCWRP,order=2,sam=nrow(data.SCCWRP),distance.type="ortho",plot=TRUE)

#Zeta diversity with respect to environmental variables.
zetaTest <- Zeta.msgdm(data.SCCWRP,env.SCCWRP,xy=NULL,sam=nrow(env.SCCWRP),order=2,rescale=FALSE)

#Computes zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
#using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
zetaDecay <- Zeta.decline.mc(data.SCCWRP,xy=NULL,orders=1:10,sam=1000)

#Analyzing trends relating to taxonomic overlap between samples.
library(ade4)
library(proxy)
library(dendextend)
labels.SCCWRP <- selected[,c("UniqueID","Watershed","Ag_2000_5K","CODE_21_2000_5K","URBAN_2000_5K","LU_2000_5K","altitude","nTaxa","Latitude","Longitude","Year","CSCI")]
labels.SCCWRP <- labels.SCCWRP[!duplicated(labels.SCCWRP[c("UniqueID")]),]
labels.SCCWRP <- data.frame(labels.SCCWRP[,-1],row.names=labels.SCCWRP[,1])
colnames(labels.SCCWRP) <- c("Watershed","Ag_2000_5K","CODE_21_2000_5K","URBAN_2000_5K","LU_2000_5K","altitude","nTaxa","Latitude","Longitude","Year","CSCI")
data.SCCWRP.tagged <- merge(data.SCCWRP,labels.SCCWRP,by="row.names")
#MDS plotting of Jaccard similarity between samples.
mds <- cmdscale(dist(data.SCCWRP,method="binary"), eig=TRUE,k=2)
mdsPlot <- data.frame(mds$points, group=data.SCCWRP.tagged$Watershed)
#For continuous factors.
ggplot(mdsPlot)+geom_point(aes(x = X1,y = X2, color = group))+ scale_colour_gradientn(colours=rainbow(4))
#For discrete factors.
ggplot(mdsPlot)+geom_point(aes(x = X1,y = X2, color = group))+ theme(legend.position="right")
