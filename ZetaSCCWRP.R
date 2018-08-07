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
#Taxa frequency table.
taxaFreq <- as.data.frame(table(GISBioData$FinalID))
colnames(taxaFreq) <- c("FinalID","Freq")
#Find the total number of taxa in the full data set.
taxaMax <- length(unique(GISBioData$FinalID))

#Generating a presence/absence matrix for California SCCWRP data.
selected <- GISBioData
selected <- subset(GISBioData,LU_2000_5K>=90)
#selected <- subset(selected,Year==2000)
#selected <- subset(selected,LU_2000_5K>=20)
#selected <- GISBioData
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

#Store presence/absence matrix of California SCCWRP data.
#write.table(data.SCCWRP,"SCCWRPCAPresenceAbsence.txt",quote=FALSE,sep="\t",row.names = FALSE)

#If the presence/absence matrix already exists load it here.
#data.SCCWRP <- read.table("SCCWRPCAPresenceAbsence.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Subset location data.
xy.SCCWRP <- selected[,c("Latitude","Longitude","UniqueID")]
xy.SCCWRP <- xy.SCCWRP[!duplicated(xy.SCCWRP[c("UniqueID")]),]
xy.SCCWRP <- xy.SCCWRP[,c("Latitude","Longitude")]

#Subset environmental factor data.
env.SCCWRP <- selected[,c("UniqueID","Watershed","Ag_2000_5K","CODE_21_2000_5K","URBAN_2000_5K","LU_2000_5K","Year")]
env.SCCWRP <- env.SCCWRP[!duplicated(env.SCCWRP[c("UniqueID")]),]
env.SCCWRP <- env.SCCWRP[,c("Watershed","Ag_2000_5K","CODE_21_2000_5K","URBAN_2000_5K","LU_2000_5K","Year")]
env.SCCWRP$Year <- as.numeric(env.SCCWRP$Year)
env.SCCWRP$Watershed <- as.factor(env.SCCWRP$Watershed)

#Plotting zeta diversity decay with distance.
Zeta.ddecay(xy.SCCWRP,data.SCCWRP,order=2,sam=nrow(data.SCCWRP),distance.type="ortho",plot=TRUE)

#Zeta diversity with respect to environmental variables.
Zeta.msgdm(data.SCCWRP,env.SCCWRP,xy=NULL,sam=1000,order=2)

#Computes zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
#using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
Zeta.decline.mc(data.SCCWRP,xy=NULL,orders=1:10,sam=1000)

#Clustering samples into groups based on maximizing overlap in taxa present.
library(ade4)
comDist <- dist.binary(data.SCCWRP,method=1,diag=FALSE,upper=FALSE)
comClust <- hclust(comDist,method="complete")
comGroups <- cutree(comClust,h=0.8)
comGroups <- as.data.frame(comGroups)
data.SCCWRP.subset <- merge(data.SCCWRP,comGroups,by="row.names",all.x=TRUE)
data.SCCWRP.subset <- subset(data.SCCWRP.subset,comGroups==1)
dev.off()
plot(as.dendrogram(comClust),horiz=TRUE)
rect.hclust(comClust,h=0.6)

library(proxy)
library(dendextend)
labels.SCCWRP <- selected[,c("UniqueID","Watershed","Ag_2000_5K","CODE_21_2000_5K","URBAN_2000_5K","LU_2000_5K","altitude","nTaxa","Latitude","Longitude")]
labels.SCCWRP <- labels.SCCWRP[!duplicated(labels.SCCWRP[c("UniqueID")]),]
labels.SCCWRP <- data.frame(labels.SCCWRP[,-1],row.names=labels.SCCWRP[,1])
colnames(labels.SCCWRP) <- c("Watershed","Ag_2000_5K","CODE_21_2000_5K","URBAN_2000_5K","LU_2000_5K","altitude","nTaxa","Latitude","Longitude")
data.SCCWRP.tagged <- merge(data.SCCWRP,labels.SCCWRP,by="row.names")

#MDS plotting of Jaccard similarity between samples.
mds <- cmdscale(dist(data.SCCWRP,method="binary"), eig=TRUE,k=2)
mdsPlot <- data.frame(mds$points, group=data.SCCWRP.tagged$Latitude)
#For continuous factors.
ggplot(mdsPlot)+geom_point(aes(x = X1,y = X2, color = group))+ scale_colour_gradientn(colours=rainbow(4))
#For discrete factors.
ggplot(mdsPlot)+geom_point(aes(x = X1,y = X2, color = group))+ theme(legend.position="none")
