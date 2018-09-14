library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)
library(MASS)
library(zetadiv)
library(magrittr)

setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage values.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Remove duplicate measures.
GISBioData <- GISBioData[!duplicated(GISBioData[,c("UniqueID","FinalID","Count")]),]
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Merge in altitude.
GISBioData <- join(GISBioData,SCCWRP[,c("UniqueID","altitude")],by=c("UniqueID"))
#Add the number of genera per sample
taxaBySample <- count(GISBioData,UniqueID)
colnames(taxaBySample) <- c("UniqueID","nTaxa")
GISBioData <- join(GISBioData,taxaBySample,by=c("UniqueID"))
#Add the number of samples per watershed
WSBySample <- as.data.frame(table(SCCWRP$Watershed))
colnames(WSBySample) <- c("Watershed","NSamples")
GISBioData <- join(GISBioData,WSBySample,by=c("Watershed"))
#Read in watershed metadata
Watersheds <- read.table("CAWatersheds.tsv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
Watersheds$Watershed <- gsub('\\s+', '', Watersheds$Watershed)
#Merge in watershed area
GISBioData <- join(GISBioData,Watersheds[,c("Watershed","CatalogingUnitAreaSqMi")],by=c("Watershed"))

#Find watersheds with a larger enough set of samples for downstream analysis.
sampleMin <- 25
GISBioDataLWS <- subset(GISBioData,NSamples>=sampleMin)
set.seed(1)
zetaAnalysis <- data.frame()
for(WS in unique(GISBioDataLWS$Watershed)){
  GISBioDataLocal <- subset(GISBioDataLWS,Watershed==WS) #Subsample by watershed.
  GISBioDataLocal <- GISBioDataLocal[GISBioDataLocal$UniqueID %in% sample(unique(GISBioDataLocal$UniqueID),sampleMin),] #Subsample to a uniform sample number.
  metadata <- GISBioDataLocal[,c("UniqueID","LU_2000_5K","altitude","CSCI")]
  metadata <- metadata[!duplicated(metadata),] #Get unique environmental parameters per watershed set of samples.
  selected <- GISBioDataLocal
  #Get all unique taxa in statewide data set.
  uniqueTaxa <- as.data.frame(unique(selected$FinalID))
  colnames(uniqueTaxa) <- c("FinalID")
  uniqueTaxa <- arrange(uniqueTaxa,FinalID)
  #Create presence/absence matrix of taxa in samples.
  #Rows for sample ID and columns 
  PresenceAbsence <- uniqueTaxa
  for(ID in unique(selected$UniqueID)){
    sampleDF <- subset(GISBioData,UniqueID == ID)
    sampleDF <- sampleDF[,c("FinalID","Count")]
    tmp <- merge(sampleDF,uniqueTaxa,by=c("FinalID"),all=TRUE)
    colnames(tmp) <- c("FinalID",ID)
    PresenceAbsence <- cbind(PresenceAbsence,tmp[,c(2)])
    colnames(PresenceAbsence)[ncol(PresenceAbsence)] <- ID
  }
  #Generate a presence/absence dataframe for zeta diversity analysis of taxa.
  #Rows for samples, columns for taxa IDs.
  PresenceAbsence[is.na(PresenceAbsence)] <- 0
  PresenceAbsence[PresenceAbsence > 0] <- 1
  data.SCCWRP <- as.data.frame(t(PresenceAbsence[,-c(1)]))
  colnames(data.SCCWRP) <- uniqueTaxa$FinalID
  #Calculate zeta diversity decay within each watershed set of samples.
  zetaDecay <- Zeta.decline.ex(data.SCCWRP,orders=1:10)
  ExpIntercept <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
  ExpExp <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
  ExpAIC <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
  PLIntercept <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
  PLExp <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
  PLAIC <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
  row <- t(as.data.frame(c(WS,mean(metadata$CSCI),mean(metadata$LU_2000_5K),sd(metadata$LU_2000_5K),mean(metadata$altitude),sd(metadata$altitude),ExpIntercept,ExpExp,ExpAIC,PLIntercept,PLExp,PLAIC)))
  zetaAnalysis <- rbind(zetaAnalysis,row)
}
colnames(zetaAnalysis) <- c("Watershed","MeanCSCI","MeanLU","SDLU","MeanAltitude","SDAltitude","ExpIntercept","ExpExp","ExpAIC","PLIntercept","PLExp","PLAIC")
zetaAnalysis[,1:12] <- sapply(zetaAnalysis[,1:12],as.character)
zetaAnalysis[,2:12] <- sapply(zetaAnalysis[,2:12],as.numeric)
rownames(zetaAnalysis) <- 1:nrow(zetaAnalysis)
WSAreas <- GISBioDataLWS[,c("Watershed","CatalogingUnitAreaSqMi")]
WSAreas <- WSAreas[!duplicated(WSAreas),]
zetaAnalysis <- join(zetaAnalysis,WSAreas,by=c("Watershed"))




#Scrap###############

#Determine land use deciles for the full state data set.
LUdf <- GISBioData[,c("UniqueID","LU_2000_5K")]
LUdf <- LUdf[!duplicated(LUdf),]
LUquantile <- quantile(LUdf$LU_2000_5K,probs=seq(0,1,0.1))
#Determine altitude deciles for the full state data set.
Altitudedf <- GISBioData[,c("UniqueID","altitude")]
Altitudedf <- Altitudedf[!duplicated(Altitudedf),]
Altitudequantile <- quantile(Altitudedf$altitude,probs=seq(0,1,0.1))

selected <- GISBioData

#Get all unique taxa in statewide data set.
uniqueTaxa <- as.data.frame(unique(selected$FinalID))
colnames(uniqueTaxa) <- c("FinalID")
uniqueTaxa <- arrange(uniqueTaxa,FinalID)

#Create presence/absence matrix of taxa in samples.
#Rows for sample ID and columns 
PresenceAbsence <- uniqueTaxa
for(ID in unique(selected$UniqueID)){
  sampleDF <- subset(GISBioData,UniqueID == ID)
  sampleDF <- sampleDF[,c("FinalID","Count")]
  tmp <- merge(sampleDF,uniqueTaxa,by=c("FinalID"),all=TRUE)
  colnames(tmp) <- c("FinalID",ID)
  PresenceAbsence <- cbind(PresenceAbsence,tmp[,c(2)])
  colnames(PresenceAbsence)[ncol(PresenceAbsence)] <- ID
}
#Generate a presence/absence dataframe for zeta diversity analysis of taxa.
#Rows for samples, columns for taxa IDs.
PresenceAbsence[is.na(PresenceAbsence)] <- 0
PresenceAbsence[PresenceAbsence > 0] <- 1
data.SCCWRP <- as.data.frame(t(PresenceAbsence[,-c(1)]))
colnames(data.SCCWRP) <- uniqueTaxa$FinalID

#Subset environmental factor data.
env.SCCWRP <- join(GISBioData,SCCWRP,by=c("UniqueID"))
env.SCCWRP <- env.SCCWRP[,c("UniqueID","Watershed","LU_2000_5K","altitude","Year")]
env.SCCWRP <- env.SCCWRP[!duplicated(env.SCCWRP[c("UniqueID")]),]
env.SCCWRP <- env.SCCWRP[,c("Watershed","LU_2000_5K","altitude","Year")]
env.SCCWRP <- env.SCCWRP[,c("Watershed","LU_2000_5K","altitude")]
env.SCCWRP$LU_2000_5K <- as.numeric(env.SCCWRP$LU_2000_5K)
env.SCCWRP$Watershed <- as.factor(env.SCCWRP$Watershed)
#env.SCCWRP$Year <- as.factor(env.SCCWRP$Year)
env.SCCWRP$altitude <- as.numeric(env.SCCWRP$altitude)

#Subset location data.
xy.SCCWRP <- SCCWRP[,c("Latitude","Longitude","UniqueID")]
xy.SCCWRP <- xy.SCCWRP[!duplicated(xy.SCCWRP[c("UniqueID")]),]
xy.SCCWRP <- xy.SCCWRP[,c("Latitude","Longitude")]

#Zeta diversity with respect to environmental variables.
zetaTest <- Zeta.msgdm(data.spec=data.SCCWRP,data.env=env.SCCWRP,xy=xy.SCCWRP,reg.type="glm",sam=nrow(env.SCCWRP),order=2,rescale=FALSE)
#Zeta.varpart returns a data frame with one column containing the variation explained by each component
#a (the variation explained by distance alone), b (the variation explained by either distance or the environment),
#c (the variation explained by the environment alone) and d (the unexplained variation).
zetaVar <- Zeta.varpart(zetaTest)

#Plotting zeta diversity decay with distance.
zetaDist <- Zeta.ddecays(xy.SCCWRP,data.SCCWRP,orders=2:10,normalize="Jaccard",sam=nrow(data.SCCWRP),distance.type="ortho",plot=TRUE)
zetaDist

#Computes the expectation of zeta diversity, the number of species shared by multiple assemblages,
#for a specific order (number of assemblages or sites) using a formula based on the occupancy of the species.
Zeta.order.ex(data.SCCWRP, order = 4, sd.correct = TRUE, rescale = FALSE,empty.row = "empty")

#Zeta diversity for the entire state
zetaState <- Zeta.decline.ex(data.SCCWRP,orders=1:10)

#Determine land use deciles for the full state data set.
LUdf <- SCCWRP[,c("UniqueID","LU","altitude")]
colnames(LUdf) <- c("UniqueID","LU_2000_5K","altitude")
LUdf <- subset(LUdf,altitude <= median(LUdf$altitude))
LUdf <- LUdf[!duplicated(LUdf),]
LUdf <- arrange(LUdf,LU_2000_5K)
#LUquantile <- quantile(LUdf$LU_2000_5K,probs=seq(0,1,0.1))
#Determine altitude deciles for the full state data set.
Altitudedf <- SCCWRP[,c("UniqueID","LU","altitude")]
colnames(Altitudedf) <- c("UniqueID","LU_2000_5K","altitude")
Altitudedf <- Altitudedf[!duplicated(Altitudedf),]
Altitudedf <- subset(Altitudedf,LU_2000_5K <= median(Altitudedf$LU_2000_5K))
Altitudedf <- arrange(Altitudedf,altitude)
#Altitudequantile <- quantile(Altitudedf$altitude,probs=seq(0,1,0.1))

#Determine zeta diversity for particular subsets of sample data by land use quantiles.
zetaAnalysis <- data.frame()
#Divide sample metadata into quantiles
Rowquantile <- quantile(1:nrow(LUdf),probs=seq(0,1,0.1))
NSamples <- 1000
for(i in 1:length(Rowquantile)){
  if(i>1){
    LUSubset <- LUdf[Rowquantile[i-1]:Rowquantile[i],]
    MidLU <- mean(LUSubset$LU_2000_5K)
    localPA <- data.SCCWRP[unique(LUSubset$UniqueID),]
    dat <- data.frame()
    #Compute zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
    #using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
    zetaDecay <- Zeta.decline.mc(localPA,xy=NULL,orders=1:10,sam=NSamples)
    dat[1,1] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    dat[1,2] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    dat[1,3] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    dat[1,4] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
    dat[1,5] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    dat[1,6] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    dat[1,7] <- MidLU
    zetaAnalysis <- rbind(zetaAnalysis,dat)
    print(paste(MidLU,length(unique(LUSubset$UniqueID))))
  }
}
colnames(zetaAnalysis) <- c("zetaExpIntercept","zetaExpExponent","zetaExpAIC","zetaPLIntercept","zetaPLExponent","zetaPLAIC","MidLU")

#Determine zeta diversity for particular subsets of sample data by altitude quantiles.
zetaAnalysis <- data.frame()
#Divide sample metadata into quantiles
Rowquantile <- quantile(1:nrow(Altitudedf),probs=seq(0,1,0.1))
NSamples <- 1000
for(i in 1:length(LUquantile)){
  if(i>1){
    AltitudeSubset <- Altitudedf[Rowquantile[i-1]:Rowquantile[i],]
    MidAltitude <- mean(AltitudeSubset$altitude)
    localPA <- data.SCCWRP[unique(AltitudeSubset$UniqueID),]
    dat <- data.frame()
    #Compute zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
    #using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
    zetaDecay <- Zeta.decline.mc(localPA,xy=NULL,orders=1:10,sam=NSamples)
    dat[1,1] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    dat[1,2] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    dat[1,3] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    dat[1,4] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
    dat[1,5] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    dat[1,6] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    dat[1,7] <- MidAltitude
    zetaAnalysis <- rbind(zetaAnalysis,dat)
    print(paste(MidAltitude,length(unique(AltitudeSubset$UniqueID))))
  }
}
