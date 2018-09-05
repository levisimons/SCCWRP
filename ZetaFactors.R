#Script to analyze environmental factors, and their significances, in shaping trends in zeta diversity.
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
library(stats)
library(CINNA)
library(fitdistrplus)

#Check for co-occurrence frequencies by watershed in the SCCWRP data set.
setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage values.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Remove duplicate measures.
GISBioData <- GISBioData[!duplicated(GISBioData[,c("UniqueID","FinalID","Count")]),]
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Read in functional feeding group for each taxon.
#Abbreviations used in denoting functional feeding groups are as follows ( http://www.safit.org/Docs/CABW_std_taxonomic_effort.pdf ):
#P= predator MH= macrophyte herbivore OM= omnivore
#PA= parasite PH= piercer herbivore XY= xylophage (wood eater)
#CG= collector-gatherer SC= scraper
#CF= collector filterer SH= shredder 
FFG <- read.table("metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
# Filter data so only known functional feeding groups are kept.
FFG <- subset(FFG, FunctionalFeedingGroup != "")
# Generate functional feeding group data frame.
FFG <- FFG[,c("FinalID","LifeStageCode","FunctionalFeedingGroup")]
FFG <- subset(FFG,LifeStageCode=="L" | LifeStageCode=="X" | FinalID=="Hydrophilidae" | FinalID=="Hydraenidae")
#Merge in functional feeding groups into sample data.
GISBioData <- join(GISBioData,FFG[,c("FinalID","FunctionalFeedingGroup")],by=c("FinalID"))
FFGCounts <- na.omit(as.data.frame(unique(GISBioData$FunctionalFeedingGroup)))
colnames(FFGCounts) <- c("FunctionalFeedingGroups")
#How many samples per watershed?
groupNum=20
#Select watersheds with a large enough set of samples for analysis.
watersheds <- as.data.frame(table(SCCWRP$Watershed))
colnames(watersheds) <- c("Watershed","Samples")
GISBioData <- join(GISBioData,watersheds,by=c("Watershed"))
#Get samples only found in more heavily sampled watersheds.
GISBioDataLargeWS <- subset(GISBioData,Samples>=groupNum)
#Add the number of genera per sample
taxaBySample <- count(GISBioDataLargeWS,UniqueID)
colnames(taxaBySample) <- c("UniqueID","nTaxa")
GISBioDataLargeWS <- join(GISBioDataLargeWS,taxaBySample,by=c("UniqueID"))

#Initialize data frames to compute zeta diversity for genera and functional feeding groups
selected <- GISBioDataLargeWS
selected <- arrange(selected,Year,UniqueID)
#Get zeta diversity decay parameters for taxonomic diversity for the same set of samples within a given land use band.
eLSAInput <- as.data.frame(unique(selected$FinalID))
colnames(eLSAInput) <- c("FinalID")
eLSAInput <- as.data.frame(eLSAInput[order(as.character(eLSAInput$FinalID)),])
colnames(eLSAInput) <- c("FinalID")
taxa <- eLSAInput
eLSAInputRand <- eLSAInput
selected <- selected[order(selected$Year,selected$UniqueID,selected$FinalID),]
#Get zeta diversity decay parameters for functional feeding group diversity for the same set of samples within a given land use band.
FFGInput <- as.data.frame(unique(selected$FunctionalFeedingGroup))
colnames(FFGInput) <- c("FunctionalFeedingGroup")
FFGInput <- as.data.frame(FFGInput[order(as.character(FFGInput$FunctionalFeedingGroup)),])
colnames(FFGInput) <- c("FunctionalFeedingGroup")
FFGInput <- na.omit(FFGInput)
FFGrand <- FFGInput
FFgroups <- FFGInput

#Generate presence/absence matrices for genera and functional feeding groups by stream sample.
for(ID in unique(selected$UniqueID)){
  #Add the relative taxa abundances by column to a new dataframe.
  #The rows are the unique taxa in a given subset of data.
  tmp <- filter(selected, UniqueID == ID)[,c("FinalID","Measurement","UniqueID")]
  tmp <- as.data.frame(tmp[order(tmp$FinalID),])
  tmp <- tmp[-c(3)]
  colnames(tmp) <- c("FinalID",ID)
  tmp <- tmp %>% group_by(FinalID) %>% summarise_if(is.numeric,mean,na.rm=TRUE)
  tmp <- join(tmp,taxa,type="full",by=c("FinalID"))
  tmp <- as.data.frame(tmp[order(tmp$FinalID),])
  eLSAInput <- cbind(eLSAInput,tmp)
  eLSAInput <- eLSAInput[,!duplicated(colnames(eLSAInput))]
  #Compute functional feeding group diversity by sample and sample grouping.
  tmp2 <- filter(selected, UniqueID == ID)[,c("FunctionalFeedingGroup","Count","UniqueID")]
  tmp2 <- as.data.frame(tmp2[order(tmp2$FunctionalFeedingGroup),])
  tmp2 <- tmp2[-c(3)]
  colnames(tmp2) <-  c("FunctionalFeedingGroup",ID)
  tmp2 <- tmp2 %>% group_by(FunctionalFeedingGroup) %>% summarise_if(is.numeric,sum,na.rm=TRUE)
  tmp2 <- join(tmp2,FFgroups,type="full",by=c("FunctionalFeedingGroup"))
  tmp2 <- as.data.frame(tmp2[order(tmp2$FunctionalFeedingGroup),])
  tmp2 <- tmp2[!is.na(tmp2$FunctionalFeedingGroup),]
  FFGInput <- cbind(FFGInput,tmp2)
  FFGInput <-  FFGInput[,!duplicated(colnames(FFGInput))]
  #Randomly assign functional feeding groups to their sample counts to eventually test how
  #far from random their relative abundances are.
  tmp3 <- tmp2[sample(nrow(tmp2)),]
  tmp3$FunctionalFeedingGroup <- tmp2$FunctionalFeedingGroup
  colnames(tmp3) <-  c("FunctionalFeedingGroup",ID)
  FFGrand <- cbind(FFGrand,tmp3)
  FFGrand <-  FFGrand[,!duplicated(colnames(FFGrand))]
  #Randomly assign genera to their sample counts to eventually test how
  #far from random their relative abundances are.
  tmp4 <- tmp[sample(nrow(tmp)),]
  tmp4$FinalID <- tmp$FinalID
  colnames(tmp4) <- c("FinalID",ID)
  eLSAInputRand <- cbind(eLSAInputRand,tmp4)
  eLSAInputRand <- eLSAInputRand[,!duplicated(colnames(eLSAInputRand))]
}

#Generate a presence/absence dataframe for zeta diversity analysis of taxa.
#Rows for samples, columns for taxa IDs.
eLSAInput[is.na(eLSAInput)] <- 0
eLSANames <- eLSAInput$FinalID
data.SCCWRP <- as.data.frame(t(eLSAInput[,-c(1)]))
colnames(data.SCCWRP) <- eLSANames
data.SCCWRP[data.SCCWRP > 0] <- 1
#Generate a presence/absence dataframe for zeta diversity analysis of functional feeding groups.
#Rows for samples, columns for functional feeding group types.
FFGInput[is.na(FFGInput)] <- 0
FFGNames <- FFGInput$FunctionalFeedingGroup
ffg.SCCWRP <- as.data.frame(t(FFGInput[,-c(1)]))
colnames(ffg.SCCWRP) <- FFGNames
ffg.SCCWRP[ffg.SCCWRP > 0] <- 1
#Generate a presence/absence dataframe for zeta diversity analysis of randomly assigned functional feeding groups.
#Rows for samples, columns for functional feeding group types.
FFGrand[is.na(FFGrand)] <- 0
FFGrandNames <- FFGrand$FunctionalFeedingGroup
ffg.rand.SCCWRP <- as.data.frame(t(FFGrand[,-c(1)]))
colnames(ffg.rand.SCCWRP) <- FFGrandNames
ffg.rand.SCCWRP[ffg.rand.SCCWRP > 0] <- 1
#Generate a presence/absence dataframe for zeta diversity analysis of randomly assigned genera.
#Rows for samples, columns for genera.
eLSAInputRand[is.na(eLSAInputRand)] <- 0
eLSAInputRandRames <- eLSAInputRand$FinalID
data.rand.SCCWRP <- as.data.frame(t(eLSAInputRand[,-c(1)]))
colnames(data.rand.SCCWRP)<- eLSAInputRandRames
data.rand.SCCWRP[data.rand.SCCWRP > 0] <- 1

#Subset environmental factor data.
env.SCCWRP <- join(GISBioDataLargeWS,SCCWRP,by=c("UniqueID"))
env.SCCWRP <- env.SCCWRP[,c("UniqueID","Watershed","LU_2000_5K","altitude","Year")]
env.SCCWRP <- env.SCCWRP[!duplicated(env.SCCWRP[c("UniqueID")]),]
env.SCCWRP <- env.SCCWRP[,c("Watershed","LU_2000_5K","altitude","Year")]
env.SCCWRP$Watershed <- as.factor(env.SCCWRP$Watershed)
env.SCCWRP$Year <- as.numeric(env.SCCWRP$Year)
env.SCCWRP$altitude <- as.numeric(env.SCCWRP$altitude)
#Zeta diversity with respect to environmental variables.
zetaTest <- Zeta.msgdm(data.spec=data.SCCWRP,data.env=env.SCCWRP,xy=NULL,sam=nrow(env.SCCWRP),order=2,rescale=FALSE)
