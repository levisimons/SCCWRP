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

#How many samples per watershed?
groupNum=20

#Look through watersheds which are sufficiently sampled and do the following:
#Divide the entire data set, composed of heavily sampled watersheds, into quantiles based on land use.
#Calculate zeta diversity of each land use band.
#Determine the histogram of the number of watersheds each co-occurrence occurs in.
#Calculate the parameters of the gamma distribution fit to this histogram.
#Select watersheds with a large enough set of samples for analysis.
zetaAnalysis <- data.frame()
watersheds <- subset(as.data.frame(table(SCCWRP$Watershed)),Freq>=groupNum)
colnames(watersheds) <- c("Watershed","Samples")
#Get samples only found in more heavily sampled watersheds.
GISBioDataLargeWS <- subset(GISBioData,Watershed %in% watersheds$Watershed)
LUquantile <- quantile(GISBioDataLargeWS$LU_2000_5K,probs=seq(0,1,0.1))#To get land use quantiles.
for(i in 1:length(LUquantile)){
  watersheds <- subset(as.data.frame(table(SCCWRP$Watershed)),Freq>=groupNum)
  colnames(watersheds) <- c("Watershed","Samples")
  #Get samples only found in more heavily sampled watersheds.
  GISBioDataLargeWS <- subset(GISBioData,Watershed %in% watersheds$Watershed)
  LULow <- as.numeric(LUquantile[i])
  if(i<length(LUquantile)){
    LUHigh <- as.numeric(LUquantile[i+1])
  }
  if(i==length(LUquantile)){
    LUHigh == 100
  }
  if(LULow == 0 & LUHigh == 0){
    GISBioDataLargeWS <- subset(GISBioDataLargeWS,LU_2000_5K==LULow) #Subset samples by aggregated land use.
    #print(paste(i-1,i,LULow,LUHigh,length(unique(GISBioDataLargeWS$UniqueID))))
  }
  if(LULow != LUHigh){
    GISBioDataLargeWS <- subset(GISBioDataLargeWS,LU_2000_5K>=LULow & LU_2000_5K < LUHigh) #Subset samples by aggregated land use.
    #print(paste(i-1,i,LULow,LUHigh,length(unique(GISBioDataLargeWS$UniqueID))))
  }
  if(i < length(LUquantile)){
    selected <- GISBioDataLargeWS
    selected <- arrange(selected,Year,UniqueID)
    #Get zeta diversity decay parameters for taxonomic diversity for the same set of samples within a given land use band.
    eLSAInput <- as.data.frame(unique(selected$FinalID))
    colnames(eLSAInput) <- c("FinalID")
    eLSAInput <- as.data.frame(eLSAInput[order(as.character(eLSAInput$FinalID)),])
    colnames(eLSAInput) <- c("FinalID")
    taxa <- eLSAInput
    selected <- selected[order(selected$Year,selected$UniqueID,selected$FinalID),]
    #Get zeta diversity decay parameters for functional feeding group diversity for the same set of samples within a given land use band.
    FFGInput <- as.data.frame(unique(selected$FunctionalFeedingGroup))
    colnames(FFGInput) <- c("FunctionalFeedingGroup")
    FFGInput <- as.data.frame(FFGInput[order(as.character(FFGInput$FunctionalFeedingGroup)),])
    colnames(FFGInput) <- c("FunctionalFeedingGroup")
    FFGInput <- na.omit(FFGInput)
    FFgroups <- FFGInput
    i=0
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
    
    dat <- data.frame()
    #Computes zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
    #using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
    zetaDecay <- Zeta.decline.mc(data.SCCWRP,xy=NULL,orders=1:10,sam=1000)
    dat[1,1] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    dat[1,2] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    dat[1,3] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    dat[1,4] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
    dat[1,5] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    dat[1,6] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    #Computes zeta diversity, the number of functional feeding groups shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
    #using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
    zetaDecay <- Zeta.decline.mc(ffg.SCCWRP,xy=NULL,orders=1:10,sam=1000)
    dat[1,7] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    dat[1,8] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    dat[1,9] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    dat[1,10] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
    dat[1,11] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    dat[1,12] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    
    #Get the frequency of pairs of taxa showing up by watershed.
    #For example, a pair showing up three times in one watershed, and eight times in a second watershed,
    #will be counted as having shown up in two unique watersheds.
    CAMatch <- data.frame()
    #Which watersheds are heavily sampled and meet a land use criterion?
    watersheds <- subset(watersheds,Watershed %in% GISBioDataLargeWS$Watershed)
    for(WS in watersheds$Watershed){
      WSSubset <- subset(GISBioDataLargeWS,Watershed==WS)
      WSMatch <- data.frame()
      for(ID in unique(WSSubset$UniqueID)){
        WSSample <- subset(WSSubset,UniqueID==ID)
        #Subset samples by aggregated land use.
        if(length(unique(WSSample$FinalID))>2){
          taxaMatch <- as.data.frame(t(combn(unique(WSSample$FinalID),2)))
          WSMatch <- rbind(WSMatch,taxaMatch)
          WSMatch <- WSMatch[!duplicated(WSMatch[,c("V1","V2")]),]
          print(paste(WS,ID,length(unique(WSSample$FinalID))))
        }
      }
      WSName <- data.frame(matrix(nrow=nrow(WSMatch),ncol=1))
      colnames(WSName) <- c("Watershed")
      WSName$Watershed <- WS
      WSMatch <- cbind(WSMatch,WSName)
      CAMatch <- rbind(CAMatch,WSMatch)
    }
    CAMatch <- ddply(CAMatch, .(CAMatch$V1,CAMatch$V2),nrow)
    colnames(CAMatch) <- c("V1","V2","NumWS")
    hist(CAMatch$NumWS,xlim=c(0,63),ylim=c(0,40000))
    histDecay <- fitdist(CAMatch$NumWS,"gamma",method="mle")
    dat[1,13] <- as.numeric(histDecay$estimate[1]) #Gamma distribution histogram fit shape parameter.
    dat[1,14] <- as.numeric(histDecay$sd[1]) #Gamma distribution histogram fit shape parameter standard error.
    dat[1,15] <- as.numeric(histDecay$estimate[2])#Gamma distribution histogram fit rate parameter.
    dat[1,16] <- as.numeric(histDecay$sd[2]) #Gamma distribution histogram fit rate parameter standard error.
    dat[1,17] <- LULow
    dat[1,18] <- LUHigh
    print(dat)
  }
  zetaAnalysis <- rbind(zetaAnalysis,dat)
}
colnames(zetaAnalysis) <- c("zetaExpIntercept","zetaExpExponent","zetaExpAIC","zetaPLIntercept","zetaPLExponent","zetaPLAIC","zetaFFGExpIntercept","zetaFFGExpExponent","zetaFFGExpAIC","zetaFFGPLIntercept","zetaFFGPLExponent","zetaFFGPLAIC","GammaShapeParameter","GammaShapeSE","GammaRateParameter","GammaRateSE","LULow","LUHigh")
zetaAnalysis <- head(zetaAnalysis,-1)
write.table(zetaAnalysis,"ZetaAndFFGLUTrends.txt",quote=FALSE,sep="\t",row.names = FALSE)


#####################################################################

#Run through analysis on SCCWRP archive on a watershed-level scale.
#Select watersheds with a large enough set of samples for analysis.
watersheds <- subset(as.data.frame(table(SCCWRP$Watershed)),Freq>=groupNum)
colnames(watersheds) <- c("Watershed","Samples")
#Get samples only found in more heavily sampled watersheds.
GISBioDataLargeWS <- subset(GISBioData,Watershed %in% watersheds$Watershed)

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
#Get the samples per watershed for watersheds with at least a certain number of samples.
LargeWatersheds <- subset(watersheds,Samples>=20)
#Taxa frequency table.
taxaFreq <- as.data.frame(table(GISBioData$FinalID))
colnames(taxaFreq) <- c("FinalID","Freq")
#Find the total number of taxa in the full data set.
taxaMax <- length(unique(GISBioData$FinalID))
#Get number of unique LU_2000_5K values.
sitesNum <- length(unique(GISBioData$UniqueID))
#Enter number of divisions for subsampling.
divisionNum = 60
#Obtain subsampling number.
sampleNum <- as.integer(sitesNum/divisionNum)
uniqueSamples <- as.data.frame(unique(GISBioData$UniqueID))
colnames(uniqueSamples) <- c("UniqueID")

zetaAnalysis <- data.frame()
for(i in 1:divisionNum){
  lowNum=(i-1)*sampleNum+1
  highNum=i*sampleNum
  GISBioData <- arrange(GISBioData,LU_2000_5K)
  uniqueSampleSubset <- as.data.frame(uniqueSamples[lowNum:highNum,1])
  colnames(uniqueSampleSubset) <- c("UniqueID")
  GISBioDataSubset <- GISBioData[GISBioData$UniqueID %in% as.vector(uniqueSampleSubset$UniqueID),]
  #Determine the average LU_2000_5K per subsample of sites.
  meanLU_2000_5K = mean(na.omit(GISBioDataSubset$LU_2000_5K))
  print(paste(lowNum,highNum,meanLU_2000_5K))
  #Initialize a data frame where the rows are all of the unique measurements for a given
  #subset of the data.
  #Order the data frame by measurement name.
  selected <- arrange(GISBioDataSubset,Year,UniqueID)
  
  #Generating a presence/absence matrix for California SCCWRP data.
  eLSAInput <- as.data.frame(unique(selected$FinalID))
  colnames(eLSAInput)<-c("FinalID")
  eLSAInput <- as.data.frame(eLSAInput[order(as.character(eLSAInput$FinalID)),])
  colnames(eLSAInput)<-c("FinalID")
  taxa <- eLSAInput
  #Add the relative taxa abundances by column to a new dataframe.
  #The rows are the unique taxa in a given subset of data.
  selected <- selected[order(selected$Year,selected$UniqueID,selected$FinalID),]
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
  }
  
  #Generate a presence/absence dataframe for zeta diversity analysis.
  #Rows for samples, columns for taxa IDs.
  eLSAInput[is.na(eLSAInput)] <- 0
  eLSANames <- eLSAInput$FinalID
  data.SCCWRP <- as.data.frame(t(eLSAInput[,-c(1)]))
  colnames(data.SCCWRP) <- eLSANames
  data.SCCWRP[data.SCCWRP > 0] <- 1
  
  #Computes zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
  #using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
  zetaDecay <- Zeta.decline.mc(data.SCCWRP,xy=NULL,orders=1:10,sam=1000)
  
  dat <- data.frame()
  dat[1,1] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
  dat[1,2] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
  dat[1,3] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
  dat[1,4] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
  dat[1,5] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
  dat[1,6] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
  dat[1,7] <- meanLU_2000_5K #Mean land use.
  
  zetaAnalysis <- rbind(zetaAnalysis,dat)
  print(dat)
}
colnames(zetaAnalysis) <- c("ZetaExponentialIntercept","ZetaExponentialExponent","ZetaExponentialAIC","ZetaPLIntercept","ZetaPLExponent","ZetaPLAIC","meanLU_2000_5K")
write.table(zetaAnalysis,"LU_2000_5KSiteSweepCAZeta.txt",quote=FALSE,sep="\t",row.names = FALSE)
