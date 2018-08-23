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

#Look through watersheds which are sufficiently sampled and do the following:
#Divide the entire data set, composed of heavily sampled watersheds, into quantiles based on land use.
#Calculate zeta diversity of each land use band.
#Determine the histogram of the number of watersheds each co-occurrence occurs in.
#Calculate the parameters of the gamma distribution fit to this histogram.
zetaAnalysis <- data.frame()
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
    FFGrand <- FFGInput
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
      #Randomly assign functional feeding groups to their sample counts to eventually test how
      #far from random their relative abundances are.
      tmp3 <- tmp2[sample(nrow(tmp2)),]
      tmp3$FunctionalFeedingGroup <- tmp2$FunctionalFeedingGroup
      colnames(tmp3) <-  c("FunctionalFeedingGroup",ID)
      FFGrand <- cbind(FFGrand,tmp3)
      FFGrand <-  FFGrand[,!duplicated(colnames(FFGrand))]
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
    
    dat <- data.frame()
    #Compute zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
    #using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
    zetaDecay <- Zeta.decline.mc(data.SCCWRP,xy=NULL,orders=1:10,sam=1000)
    dat[1,1] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    dat[1,2] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    dat[1,3] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    dat[1,4] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
    dat[1,5] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    dat[1,6] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    #Compute zeta diversity, the number of functional feeding groups shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
    #using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
    zetaDecay <- Zeta.decline.mc(ffg.SCCWRP,xy=NULL,orders=1:10,sam=1000)
    dat[1,7] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    dat[1,8] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    dat[1,9] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    dat[1,10] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
    dat[1,11] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    dat[1,12] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    #Compute zeta diversity, the number of functional feeding groups shared by multiple assemblages, for a range of orders (number of assemblages or sites), 
    #using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
    zetaDecay <- Zeta.decline.mc(ffg.rand.SCCWRP,xy=NULL,orders=1:10,sam=1000)
    dat[1,13] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    dat[1,14] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    dat[1,15] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    dat[1,16] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
    dat[1,17] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    dat[1,18] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    
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
    dat[1,19] <- as.numeric(histDecay$estimate[1]) #Gamma distribution histogram fit shape parameter.
    dat[1,20] <- as.numeric(histDecay$sd[1]) #Gamma distribution histogram fit shape parameter standard error.
    dat[1,21] <- as.numeric(histDecay$estimate[2])#Gamma distribution histogram fit rate parameter.
    dat[1,22] <- as.numeric(histDecay$sd[2]) #Gamma distribution histogram fit rate parameter standard error.
    dat[1,23] <- LULow
    dat[1,24] <- LUHigh
    #Get relative abundances of taxa by functional feeding groups across a set of samples.
    FFGTotals <- t(as.data.frame(rowSums(FFGInput[,2:ncol(FFGInput)]) / sum(rowSums(FFGInput[,2:ncol(FFGInput)]))))
    rownames(FFGTotals) <- 1:nrow(FFGTotals)
    dat <- cbind(dat,FFGTotals)
    print(dat)
  }
  zetaAnalysis <- rbind(zetaAnalysis,dat)
}
colnames(zetaAnalysis) <- c("zetaExpIntercept","zetaExpExponent","zetaExpAIC","zetaPLIntercept","zetaPLExponent","zetaPLAIC","zetaFFGExpIntercept","zetaFFGExpExponent","zetaFFGExpAIC","zetaFFGPLIntercept","zetaFFGPLExponent","zetaFFGPLAIC","zetaFFGrandExpIntercept","zetaFFGrandExpExponent","zetaFFGrandExpAIC","zetaFFGrandPLIntercept","zetaFFGrandPLExponent","zetaFFGrandPLAIC","GammaShapeParameter","GammaShapeSE","GammaRateParameter","GammaRateSE","LULow","LUHigh","CFra","CGra","MHra","OMra","Pra","PHra","SCra","SHra")
zetaAnalysis <- head(zetaAnalysis,-1)
write.table(zetaAnalysis,"ZetaAndFFGLUTrends.txt",quote=FALSE,sep="\t",row.names = FALSE)
zetaAnalysis <- read.table("ZetaAndFFGLUTrends.txt",header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#####################################################################

#This portion focuses on generating co-occurrence networks on a HUC-8 watershed scale within the SCCWRP archive.
#Subsetting waterhedss by land use bands to check for uniformity of co-occurrence network formation
#within similar watersheds to changes in land use.
sampleNum <- 20 #Number of samples per watershed by land use band to use to generate a co-occurrence network.
LUquantile <- quantile(GISBioDataLargeWS$LU_2000_5K,probs=seq(0,1,0.1))#To get land use quantiles.
taxa <- as.data.frame(sort(unique(GISBioDataLargeWS$FinalID)))#Get unique taxa in full data set.
colnames(taxa) <- c("FinalID")
eLSAInput <- taxa
FFgroups <- as.data.frame(sort(unique(GISBioDataLargeWS$FunctionalFeedingGroup)))#Get unique functional feeding groups in full data set.
FFgroups <- na.omit(FFgroups)
colnames(FFgroups) <- c("FunctionalFeedingGroup")
FFGInput <- FFgroups
ShellCommand <- as.data.frame(matrix(nrow=0,ncol=1)) #Bind eLSA commands into a data frame to write as a series of shell commands for running each analysis on a cluster.
for(WS in unique(GISBioDataLargeWS$Watershed)){
  for(i in 1:length(LUquantile)){
    LULow <- as.numeric(LUquantile[i])
    if(i<length(LUquantile)){
      LUHigh <- as.numeric(LUquantile[i+1])
    }
    if(i==length(LUquantile)){
      LUHigh == 100
    }
    if(LULow == 0 & LUHigh == 0){
      LUSubset <- subset(GISBioDataLargeWS,LU_2000_5K==LULow & Watershed==WS) #Subset samples by aggregated land use and watershed.
      #print(paste(i-1,i,LULow,LUHigh,WS,length(unique(LUSubset$UniqueID))))
    }
    if(LULow != LUHigh){
      LUSubset <- subset(GISBioDataLargeWS,LU_2000_5K>=LULow & LU_2000_5K < LUHigh & Watershed==WS) #Subset samples by aggregated land use and watershed.
      #print(paste(i-1,i,LULow,LUHigh,WS,length(unique(LUSubset$UniqueID))))
    }
    if(i < length(LUquantile) & length(unique(LUSubset$UniqueID)) >= sampleNum){
      sampleNames <- sample(unique(LUSubset$UniqueID),sampleNum)
      selected <- subset(LUSubset, UniqueID %in% sampleNames)
      selected <- arrange(selected,Year,UniqueID)
      #Get taxonomic diversity for the same set of samples within a given land use band.
      eLSAInput <- as.data.frame(sort(unique(selected$FinalID)))
      colnames(eLSAInput) <- c("FinalID")
      #Get functional feeding group counts for the same set of samples within a given land use band.
      FFGInput <- as.data.frame(sort(unique(selected$FunctionalFeedingGroup)))
      colnames(FFGInput) <- c("FunctionalFeedingGroup")
      FFGInput <- na.omit(FFGInput)
      #Generate input data frames for co-occurrence network generation with eLSA.
      for(ID in unique(selected$UniqueID)){
        #Add the relative taxa abundances by column to a new dataframe.
        #The rows are the unique taxa in a given subset of data.
        tmp <- filter(selected, UniqueID == ID)[,c("FinalID","Measurement","UniqueID")]
        tmp <- as.data.frame(tmp[order(tmp$FinalID),])
        tmp <- tmp[,c("FinalID","Measurement")]
        colnames(tmp) <- c("FinalID",ID)
        tmp <- tmp %>% group_by(FinalID) %>% summarise_if(is.numeric,mean,na.rm=TRUE)
        tmp <- join(tmp,taxa,type="full",by=c("FinalID"))
        tmp <- as.data.frame(tmp[order(tmp$FinalID),])
        eLSAInput <- join(eLSAInput,tmp,by=c("FinalID"))
        eLSAInput <- eLSAInput[,!duplicated(colnames(eLSAInput))]
        #Compute functional feeding group diversity by sample and sample grouping.
        tmp2 <- filter(selected, UniqueID == ID)[,c("FunctionalFeedingGroup","Count","UniqueID")]
        tmp2 <- as.data.frame(tmp2[order(tmp2$FunctionalFeedingGroup),])
        tmp2 <- tmp2[!is.na(tmp2$FunctionalFeedingGroup),]
        tmp2 <- tmp2[,c("FunctionalFeedingGroup","Count")]
        colnames(tmp2) <-  c("FunctionalFeedingGroup",ID)
        tmp2 <- tmp2 %>% group_by(FunctionalFeedingGroup) %>% summarise_if(is.numeric,sum,na.rm=TRUE)
        tmp2[,2] <- tmp2[,2]/sum(na.omit(tmp2[,2]))
        tmp2 <- join(tmp2,FFgroups,type="full",by=c("FunctionalFeedingGroup"))
        tmp2 <- as.data.frame(tmp2[order(tmp2$FunctionalFeedingGroup),])
        FFGInput <- join(FFGInput,tmp2,by=c("FunctionalFeedingGroup"))
        FFGInput <-  FFGInput[,!duplicated(colnames(FFGInput))]
      }
      
      #To generate genera relative abundances data for eLSA.
      #How many years are in each set of samples, and how many samples were taken by year?
      SamplesByYear <- as.data.frame(setDT(selected)[, .(count = uniqueN(UniqueID)), by = Year])
      #Determine the number of time points in the eLSA input file for genera relative abundance data.
      spotNum = nrow(SamplesByYear)
      #Determine the number of replicates per time point in the eLSA input file.
      #In order to ensure a uniform number of replicates per year this needs to
      #be the maximum number of replicates for all of the years available.
      repNum = max(SamplesByYear$count)
      #Now insert the replicates with actual data in between the "NA" dummy columns
      #which ensure that the final eLSA input file has an even number of replicates
      #per year regardless of the variations in the actual number of sites (replicates)
      #sampled per year.
      eLSAtmp <- as.data.frame(eLSAInput[,1])
      colnames(eLSAtmp) <- c("FinalID")
      j=1
      k=1
      nulCol <- data.frame(matrix(ncol=repNum*spotNum+1,nrow=length(na.omit(unique(selected$FinalID)))))
      nulCol[,1] <- eLSAInput[,1]
      eLSANames <- c("FinalID")
      for(year in unique(selected$Year)){
        tmp <- filter(selected, Year == year)
        rep = length(unique(tmp$UniqueID))
        for(i in 1:repNum){
          if(i <= rep){
            repLabel = paste(year,"DoneRep",i,sep="")
            eLSANames <- c(eLSANames,repLabel)
            j=j+1
            k=k+1
            eLSAtmp[,k] <- eLSAInput[,j]
          }
          else{
            repLabel = as.character(paste(year,"Rep",i,sep=""))
            eLSANames <- c(eLSANames,repLabel)
            k=k+1
            eLSAtmp[,k] <- NA
          }
        }
      }
      eLSAInput <- eLSAtmp
      colnames(eLSAInput) <- eLSANames
      
      #To generate functional feeding group relative abundances data for eLSA.
      #How many years are in each set of samples, and how many samples were taken by year?
      SamplesByYear <- as.data.frame(setDT(selected)[, .(count = uniqueN(UniqueID)), by = Year])
      #Determine the number of time points in the eLSA input file for genera relative abundance data.
      spotNum = nrow(SamplesByYear)
      #Determine the number of replicates per time point in the eLSA input file.
      #In order to ensure a uniform number of replicates per year this needs to
      #be the maximum number of replicates for all of the years available.
      repNum = max(SamplesByYear$count)
      #Now insert the replicates with actual data in between the "NA" dummy columns
      #which ensure that the final eLSA input file has an even number of replicates
      #per year regardless of the variations in the actual number of sites (replicates)
      #sampled per year.
      FFGtmp <- as.data.frame(FFGInput[,1])
      colnames(FFGtmp) <- c("FunctionalFeedingGroup")
      j=1
      k=1
      nulCol <- data.frame(matrix(ncol=repNum*spotNum+1,nrow=length(na.omit(unique(selected$FunctionalFeedingGroup)))))
      nulCol[,1] <- FFGInput[,1]
      FFGNames <- c("FunctionalFeedingGroup")
      for(year in unique(selected$Year)){
        tmp <- filter(selected, Year == year)
        rep = length(unique(tmp$UniqueID))
        for(i in 1:repNum){
          if(i <= rep){
            repLabel = paste(year,"DoneRep",i,sep="")
            FFGNames <- c(FFGNames,repLabel)
            j=j+1
            k=k+1
            FFGtmp[,k] <- FFGInput[,j]
          }
          else{
            repLabel = as.character(paste(year,"Rep",i,sep=""))
            FFGNames <- c(FFGNames,repLabel)
            k=k+1
            FFGtmp[,k] <- NA
          }
        }
      }
      FFGInput <- FFGtmp
      colnames(FFGInput) <- FFGNames
      
      #Output files for co-occurrence network generation with eLSA for both relative abundances
      #of genera by sample and for functional feeding groups by sample.
      eLSAFilename <- paste("GeneraAbundancesWatershed",WS,"LU",LULow,"to",LUHigh,"SampleNum",sampleNum,sep="")
      FFGFilename <- paste("FFGAbundancesWatershed",WS,"LU",LULow,"to",LUHigh,"SampleNum",sampleNum,sep="")
      write.table(eLSAInput,paste(eLSAFilename,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
      write.table(FFGInput,paste(FFGFilename,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
      eLSACommand = paste("lsa_compute ",eLSAFilename,".txt -r ",repNum," -s ",spotNum," ",eLSAFilename,"Network.txt;",sep="")
      print(eLSACommand)
      FFGCommand = paste("lsa_compute ",FFGFilename,".txt -r ",repNum," -s ",spotNum," ",FFGFilename,"Network.txt;",sep="")
      print(FFGCommand)
      ShellCommand <- rbindlist(list(ShellCommand,data.table(eLSACommand),data.table(FFGCommand)), use.names=FALSE)
      #print(paste(i-1,i,LULow,LUHigh,WS,length(unique(LUSubset$UniqueID)),nrow(eLSAInput),nrow(FFGInput)))
    }
  }
}

#Split command outputs into parts to make smaller lists of shell commands.
div <- 4 #Number of shell scripts to run eLSA commands.
n <- nrow(ShellCommand)/div #Number of eLSA commands per script.
nr <- nrow(ShellCommand)
test <- split(ShellCommand, rep(1:ceiling(nr/n), each=n, length.out=nr))
for(i in 1:div){
  write.table(test[i],paste(i,".sh",sep=""),quote=FALSE,sep="\t",row.names = FALSE,col.names=FALSE)
}

##############################################################################
#Read in eLSA output.
#Compute network statistics of the likeliest association networks between taxa.
library(igraph)
library(network)
library(stringr)
choice=1# 1 for functional feeding group networks and 2 for genera networks.
#Read in co-occurrence networks and analyze their topologies.
if(choice==1){networkfiles <- Sys.glob("FFGAbundances*Network.txt")}#Functional feeding group networks.
FFGTallyCov <- data.frame() #Used to tally organism counts by functional feeding groups per set of samples used to generate co-occurrence networks.
FFGTallyCon <- data.frame()
if(choice==2){networkfiles <- Sys.glob("GeneraAbundances*Network.txt")}#Genera networks.
networkAnalysis <- data.frame()
networkConTaxa <- data.frame()
networkCovTaxa <- data.frame()
#Define a 'not in' function.
'%!in%' <- function(x,y)!('%in%'(x,y))
for(networkFile in networkfiles){
  networkdata <- read.delim(networkFile,header=TRUE, sep="\t",as.is=T,check.names=FALSE)
  #Filter out association network data based on P and Q scores, for the local similarity
  #between two factors, with values less than a particuar threshold.
  networkdata <- filter(networkdata, P <= 5e-2)
  networkdata <- filter(networkdata, Q <= 5e-2)
  names(networkdata)[names(networkdata)=="LS"]<-"weight"
  MidLU <- 0.5*(as.numeric(str_match(networkFile,"(.*?)AbundancesWatershed(.*?)LU(.*?)to(.*?)SampleNum20Network.txt")[4])+as.numeric(str_match(networkFile,"(.*?)AbundancesWatershed(.*?)LU(.*?)to(.*?)SampleNum20Network.txt")[5]))
  if(nrow(networkdata) > 0){
    #Generate network graph and begin calculating network parameters.
    networkgraph=graph.data.frame(networkdata,directed=FALSE)
    Network_size<-network.size(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
    #Get the full weighted adjacency matrix.
    networkmatrix <- as.matrix(get.adjacency(networkgraph,attr='weight'))
    #Mean interaction strength
    meanStrength <- mean(abs(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_network <- eigen(networkmatrix)
    #Get the real component first eigenvalue.
    lambda_network_m <- Re(lambda_network$values[1])
    #Generate randomized version of full weighted adjacency matrix.
    set.seed(1)
    randnetworkmatrix <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_rand <- eigen(randnetworkmatrix)
    #Get the real component of the first eigenvalue.
    lambda_rand_m <- Re(lambda_rand$values[1])
    #Calculate stability parameter.
    gamma <- lambda_network_m/lambda_rand_m
    #Calculate the degree heterogeneity.
    networkmatrix[upper.tri(networkmatrix)] <- 0
    networkmatrix <- ifelse(networkmatrix!=0,1,networkmatrix)
    zeta <- mean(colSums(networkmatrix)^2)/mean(colSums(networkmatrix))^2
    #Calculate modularity
    networkModularity <- modularity(cluster_edge_betweenness(networkgraph, weights=NULL,directed=FALSE))
    M <- networkModularity
    networkNodecount <-network.size(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
    # Get the number of unique network edges
    networkEdgecount <- network.edgecount(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
    # Get the number of nodes
    networkNodecount <- network.size(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
    # Get the average degree per node.
    k <- (2*networkEdgecount)/networkNodecount
    # Calculate the modularity of the random network.
    networkRandModularity <- (1-(2/sqrt(networkNodecount)))*(2/k)^(2/3)
    # Calculate the log ratio of the modularities.
    l_rM <- log(networkModularity/networkRandModularity)
    # Calculate connectance
    adj <- as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
    C <- network.density(adj)
    # Calculate the average network path length
    L <- mean_distance(networkgraphCon,directed=FALSE)
    # Calculate the log-ratio of the average network path length to that of its corresponding random network.
    networkRandLength <- 0.5+((log(networkNodecount)-0.5772156649)/log(k))
    l_rL <- log(L/networkRandLength)
    # Calculate the clustering coefficient
    Cl <- transitivity(networkgraphCon,type="globalundirected",isolate="zero")
    # Calculate the log-ratio of the clustering coefficient to that of its corresponding random network.
    l_rCl <- log(Cl/(k/networkNodecount))
  }
  FFGTallyCon <- data.frame()
  #Filter contravariant network data based on local similarity scores.
  networkdataCon <- subset(networkdata,networkdata$weight<0)
  if(nrow(networkdataCon)>0){
    #Aggregate significantly contravarying taxa.
    networkdataConTemp <- networkdataCon[,c("X","Y","weight")]
    #Count number of edges per functional feeding group.
    if(choice==1){
      FFGTallyCon <- as.data.table(table(unlist(networkdataCon[,c("X","Y")])))
      colnames(FFGTallyCon) <- c("FunctionalFeedingGroup","N_con")
      FFGTallyCon <- merge(FFGTallyCon,FFgroups,all=TRUE)
      FFGTallyCon$FunctionalFeedingGroup <- as.character(FFGTally$FunctionalFeedingGroup)
      FFGTallyCon$N_con <- as.numeric(as.character(FFGTallyCon$N_con))
      FFGTallyCon <- arrange(FFGTallyCon,FunctionalFeedingGroup) 
    }
    networkdataConTemp <- as.data.frame(table(append(networkdataConTemp$X,networkdataConTemp$Y,after=length(networkdataConTemp$X))))
    networkdataConTemp$MidLU <- MidLU
    networkConTaxa <- rbind(networkConTaxa,networkdataConTemp)
    #Generate network graph and begin calculating network parameters.
    networkgraphCon=graph.data.frame(networkdataCon,directed=FALSE)
    #Get the full weighted adjacency matrix.
    networkmatrix <- as.matrix(get.adjacency(networkgraphCon,attr='weight'))
    #Mean interaction strength
    meanStrength_Con <- mean(abs(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_network <- eigen(networkmatrix)
    #Get the real component first eigenvalue.
    lambda_network_m_Con <- Re(lambda_network$values[1])
    #Generate randomized version of full weighted adjacency matrix.
    set.seed(1)
    randnetworkmatrix <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_rand_Con <- eigen(randnetworkmatrix)
    #Get the real component of the first eigenvalue.
    lambda_rand_Con <- Re(lambda_rand_Con$values[1])
    #Calculate stability parameter.
    gamma_Con <- lambda_network_m_Con/lambda_rand_Con
    #Calculate the degree heterogeneity.
    networkmatrixCon <- networkmatrix
    networkmatrixCon[upper.tri(networkmatrixCon)] <- 0
    networkmatrixCon <- ifelse(networkmatrixCon!=0,1,networkmatrixCon)
    zeta_Con <- mean(colSums(networkmatrixCon)^2)/mean(colSums(networkmatrixCon))^2
    #Calculate the degree heterogeneity of the corresponding random network.
    randnetworkmatrixCon <- randnetworkmatrix
    randnetworkmatrixCon[upper.tri(randnetworkmatrixCon)] <- 0
    randnetworkmatrixCon <- ifelse(randnetworkmatrixCon!=0,1,randnetworkmatrixCon)
    zeta_rand_Con <- mean(colSums(randnetworkmatrixCon)^2)/mean(colSums(randnetworkmatrixCon))^2
    # Log response ratio of degree heterogeneity.
    l_con_rzeta <- log(zeta_Con/zeta_rand_Con)
    # Generate adjacency matrix of relative taxa abundance correlations
    adj= as.network(get.adjacency(networkgraphCon,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
    # Get the number of unique network edges
    networkEdgecount <- network.edgecount(adj)
    networkEdgecountCon <- networkEdgecount
    # Get the number of nodes
    networkNodecount <- network.size(adj)
    Network_sizeCon <- networkNodecount
    # Get the average degree per node.
    k <- (2*networkEdgecount)/networkNodecount
    # Get the random characteristic path length.
    networkRandLength <- 0.5+((log(networkNodecount)-0.5772156649)/log(k))
    # Get the random clustering coefficient.
    networkRandClustering <- k/networkNodecount
    # Get the network density.
    networkDensity <- network.density(adj)
    con_C <- networkDensity
    # Calculate the modularity of the network.
    networkModularity <- modularity(cluster_edge_betweenness(networkgraphCon, weights=NULL,directed=FALSE))
    con_M <- networkModularity
    # Calculate the number of groups related to the modularity value.
    networkModGroups <- length(cluster_edge_betweenness(networkgraphCon, weights=NULL,directed=FALSE))
    # Calculate the average network path length
    networkLength <- mean_distance(networkgraphCon,directed=FALSE)
    con_L <- networkLength
    # Calculate the clustering coefficient
    networkClustering <- transitivity(networkgraphCon,type="globalundirected",isolate="zero")
    con_Cl <- networkClustering
    # Calcuate the log ratio of clustering coefficients.
    l_con_rCl <- log(networkClustering/networkRandClustering)
    # Calculate the modularity of the random network.
    networkRandModularity <- (1-(2/sqrt(networkNodecount)))*(2/k)^(2/3)
    # Calculate the log ratio of the modularities.
    l_con_rM <- log(networkModularity/networkRandModularity)
    # Get log ratio of characteristic path lengths.
    l_con_rL <- log(networkLength/networkRandLength)
  }
  #Filter covariant network data based on local similarity scores.
  networkdataCov <- subset(networkdata,networkdata$weight>0)
  if(nrow(networkdataCov)>0){
    #Aggregate significantly contravarying taxa.
    networkdataCovTemp <- networkdataCov[,c("X","Y","weight")]
    #Count number of edges per functional feeding group.
    if(choice==1){
      FFGTallyCov <- as.data.table(table(unlist(networkdataCov[,c("X","Y")])))
      colnames(FFGTallyCov) <- c("FunctionalFeedingGroup","N_cov")
      FFGTallyCov <- merge(FFGTallyCov,FFgroups,all=TRUE)
      FFGTallyCov$FunctionalFeedingGroup <- as.character(FFGTally$FunctionalFeedingGroup)
      FFGTallyCov$N_cov <- as.numeric(as.character(FFGTallyCov$N_cov))
      FFGTallyCov <- arrange(FFGTallyCov,FunctionalFeedingGroup) 
    }
    networkdataCovTemp <- as.data.frame(table(append(networkdataCovTemp$X,networkdataCovTemp$Y,after=length(networkdataCovTemp$X))))
    networkdataCovTemp$MidLU <- MidLU
    networkCovTaxa <- rbind(networkCovTaxa,networkdataCovTemp)
    #Generate network graph and begin calculating network parameters.
    networkgraphCov=graph.data.frame(networkdataCov,directed=FALSE)
    #Get the full weighted adjacency matrix.
    networkmatrix <- as.matrix(get.adjacency(networkgraphCov,attr='weight'))
    #Mean interaction strength
    meanStrength_Cov <- mean(abs(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_network <- eigen(networkmatrix)
    #Get the real component first eigenvalue.
    lambda_network_m_Cov <- Re(lambda_network$values[1])
    #Generate randomized version of full weighted adjacency matrix.
    set.seed(1)
    randnetworkmatrix <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_rand_Cov <- eigen(randnetworkmatrix)
    #Get the real component of the first eigenvalue.
    lambda_rand_Cov <- Re(lambda_rand_Cov$values[1])
    #Calculate stability parameter.
    gamma_Cov <- lambda_network_m_Cov/lambda_rand_Cov
    #Calculate the degree heterogeneity.
    networkmatrixCov <- networkmatrix
    networkmatrixCov[upper.tri(networkmatrixCov)] <- 0
    networkmatrixCov <- ifelse(networkmatrixCov!=0,1,networkmatrixCov)
    zeta_Cov <- mean(colSums(networkmatrixCov)^2)/mean(colSums(networkmatrixCov))^2
    #Calculate the degree heterogeneity of the corresponding random network.
    randnetworkmatrixCov <- randnetworkmatrix
    randnetworkmatrixCov[upper.tri(randnetworkmatrixCov)] <- 0
    randnetworkmatrixCov <- ifelse(randnetworkmatrixCov!=0,1,randnetworkmatrixCov)
    zeta_rand_Cov <- mean(colSums(randnetworkmatrixCov)^2)/mean(colSums(randnetworkmatrixCov))^2
    # Log response ratio of degree heterogeneity.
    l_cov_rzeta <- log(zeta_Cov/zeta_rand_Cov)
    # Generate adjacency matrix of relative taxa abundance correlations
    adj= as.network(get.adjacency(networkgraphCov,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
    # Get the number of unique network edges
    networkEdgecount <- network.edgecount(adj)
    networkEdgecountCov <- networkEdgecount
    # Get the number of nodes
    networkNodecount <- network.size(adj)
    Network_sizeCov <-  networkNodecount
    # Get the average degree per node.
    k <- (2*networkEdgecount)/networkNodecount
    # Get the random characteristic path length.
    networkRandLength <- 0.5+((log(networkNodecount)-0.5772156649)/log(k))
    # Get the random clustering coefficient.
    networkRandClustering <- k/networkNodecount
    # Get the network density.
    networkDensity <- network.density(adj)
    cov_C <- networkDensity
    # Calculate the modularity of the network.
    networkModularity <- modularity(cluster_edge_betweenness(networkgraphCov, weights=NULL,directed=FALSE))
    cov_M <- networkModularity
    # Calculate the number of groups related to the modularity value.
    networkModGroups <- length(cluster_edge_betweenness(networkgraphCov, weights=NULL,directed=FALSE))
    # Calculate the average network path length
    networkLength <- mean_distance(networkgraphCov,directed=FALSE)
    cov_L <- networkLength
    # Calculate the clustering coefficient
    networkClustering <- transitivity(networkgraphCov,type="globalundirected",isolate="zero")
    cov_Cl <- networkClustering
    # Calcuate the log ratio of clustering coefficients.
    l_cov_rCl <- log(networkClustering/networkRandClustering)
    # Calculate the modularity of the random network.
    networkRandModularity <- (1-(2/sqrt(networkNodecount)))*(2/k)^(2/3)
    # Calculate the log ratio of the modularities.
    l_cov_rM <- log(networkModularity/networkRandModularity)
    # Get log ratio of characteristic path lengths.
    l_cov_rL <- log(networkLength/networkRandLength)
  }
  Watershed <- str_match(networkFile,"(.*?)AbundancesWatershed(.*?)LU(.*?)to(.*?)SampleNum20Network.txt")[3]
  if(nrow(networkdata) > 0 & nrow(networkdataCov) > 0 & nrow(networkdataCon) > 0){
    Pm <- networkEdgecountCov/(networkEdgecountCov+networkEdgecountCon)
    if(choice==1){
      CF_cov <- FFGTallyCov[FFGTallyCov$FunctionalFeedingGroup=="CF",2]
      CG_cov <- FFGTallyCov[FFGTallyCov$FunctionalFeedingGroup=="CG",2]
      MH_cov <- FFGTallyCov[FFGTallyCov$FunctionalFeedingGroup=="MH",2]
      OM_cov <- FFGTallyCov[FFGTallyCov$FunctionalFeedingGroup=="OM",2]
      P_cov <- FFGTallyCov[FFGTallyCov$FunctionalFeedingGroup=="P",2]
      PH_cov <- FFGTallyCov[FFGTallyCov$FunctionalFeedingGroup=="PH",2]
      SC_cov <- FFGTallyCov[FFGTallyCov$FunctionalFeedingGroup=="SC",2]
      SH_cov <- FFGTallyCov[FFGTallyCov$FunctionalFeedingGroup=="SH",2]
      CF_con <- FFGTallyCon[FFGTallyCon$FunctionalFeedingGroup=="CF",2]
      CG_con <- FFGTallyCon[FFGTallyCon$FunctionalFeedingGroup=="CG",2]
      MH_con <- FFGTallyCon[FFGTallyCon$FunctionalFeedingGroup=="MH",2]
      OM_con <- FFGTallyCon[FFGTallyCon$FunctionalFeedingGroup=="OM",2]
      P_con <- FFGTallyCon[FFGTallyCon$FunctionalFeedingGroup=="P",2]
      PH_con <- FFGTallyCon[FFGTallyCon$FunctionalFeedingGroup=="PH",2]
      SC_con <- FFGTallyCon[FFGTallyCon$FunctionalFeedingGroup=="SC",2]
      SH_con <- FFGTallyCon[FFGTallyCon$FunctionalFeedingGroup=="SH",2]
      dat <- list(networkFile,Watershed,MidLU,Network_size,Network_sizeCov,Network_sizeCon,C,cov_C,con_C,Cl,cov_Cl,con_Cl,l_rCl,l_cov_rCl,l_con_rCl,L,cov_L,con_L,l_rL,l_cov_rL,l_con_rL,M,cov_M,con_M,zeta,zeta_Cov,zeta_Con,meanStrength,meanStrength_Cov,meanStrength_Con,Pm,CF_cov,CG_cov,MH_cov,OM_cov,P_cov,PH_cov,SC_cov,SH_cov,CF_con,CG_con,MH_con,OM_con,P_con,PH_con,SC_con,SH_con)
    }
    if(choice==2){
      dat <- list(networkFile,Watershed,MidLU,Network_size,Network_sizeCov,Network_sizeCon,C,cov_C,con_C,Cl,cov_Cl,con_Cl,l_rCl,l_cov_rCl,l_con_rCl,L,cov_L,con_L,l_rL,l_cov_rL,l_con_rL,M,cov_M,con_M,zeta,zeta_Cov,zeta_Con,meanStrength,meanStrength_Cov,meanStrength_Con,Pm)
    }
    dat <- data.frame(t(unlist(dat)))
    networkAnalysis <- rbind(networkAnalysis,dat)
    print(dat)
  }
}
if(choice==1){
  colnames(networkAnalysis) <- c("networkFile","Watershed","MidLU","Network_size","Network_sizeCov","Network_sizeCon","C","cov_C","con_C","Cl","cov_Cl","con_Cl","l_rCl","l_cov_rCl","l_con_rCl","L","cov_L","con_L","l_rL","l_cov_rL","l_con_rL","M","cov_M","con_M","zeta","zeta_Cov","zeta_Con","meanStrength","meanStrength_Cov","meanStrength_Con","Pm","CF_cov","CG_cov","MH_cov","OM_cov","P_cov","PH_cov","SC_cov","SH_cov","CF_con","CG_con","MH_con","OM_con","P_con","PH_con","SC_con","SH_con")
}
if(choice==2){
  colnames(networkAnalysis) <- c("networkFile","Watershed","MidLU","Network_size","Network_sizeCov","Network_sizeCon","C","cov_C","con_C","Cl","cov_Cl","con_Cl","l_rCl","l_cov_rCl","l_con_rCl","L","cov_L","con_L","l_rL","l_cov_rL","l_con_rL","M","cov_M","con_M","zeta","zeta_Cov","zeta_Con","meanStrength","meanStrength_Cov","meanStrength_Con","Pm")
}

networkAnalysis[networkAnalysis=="-Inf"] <- NA
networkAnalysis[networkAnalysis=="Inf"] <- NA
networkAnalysis[networkAnalysis=="NaN"] <- NA
cols <- 3:ncol(networkAnalysis)
networkAnalysis[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
networkAnalysis <- arrange(networkAnalysis,Watershed,MidLU)

dev.off()
for(WS in unique(networkAnalysis$Watershed)){
  test <- subset(networkAnalysis,Watershed==WS)
  plot(networkAnalysis$MidLU,networkAnalysis$CF_cov,type="l",xlim=c(0,100),ylim=c(0,3))
  par(new=T)
}

#Regression between network parameters.
library(Hmisc)
library(corrplot)
library("PerformanceAnalytics")
#Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
chart.Correlation(networkAnalysis[,c("MidLU","zeta_Cov","cov_M","cov_C","meanStrength_Cov","l_cov_rL","l_cov_rCl","Network_sizeCov")], histogram=FALSE, method=c("spearman"))
chart.Correlation(networkAnalysis[,c("MidLU","zeta_Con","con_M","con_C","meanStrength_Con","l_con_rL","l_con_rCl","Network_sizeCon")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("MidLU","zeta","M","C","meanStrength","l_rL","l_rCl","Network_size","Pm")], histogram=FALSE, method=c("spearman"))

chart.Correlation(networkAnalysis[,c("MidLU","CF_cov","CG_cov","MH_cov","OM_cov","P_cov","PH_cov","SC_cov","SH_cov")], histogram=FALSE, method=c("spearman"))
chart.Correlation(networkAnalysis[,c("MidLU","CF_con","CG_con","MH_con","OM_con","P_con","PH_con","SC_con","SH_con")], histogram=FALSE, method=c("spearman"))
res.aov2 <- aov(zeta ~ Watershed + MidLU + Watershed:MidLU, data=networkAnalysis)
summary(res.aov2)
plot(res.aov2,1)
plot(res.aov2,2)
