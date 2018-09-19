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
library(relaimpo)

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
#HUC2 = State level, HUC4 = super-regional level, HUC6 = regional level, HUC8 = local level
Watersheds <- read.table("SCCWRPWatersheds.tsv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
Watersheds$Watershed <- gsub('\\s+', '', Watersheds$Watershed)
#Merge in watershed area
GISBioData <- join(GISBioData,Watersheds,by=c("Watershed"))
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
FFG <- FFG[!duplicated(FFG$FinalID),]
#Merge in functional feeding groups into sample data.
GISBioData <- join(GISBioData,FFG[,c("FinalID","FunctionalFeedingGroup")],by=c("FinalID"))
FFGCounts <- na.omit(as.data.frame(unique(GISBioData$FunctionalFeedingGroup)))
colnames(FFGCounts) <- c("FunctionalFeedingGroup")
FFGCounts$FunctionalFeedingGroup <- as.character(as.factor(FFGCounts$FunctionalFeedingGroup))
FFGCounts <- arrange(FFGCounts,FunctionalFeedingGroup)
#Add column containing the sum of taxa, by functional feeding groups, within each sample.
tmp <- GISBioData[,c("UniqueID","FunctionalFeedingGroup","Count")]
colnames(tmp) <- c("UniqueID","FunctionalFeedingGroup","FFGCount")
tmp <- aggregate(tmp$FFGCount, by=list(Category=tmp$UniqueID,tmp$FunctionalFeedingGroup), FUN=sum)
colnames(tmp) <- c("UniqueID","FunctionalFeedingGroup","FFGCount")
tmp <- arrange(tmp,UniqueID,FunctionalFeedingGroup)
GISBioData <- join(GISBioData,tmp,by=c("UniqueID","FunctionalFeedingGroup"))

#Find watersheds with a larger enough set of samples for downstream analysis.
sampleMin <- 25
GISBioDataLWS <- subset(GISBioData,NSamples>=sampleMin)
zetaAnalysis <- data.frame()
#Determine zeta diversity per subsample group per watershed.  Run repeated subsampling
#to help with determinining significance of correlations between environmental factors
#and zeta diversity decay parameters.
for(j in 1:100){
  for(WS in unique(GISBioDataLWS$HUC8)){
    GISBioDataLocal <- subset(GISBioDataLWS,HUC8==WS) #Subsample by watershed.
    GISBioDataLocal <- GISBioDataLocal[GISBioDataLocal$UniqueID %in% sample(unique(GISBioDataLocal$UniqueID),sampleMin),] #Subsample to a uniform sample number.
    GISBioDataLocal <- GISBioDataLocal[!is.na(GISBioDataLocal$FunctionalFeedingGroup),]
    selected <- GISBioDataLocal
    metadata <- GISBioDataLocal[,c("UniqueID","LU_2000_5K","altitude","CSCI")]
    metadata <- metadata[!duplicated(metadata),] #Get unique environmental parameters per watershed set of samples.
    #Get all unique taxa in statewide data set.
    uniqueTaxa <- as.data.frame(unique(selected$FinalID))
    colnames(uniqueTaxa) <- c("FinalID")
    uniqueTaxa <- arrange(uniqueTaxa,FinalID)
    #Get all unique taxa in statewide data set.
    uniqueFFG <- as.data.frame(unique(selected$FunctionalFeedingGroup))
    colnames(uniqueFFG) <- c("FunctionalFeedingGroup")
    uniqueFFG <- arrange(uniqueFFG,FunctionalFeedingGroup)
    #Create presence/absence matrix of taxa in samples.
    #Rows for sample ID and columns 
    PresenceAbsence <- uniqueTaxa
    #Create presence/absence matrix of FFGs in samples.
    #Rows for sample ID and columns
    PAFFG <- uniqueFFG
    for(ID in unique(selected$UniqueID)){
      #Presence/Absence matrix for taxa.
      sampleDF <- subset(selected,UniqueID == ID)
      sampleDF <- sampleDF[,c("FinalID","Count")]
      tmp <- merge(sampleDF,uniqueTaxa,by=c("FinalID"),all=TRUE)
      colnames(tmp) <- c("FinalID",ID)
      PresenceAbsence <- cbind(PresenceAbsence,tmp[,c(2)])
      colnames(PresenceAbsence)[ncol(PresenceAbsence)] <- ID
      #Presence/Absence matrix for FFGs.
      sampleDF <- subset(selected,UniqueID == ID)
      sampleDF <- sampleDF[,c("FunctionalFeedingGroup","FFGCount")]
      sampleDF <- sampleDF[!duplicated(sampleDF),]
      tmp <- merge(sampleDF,uniqueFFG,by=c("FunctionalFeedingGroup"),all=TRUE)
      colnames(tmp) <- c("FunctionalFeedingGroup",ID)
      PAFFG <- cbind(PAFFG,tmp[,c(2)])
      colnames(PAFFG)[ncol(PAFFG)] <- ID
    }
    #Generate a presence/absence dataframe for zeta diversity analysis of taxa.
    #Rows for samples, columns for taxa IDs.
    PresenceAbsence[is.na(PresenceAbsence)] <- 0
    PresenceAbsence[PresenceAbsence > 0] <- 1
    data.SCCWRP <- as.data.frame(t(PresenceAbsence[,-c(1)]))
    colnames(data.SCCWRP) <- uniqueTaxa$FinalID
    #Generate a presence/absence dataframe for zeta diversity analysis of FFGs.
    #Rows for samples, columns for FFGs IDs.
    PAFFG[is.na(PAFFG)] <- 0
    PAFFG[PAFFG > 0] <- 1
    ffg.SCCWRP <- as.data.frame(t(PAFFG[,-c(1)]))
    colnames(ffg.SCCWRP) <- uniqueFFG$FunctionalFeedingGroup
    #Generate a presence/absence dataframe for zeta diversity analysis of randomly assigned FFGs.
    #Rows for samples, columns for FFGs IDs.
    ffg.rand.SCCWRP <- as.data.frame(matrix(sample(unlist(ffg.SCCWRP),ncol(ffg.SCCWRP)*nrow(ffg.SCCWRP)),nrow=nrow(ffg.SCCWRP), byrow=T))
    colnames(ffg.rand.SCCWRP) <- colnames(ffg.SCCWRP)
    rownames(ffg.rand.SCCWRP) <- rownames(ffg.SCCWRP)
    #Calculate zeta diversity decay for taxa within each watershed set of samples.
    zetaDecay <- Zeta.decline.ex(data.SCCWRP,orders=1:10,plot=FALSE,rescale=FALSE)
    ExpIntercept <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    ExpExp <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    ExpAIC <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    PLIntercept <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
    PLExp <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    PLAIC <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    zeta_1 <- Zeta.order.ex(data.SCCWRP,order=1,rescale=FALSE) #Zeta order 1 for taxa.
    MeanAlpha <- zeta_1$zeta.val #Mean alpha diversity for taxa.
    SDAlpha <- zeta_1$zeta.val.sd #Standard deviation on alpha diversity for taxa.
    zeta_2 <- Zeta.order.ex(data.SCCWRP,order=2,rescale=FALSE) #Zeta order 2 for taxa.
    Beta <- zeta_2$zeta.val/(2*zeta_1$zeta.val-zeta_2$zeta.val) #Beta diversity for taxa based on Jaccard dissimilarity.
    #Calculate zeta diversity decay for FFGs within each watershed set of samples.
    zetaDecay <- Zeta.decline.ex(ffg.SCCWRP,orders=1:10,plot=FALSE)
    ExpInterceptFFG <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    ExpExpFFG <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    ExpAICFFG <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    PLInterceptFFG <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
    PLExpFFG <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    PLAICFFG <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    zeta_1 <- Zeta.order.ex(ffg.SCCWRP,order=1,rescale=FALSE) #Zeta order 1 for FFGs.
    MeanAlphaFFG <- zeta_1$zeta.val #Mean alpha diversity for FFGs.
    SDAlphaFFG <- zeta_1$zeta.val.sd #Standard deviation on alpha diversity for FFGs.
    zeta_2 <- Zeta.order.ex(ffg.SCCWRP,order=2,rescale=FALSE) #Zeta order 2 for FFGs.
    BetaFFG <- zeta_2$zeta.val/(2*zeta_1$zeta.val-zeta_2$zeta.val) #Beta diversity for FFGs based on Jaccard dissimilarity.
    #Calculate zeta diversity decay for randomly reassigned FFGs within each watershed set of samples.
    zetaDecay <- Zeta.decline.ex(ffg.rand.SCCWRP,orders=1:10,plot=FALSE,rescale=FALSE)
    ExpInterceptFFGRand <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    ExpExpFFGRand <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    ExpAICFFGRand <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    PLInterceptFFGRand <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
    PLExpFFGRand <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    PLAICFFGRand <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    zeta_1 <- Zeta.order.ex(ffg.rand.SCCWRP,order=1,rescale=FALSE) #Zeta order 1 for FFGs.
    MeanAlphaFFGRand <- zeta_1$zeta.val #Mean alpha diversity for FFGs.
    SDAlphaFFGRand <- zeta_1$zeta.val.sd #Standard deviation on alpha diversity for FFGs.
    zeta_2 <- Zeta.order.ex(ffg.rand.SCCWRP,order=2,rescale=FALSE) #Zeta order 2 for FFGs.
    BetaFFGRand <- zeta_2$zeta.val/(2*zeta_1$zeta.val-zeta_2$zeta.val) #Beta diversity for FFGs based on Jaccard dissimilarity.
    row <- t(as.data.frame(c(WS,mean(metadata$CSCI),mean(metadata$LU_2000_5K),sd(metadata$LU_2000_5K),mean(metadata$altitude),sd(metadata$altitude),MeanAlpha,SDAlpha,Beta,ExpIntercept,ExpExp,ExpAIC,PLIntercept,PLExp,PLAIC,MeanAlphaFFG,SDAlphaFFG,BetaFFG,ExpInterceptFFG,ExpExpFFG,ExpAICFFG,PLInterceptFFG,PLExpFFG,PLAICFFG,MeanAlphaFFGRand,SDAlphaFFGRand,BetaFFGRand,ExpInterceptFFGRand,ExpExpFFGRand,ExpAICFFGRand,PLInterceptFFGRand,PLExpFFGRand,PLAICFFGRand)))
    zetaAnalysis <- rbind(zetaAnalysis,row)
    print(paste(j,row))
  }
}
colnames(zetaAnalysis) <- c("Watershed","MeanCSCI","MeanLU","SDLU","MeanAltitude","SDAltitude","MeanAlpha","SDAlpha","Beta","ExpIntercept","ExpExp","ExpAIC","PLIntercept","PLExp","PLAIC","MeanAlphaFFG","SDAlphaFFG","BetaFFG","ExpInterceptFFG","ExpExpFFG","ExpAICFFG","PLInterceptFFG","PLExpFFG","PLAICFFG","MeanAlphaFFGRand","SDAlphaFFGRand","BetaFFGRand","ExpInterceptFFGRand","ExpExpFFGRand","ExpAICFFGRand","PLInterceptFFGRand","PLExpFFGRand","PLAICFFGRand")
zetaAnalysis[,1:ncol(zetaAnalysis)] <- sapply(zetaAnalysis[,1:ncol(zetaAnalysis)],as.character)
zetaAnalysis[,2:ncol(zetaAnalysis)] <- sapply(zetaAnalysis[,2:ncol(zetaAnalysis)],as.numeric)
rownames(zetaAnalysis) <- 1:nrow(zetaAnalysis)
#write.table(zetaAnalysis,"zetaAnalysis100PermutationsHUC8v2.txt",quote=FALSE,sep="\t",row.names = FALSE)
zetaAnalysis <- read.table("zetaAnalysis100PermutationsHUC8.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#zetaModel <- lm(MeanCSCI ~ MeanAlpha+MeanAlphaFFG+Beta+BetaFFG+PLExp+PLExpFFG, data=zetaAnalysis)
zetaModel <- lm(MeanCSCI ~ MeanAlpha+MeanAlphaFFG+Beta+BetaFFG+PLExpFFG, data=zetaAnalysis)
summary(zetaModel)
zetaModel$coefficients
anova(zetaModel)
calc.relimp(zetaModel,type="lmg",rela=TRUE)

#calculate a linear model fit.
#zetaAnalysis$fit1 <- zetaModel$coefficients[1]+zetaAnalysis$MeanAlpha*zetaModel$coefficients[2]+zetaAnalysis$MeanAlphaFFG*zetaModel$coefficients[3]+zetaAnalysis$Beta*zetaModel$coefficients[4]+zetaAnalysis$BetaFFG*zetaModel$coefficients[5]+zetaAnalysis$PLExp*zetaModel$coefficients[6]+zetaAnalysis$PLExpFFG*zetaModel$coefficients[7]
zetaAnalysis$fit1 <- zetaModel$coefficients[1]+zetaAnalysis$MeanAlpha*zetaModel$coefficients[2]+zetaAnalysis$MeanAlphaFFG*zetaModel$coefficients[3]+zetaAnalysis$Beta*zetaModel$coefficients[4]+zetaAnalysis$BetaFFG*zetaModel$coefficients[5]+zetaAnalysis$PLExpFFG*zetaModel$coefficients[6]

plot(zetaAnalysis$MeanCSCI,zetaAnalysis$fit1)
cor.test(zetaAnalysis$MeanCSCI,zetaAnalysis$fit1,method="pearson")

#Regression between network parameters.
library(Hmisc)
library(corrplot)
library("PerformanceAnalytics")
#Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
chart.Correlation(zetaAnalysis[,c("MeanCSCI","MeanAltitude","MeanLU","MeanAlpha","MeanAlphaFFG","SDAlpha","SDAlphaFFG","Beta","BetaFFG","PLExp","PLExpFFG","fit1")], histogram=FALSE, method="pearson")
