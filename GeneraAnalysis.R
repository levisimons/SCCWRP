library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)
library(ggplot2)
library(vegan)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(reshape2)

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

selected <- GISBioDataLWS

#Get unique genera for the heavily sampled statewide data set.
GeneraCounts <- na.omit(as.data.frame(unique(selected$FinalID)))
colnames(GeneraCounts) <- c("FinalID")
tmp <- selected[,c("UniqueID","FinalID","Count")]

#Construct a table of counts by genera in a Phyloseq otu table format
otudata <- arrange(GeneraCounts,FinalID)
for(ID in unique(selected$UniqueID)){
  tmp2 <- subset(tmp,UniqueID==ID)
  tmp2 <- tmp2[,c("FinalID","Count")]
  colnames(tmp2) <- c("FinalID",ID)
  tmp2 <- merge(tmp2,GeneraCounts,by=c("FinalID"),all=TRUE)
  otudata <- cbind(otudata,tmp2)
  otudata <- otudata[,!duplicated(colnames(otudata))]
}

#Create Phyloseq object with the OTU table, sample factors, and taxonomic data.
otumat <- t(otudata)
colnames(otumat) <- otumat[c(1),]
otumat <- otumat[-c(1),]
otumat[is.na(otumat)] <- 0
otumat <- as.data.frame(otumat)
sampleNames <- rownames(otumat)
otumat <- sapply(otumat,as.numeric)
rownames(otumat) <- sampleNames
OTU = otu_table(otumat,taxa_are_rows = FALSE)
taxmat <- as.matrix(GeneraCounts)
rownames(taxmat) <- as.data.frame(taxmat)$FinalID
TAX = tax_table(taxmat)
samplemat <- selected[,c("UniqueID","HUC2","HUC4","HUC6","HUC8","LU_2000_5K","Year","altitude")]
samplemat <- samplemat[!duplicated(samplemat),]
samplemat$HUC2 <- as.factor(samplemat$HUC2)
samplemat$HUC4 <- as.factor(samplemat$HUC4)
samplemat$HUC6 <- as.factor(samplemat$HUC6)
samplemat$HUC8 <- as.factor(samplemat$HUC8)
samplemat$LU_2000_5K <- as.numeric(as.character(samplemat$LU_2000_5K))
samplemat$Year <- as.factor(samplemat$Year)
samplemat$altitude <- as.numeric(as.character(samplemat$altitude))
row.names(samplemat) <- samplemat$UniqueID
sampledata <- sample_data(samplemat)
physeq <- phyloseq(OTU,TAX,sampledata)

#What factors are driving variations in alpha and beta diversity?
#Samples are comprised of the relative abundances of taxa by genera.
set.seed(1)
test <- subset_samples(physeq)
#test <- transform_sample_counts(test, function(x) x/sum(x))
testDF <- as(sample_data(test), "data.frame")
testDF$HUC2 <- as.factor(testDF$HUC2)
testDF$HUC4 <- as.factor(testDF$HUC4)
testDF$HUC6 <- as.factor(testDF$HUC6)
testDF$HUC8 <- as.factor(testDF$HUC8)
levels(alphaDF$HUC8) <- unique(alphaDF$HUC8)
alphaDF <- GISBioData[,c("UniqueID","HUC2","HUC4","HUC6","HUC8","altitude","LU_2000_5K","Year","nTaxa","CSCI")]
alphaDF <- alphaDF[!duplicated(alphaDF),]
alphaDF$HUC2 <- as.factor(alphaDF$HUC2)
alphaDF$HUC4 <- as.factor(alphaDF$HUC4)
alphaDF$HUC6 <- as.factor(alphaDF$HUC6)
alphaDF$HUC8 <- as.factor(alphaDF$HUC8)
levels(alphaDF$HUC8) <- unique(alphaDF$HUC8)
summary(aov(nTaxa~HUC2+HUC4+HUC6+HUC8+altitude+LU_2000_5K+Year,data=alphaDF))
testBeta <- adonis(distance(test,method="jaccard") ~ HUC2+HUC4+HUC6+HUC8+altitude+LU_2000_5K+Year, data=testDF, permutations=1000)
testBeta

#Determine land use deciles for the full state data set.
LUdf <- GISBioData[,c("UniqueID","LU_2000_5K")]
LUdf <- LUdf[!duplicated(LUdf),]
LUquantile <- quantile(LUdf$LU_2000_5K,probs=seq(0,1,0.1))
#Determine altitude deciles for the full state data set.
Altitudedf <- GISBioData[,c("UniqueID","altitude")]
Altitudedf <- Altitudedf[!duplicated(Altitudedf),]
Altitudequantile <- quantile(Altitudedf$altitude,probs=seq(0,1,0.1))
#Compute trends in beta diversity groups of samples set by
#altitude, year, and land use.
GeneraTrends <- data.frame()
for(i in 1:length(LUquantile)){
  if(i>1){
    LULow <- LUquantile[i-1]
    LUHigh <- LUquantile[i]
    MidLU <- 0.5*(LULow+LUHigh)
    for(j in 1:length(Altitudequantile)){
      if(j>1){
        AltitudeLow <- Altitudequantile[j-1]
        AltitudeHigh <- Altitudequantile[j]
        MidAltitude <- 0.5*(AltitudeLow+AltitudeHigh)
        sampleDF <- subset(GISBioData,LU_2000_5K >= LULow & LU_2000_5K <= LUHigh & altitude >= AltitudeLow & altitude <= AltitudeHigh)
        if(length(unique(sampleDF$UniqueID))>2){
          physeqLUYear <- subset_samples(physeq,LU_2000_5K >= LULow & LU_2000_5K <= LUHigh & altitude >= AltitudeLow & altitude <= AltitudeHigh)
          if(nrow(otu_table(physeqLUYear)) > 2){
            sampleBeta <- distance(physeqLUYear,method="bray")
            print(paste(MidLU,MidAltitude,length(unique(sampleDF$UniqueID)),mean(sampleBeta),sd(sampleBeta)))
            row <- t(as.data.frame(c(MidLU,MidAltitude,length(unique(sampleDF$UniqueID)),mean(sampleBeta),sd(sampleBeta))))
            GeneraTrends <- rbind(row,GeneraTrends)
          }
        }
      }
    }
  }
}
colnames(GeneraTrends) <- c("MidLU","MidAltitude","NSamples","MeanBeta","SDBeta")
rownames(GeneraTrends) <- 1:nrow(GeneraTrends)
