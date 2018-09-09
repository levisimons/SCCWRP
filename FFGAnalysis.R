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
library(devtools)
library(microbiomeSeq)

#Analysis of trends involving functional feeding groups per sample in the SCCWRP data set.
setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage values.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Remove duplicate measures.
GISBioData <- GISBioData[!duplicated(GISBioData[,c("UniqueID","FinalID","Count")]),]
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Merge in altitude and year readings.
tmp <- SCCWRP[,c("UniqueID","altitude","Year")]
GISBioData <- join(GISBioData,tmp,by=c("UniqueID"))
#Determine taxonomic richness per sample
tmp <- data.frame(table(GISBioData$UniqueID))
colnames(tmp) <- c("UniqueID","nTaxa")
GISBioData <- join(GISBioData,tmp,by=c("UniqueID"))
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
#Determine relative abundance of funcional feeding groups
GISBioData$FFGRA <- GISBioData$FFGCount/GISBioData$ActualOrganismCount

#Construct a table of counts by functional feeding group in a Phyloseq otu table format
otudata <- arrange(FFGCounts,FunctionalFeedingGroup)
for(ID in unique(GISBioData$UniqueID)){
  tmp2 <- subset(tmp,UniqueID==ID)
  tmp2 <- tmp2[,c("FunctionalFeedingGroup","FFGCount")]
  colnames(tmp2) <- c("FunctionalFeedingGroup",ID)
  tmp2 <- merge(tmp2,FFGCounts,by=c("FunctionalFeedingGroup"),all=TRUE)
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
taxmat <- as.matrix(FFGCounts)
rownames(taxmat) <- as.data.frame(taxmat)$FunctionalFeedingGroup
TAX = tax_table(taxmat)
samplemat <- GISBioData[,c("UniqueID","Watershed","LU_2000_5K","Year","altitude")]
samplemat <- samplemat[!duplicated(samplemat),]
#samplemat$Watershed <- as.factor(samplemat$Watershed)
samplemat$LU_2000_5K <- as.numeric(as.character(samplemat$LU_2000_5K))
samplemat$Year <- as.numeric(as.character(samplemat$Year))
samplemat$altitude <- as.numeric(as.character(samplemat$altitude))
row.names(samplemat) <- samplemat$UniqueID
sampledata <- sample_data(samplemat)
physeq <- phyloseq(OTU,TAX,sampledata)

#What factors are driving variations in beta diversity?
#Samples are comprised of the relative abundances of taxa by functional feeding groups.
set.seed(1)
test <- subset_samples(physeq)
test <- transform_sample_counts(test, function(x) x/sum(x))
testDF <- as(sample_data(test), "data.frame")
testAdonis <- adonis(distance(test,method="bray") ~ Watershed+altitude+LU_2000_5K, data=testDF, permutations=1000, strata=testDF$Year)
testAdonis

#Determine land use deciles for the full state data set.
LUquantile <- quantile(GISBioData$LU_2000_5K,probs=seq(0,1,0.1))
#Determine altitude deciles for the full state data set.
Altitudequantile <- quantile(GISBioData$altitude,probs=seq(0,1,0.1))
#Compute trends in beta diversity groups of samples set by
#altitude, year, and land use.
FFGTrends <- data.frame()
for(year in unique(GISBioData$Year)){
  annualSubset <- subset(GISBioData,Year==year)
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
          sampleDF <- subset(GISBioData,Year==year & LU_2000_5K >= LULow & LU_2000_5K <= LUHigh & altitude >= AltitudeLow & altitude <= AltitudeHigh)
        }
        if(length(unique(sampleDF$UniqueID))>2){
          physeqLUYear <- subset_samples(physeq,Year==year & LU_2000_5K >= LULow & LU_2000_5K <= LUHigh & altitude >= AltitudeLow & altitude <= AltitudeHigh)
          if(nrow(otu_table(physeqLUYear)) > 2){
            sampleBeta <- distance(physeqLUYear,method="bray")
            print(paste(year,MidLU,MidAltitude,length(unique(sampleDF$UniqueID)),mean(sampleBeta),sd(sampleBeta)))
            row <- t(as.data.frame(c(year,MidLU,MidAltitude,length(unique(sampleDF$UniqueID)),mean(sampleBeta),sd(sampleBeta))))
            FFGTrends <- rbind(row,FFGTrends)
          }
        }
      }
    }
  }
}
colnames(FFGTrends) <- c("Year","MidLU","MidAltitude","NSamples","MeanBeta","SDBeta")
rownames(FFGTrends) <- 1:nrow(FFGTrends)

#Plotting a NMDS plot of beta diversity of functional feeding group by sample data.
test.ord <- ordinate(test,"NMDS","bray")
p = plot_ordination(test,test.ord,type="samples", color="Year", title="Samples by functional feed group")+ theme(legend.position="none")
p

#Plot eigenvalues of PCA.
Dist = distance(test, method = "bray")
ord = ordinate(test, method = "PCoA", distance = Dist)
plot_scree(ord, "Scree Plot: Bray-Curtis MDS")

#Plot differences in bray-curtist dissimilarity with differences in environmental parameters.
distmat <- as.matrix(Dist)
distBeta <- setNames(melt(distmat), c('Sample1', 'Sample2', 'BCdist'))
LUData <- GISBioData[,c("UniqueID","LU_2000_5K")]
LUData <- LUData[!duplicated(LUData),]
distBeta <- merge(distBeta,LUData,by.x=c("Sample1"),by.y=c("UniqueID"))
colnames(distBeta)[which(names(distBeta) == "LU_2000_5K")] <- "LU1"
distBeta <- merge(distBeta,LUData,by.x=c("Sample2"),by.y=c("UniqueID"))
colnames(distBeta)[which(names(distBeta) == "LU_2000_5K")] <- "LU2"
distBeta$deltaLU <- abs(distBeta$LU1-distBeta$LU2)
altitudeData <- GISBioData[,c("UniqueID","altitude")]
altitudeData <- altitudeData[!duplicated(altitudeData),]
distBeta <- merge(distBeta,altitudeData,by.x=c("Sample1"),by.y=c("UniqueID"))
colnames(distBeta)[which(names(distBeta) == "altitude")] <- "altitude1"
distBeta <- merge(distBeta,altitudeData,by.x=c("Sample2"),by.y=c("UniqueID"))
colnames(distBeta)[which(names(distBeta) == "altitude")] <- "altitude2"
distBeta$deltaAltitude <- abs(distBeta$altitude1-distBeta$altitude2)

#Plotting trends in functional feeding group beta diversity
plot(FFGTrends$MidLU,FFGTrends$MeanBeta)
lines(smooth.spline(FFGTrends$MidLU,FFGTrends$MeanBeta,df=3),col="red")
cor.test(FFGTrends$MidLU,FFGTrends$MeanBeta,method="spearman")
