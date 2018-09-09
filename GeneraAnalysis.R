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
#Get unique genera for the statewide data set.
GeneraCounts <- na.omit(as.data.frame(unique(GISBioData$FinalID)))
colnames(GeneraCounts) <- c("FinalID")
tmp <- GISBioData[,c("UniqueID","FinalID","Count")]

#Construct a table of counts by genera in a Phyloseq otu table format
otudata <- arrange(GeneraCounts,FinalID)
for(ID in unique(GISBioData$UniqueID)){
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
#Samples are comprised of the relative abundances of taxa by genera.
set.seed(1)
test <- subset_samples(physeq)
test <- transform_sample_counts(test, function(x) x/sum(x))
testDF <- as(sample_data(test), "data.frame")
testAdonis <- adonis(distance(test,method="bray") ~ Watershed+altitude+LU_2000_5K, data=testDF, permutations=1000, strata=testDF$Year)
testAdonis
