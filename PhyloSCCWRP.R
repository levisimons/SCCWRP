library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)
library("phyloseq")
library(lattice)

setwd("~/Desktop/SCCWRP")

#The following files are read in to generate a single SCCWRP taxonomy file.

#Read in algae data from CEDEN sites.
algaeDataCEDENRaw <- read.table("AlgaeTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
algaeDataCEDEN <- filter(algaeDataCEDENRaw, CollectionReplicate==1)
algaeTaxaCEDEN <- algaeDataCEDEN[,c("Phylum","Class","Orders","Family","Genus","FinalID")]
names(algaeTaxaCEDEN)[names(algaeTaxaCEDEN)=="Orders"]<-"Order"

#Read in insect data from CEDEN sites.
insectDataCEDENRAW <- read.table("BugTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
insectDataCEDEN <- filter(insectDataCEDENRAW, CollectionReplicate==1)
insectTaxaCEDEN <- insectDataCEDEN[,c("Phylum","Class","Orders","Family","Genus","FinalID")]
names(insectTaxaCEDEN)[names(insectTaxaCEDEN)=="Orders"]<-"Order"

#Merge to make a unified taxonomy table
CEDENTaxa <- rbind(algaeTaxaCEDEN,insectTaxaCEDEN)
CEDENHighTaxa <- CEDENTaxa[,c("Phylum","Class","Order")]

SMCDataRaw <- read.table("SMCAllSites.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
SMCData <- SMCDataRaw[,c("Order","Family","Genus","FinalID")]
SMCData <- SMCData[!duplicated(SMCData$FinalID),]
SMCTaxa <- join(SMCData,CEDENHighTaxa,by=c("Order"))
SMCTaxa <- SMCTaxa[!duplicated(SMCTaxa$FinalID),]
SMCTaxa <- SMCTaxa[,c("Phylum","Class","Order","Family","Genus","FinalID")]

#Merge to make a unified taxonomy table
taxa <- rbind(SMCTaxa,CEDENTaxa)
taxa <- taxa[match(unique(taxa$FinalID), taxa$FinalID),]
taxa <- taxa[,c("Phylum","Class","Order","Family","Genus","FinalID")]
taxa[taxa==""] <- NA
taxa <- filter(taxa,taxa$FinalID!="NA")
taxaID <- unique(taxa$FinalID)

#Read in aggregated biological data organized by taxonomic counts.
#This step is to generate an OTU table for Phyloseq
bioDataRaw <- read.table("CombinedTaxonomy.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
bioData <- filter(bioDataRaw, Replicate==1)
#Filter out entries without a FinalID
bioData <- filter(bioData,FinalID!="")
#Filter out entries without a full taxonomy
bioData <- subset(bioData,bioData$FinalID %in% taxaID)

#Read in older merged data set to get unique algal and invertebrate taxa IDs.
oldSet <- read.table("GISBioDataSmallSet.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)

#Separate algae counts.
algae <- filter(oldSet,MeasurementType=="Soft-bodied algal relative abundance" | MeasurementType=="Benthic algal relative abundance")
algaeID <- unique(algae$FinalID)
algaeData <- subset(bioData,bioData$FinalID %in% algaeID)
#Separate out benthic algae counts.
benthicAlgaeData <- algaeData[!is.na(algaeData$BAResult),]
benthicAlgaeData <- benthicAlgaeData[,c("StationCode","SampleDate","FinalID","BAResult")]
benthicAlgaeID <- unique(benthicAlgaeData$FinalID)
#Separate out soft-bodied algae counts.
softAlgaeData <- algaeData[!is.na(algaeData$Result),]
softAlgaeData <- softAlgaeData[,c("StationCode","SampleDate","FinalID","Result")]
softAlgaeID <- unique(softAlgaeData$FinalID)

#Separate benthic macroinvertebrate counts.
BMI <- filter(oldSet,MeasurementType=="Invertebrate relative abundance" | MeasurementType=="Invertebrate relative abundances")
BMIID <- unique (BMI$FinalID)
BMIData <- subset(bioData,bioData$FinalID %in% BMIID)
BMIData <- BMIData[,c("StationCode","SampleDate","FinalID","BAResult")]
tmp1<-data.table(BMIData)
#Sum insect counts for matching IDs, but with different life stages.
BMIData<-tmp1[,.(BAResult=sum(BAResult)),by=.(StationCode,SampleDate,FinalID)]

#Determine the benthic algae totals count column and make it a temporary dataframe.
tmp1<-data.table(benthicAlgaeData)
tmp1<-tmp1[,.(BAResult=sum(BAResult)),by=.(StationCode,SampleDate)]
colnames(tmp1) <- c("StationCode","SampleDate","ActualOrganismCount")
#Join the total benthic algal count by sample to each unique sample.
benthicAlgaeData <- join(benthicAlgaeData,tmp1,by=c("StationCode","SampleDate"))
#Calculate the relative abundance of benthic algae data.
benthicAlgaeData$Measurement <- with(benthicAlgaeData,BAResult/ActualOrganismCount)
#Add measurement type label.
benthicAlgaeData$MeasurementType <- "benthic alage relative abundance"
#Add unique ID per sample.
benthicAlgaeData$UniqueID <- paste(benthicAlgaeData$StationCode,benthicAlgaeData$SampleDate)
#Subset key columns.
benthicAlgaeData <- benthicAlgaeData[,c("StationCode","SampleDate","UniqueID","FinalID","MeasurementType","Measurement","BAResult","ActualOrganismCount")]
names(benthicAlgaeData)[names(benthicAlgaeData)=="BAResult"]<-"Count"

#Determine the soft-bodied algae totals count column and make it a temporary dataframe.
tmp1<-data.table(softAlgaeData)
tmp1<-tmp1[,.(Result=sum(Result)),by=.(StationCode,SampleDate)]
colnames(tmp1) <- c("StationCode","SampleDate","ActualOrganismCount")
#Join the total benthic algal count by sample to each unique sample.
softAlgaeData <- join(softAlgaeData,tmp1,by=c("StationCode","SampleDate"))
#Calculate the relative abundance of benthic algae data.
softAlgaeData$Measurement <- with(softAlgaeData,Result/ActualOrganismCount)
#Add measurement type label.
softAlgaeData$MeasurementType <- "soft-bodied alage relative abundance"
#Add unique ID per sample.
softAlgaeData$UniqueID <- paste(softAlgaeData$StationCode,softAlgaeData$SampleDate)
#Subset key columns.
softAlgaeData <- softAlgaeData[,c("StationCode","SampleDate","UniqueID","FinalID","MeasurementType","Measurement","Result","ActualOrganismCount")]
names(softAlgaeData)[names(softAlgaeData)=="Result"]<-"Count"

#Determine the benthic algae totals count column and make it a temporary dataframe.
tmp1<-data.table(BMIData)
tmp1<-tmp1[,.(BAResult=sum(BAResult)),by=.(StationCode,SampleDate)]
colnames(tmp1) <- c("StationCode","SampleDate","ActualOrganismCount")
#Join the total benthic macroinvertebrate count by sample to each unique sample.
BMIData <- join(BMIData,tmp1,by=c("StationCode","SampleDate"))
#Calculate the relative abundance of benthic algae data.
BMIData$Measurement <- with(BMIData,BAResult/ActualOrganismCount)
#Add measurement type label.
BMIData$MeasurementType <- "benthic macroinvertebrate relative abundance"
#Add unique ID per sample.
BMIData$UniqueID <- paste(BMIData$StationCode,BMIData$SampleDate)
#Subset key columns.
BMIData <- BMIData[,c("StationCode","SampleDate","UniqueID","FinalID","MeasurementType","Measurement","BAResult","ActualOrganismCount")]
names(BMIData)[names(BMIData)=="BAResult"]<-"Count"

#Select taxa group for downstream Phyloseq work and network generation.
otudata <- BMIData

#Read in CSCI and land usage data by sample.
GISDataRAW <- read.csv("CSCI.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
GISData <- subset(GISDataRAW,FieldReplicate==1)
GISData$UniqueID <- paste(GISData$StationCode,GISData$SampleDate)
GISData <- subset(GISData,GISData$UniqueID %in% unique(otudata$UniqueID))
GISData <- GISData[!duplicated(GISData$UniqueID),]

#Add in aggregated land coverage at a 5km watershed buffer.
GISData$LU_2000_5K <- with(GISData,Ag_2000_5K+CODE_21_2000_5K+URBAN_2000_5K)

#Final step to ensure that both the factors and biological data have the same set of sample IDs.
otudata <- subset(otudata,otudata$UniqueID %in% unique(GISData$UniqueID))

#Generate an OTU table of SCCWRP data for use in Phyloseq.
OTUID <- as.data.frame(BMIID)
colnames(OTUID) <- c("FinalID")
for(sample in unique(otudata$UniqueID)){
  tmp1 <- subset(otudata,UniqueID==sample)
  tmp1 <- join(OTUID,tmp1,by=c("FinalID"))[,c("FinalID","Count")]
  OTUID <- cbind(OTUID,tmp1$Count)
  names(OTUID)[names(OTUID)=="tmp1$Count"]<-sample
  print(sample)
}

#Create Phyloseq object with the OTU table, sample factors, and taxonomic data.
otumat <- as.matrix(OTUID[,-c(1)])
otumat[is.na(otumat)] <- 0
rownames(otumat) <- OTUID$FinalID
OTU = otu_table(otumat,taxa_are_rows = TRUE)
taxmat <- as.matrix(subset(taxa,taxa$FinalID %in% BMIID))
rownames(taxmat) <- as.data.frame(taxmat)$FinalID
taxmat <- taxmat[,-c(6)]
TAX = tax_table(taxmat)
samplemat <- as.matrix(GISData)
row.names(samplemat) <- GISData$UniqueID
sampledata <- sample_data(as.data.frame(samplemat))
physeq <- phyloseq(OTU,TAX,sampledata)

test<-subset_samples(physeq,SMCShed=="SanGabriel")
#test <- transform_sample_counts(test, function(x) x/sum(x))
test<-tax_glom(test,taxrank=rank_names(test)[2],NArm=TRUE,bad_empty=NA)

Dist = distance(test, method = "chao")
ord = ordinate(test, method = "PCoA", distance = Dist)
plot_scree(ord, "Scree Plot: Bray-Curtis MDS")

levelplot(as.matrix(Dist))
hist(as.vector(as.matrix(Dist)))
#Write out merged data set to read back in the future as opposed to 
#generating it each time from component data files.
#write.csv(GISBioData,file="CAGISBioData.csv",row.names=FALSE)
