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
library(reshape2)
library(geosphere)
library(RADanalysis)
library(MASS)
library(fitdistrplus)

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
#Add unique ID per sample.
benthicAlgaeData$UniqueID <- paste(benthicAlgaeData$StationCode,benthicAlgaeData$SampleDate)
benthicAlgaeData <- benthicAlgaeData[,c("StationCode","SampleDate","FinalID","BAResult",'UniqueID')]
benthicAlgaeID <- unique(benthicAlgaeData$FinalID)
#Separate out soft-bodied algae counts.
softAlgaeData <- algaeData[!is.na(algaeData$Result),]
#Add unique ID per sample.
softAlgaeData$UniqueID <- paste(softAlgaeData$StationCode,softAlgaeData$SampleDate)
softAlgaeData <- softAlgaeData[,c("StationCode","SampleDate","FinalID","Result","UniqueID")]
softAlgaeID <- unique(softAlgaeData$FinalID)

#Determine the benthic algae totals count column and make it a temporary dataframe.
tmp1<-data.table(benthicAlgaeData)
tmp1<-tmp1[,.(BAResult=sum(BAResult)),by=.(UniqueID,FinalID)]
benthicAlgaeData <- tmp1
names(benthicAlgaeData)[names(benthicAlgaeData)=="BAResult"]<-"Count"

#Determine the soft-bodied algae totals count column and make it a temporary dataframe.
tmp1<-data.table(softAlgaeData)
tmp1<-tmp1[,.(Result=sum(Result)),by=.(UniqueID,FinalID)]
softAlgaeData <- tmp1
names(softAlgaeData)[names(softAlgaeData)=="Result"]<-"Count"

#Separate benthic macroinvertebrate counts.
BMI <- filter(oldSet,MeasurementType=="Invertebrate relative abundance" | MeasurementType=="Invertebrate relative abundances")
BMIID <- unique (BMI$FinalID)
BMIData <- subset(bioData,bioData$FinalID %in% BMIID)
#Add unique ID per sample.
BMIData$UniqueID <- paste(BMIData$StationCode,BMIData$SampleDate)
BMIData <- BMIData[,c("UniqueID","FinalID","BAResult")]
#Determine the benthic invertebrate totals count column and make it a temporary dataframe.
tmp1<-data.table(BMIData)
tmp1<-tmp1[,.(BAResult=sum(BAResult)),by=.(UniqueID,FinalID)]
BMIData <- tmp1
names(BMIData)[names(BMIData)=="BAResult"]<-"Count"

#Select taxa group for downstream Phyloseq work and network generation.
otudata <- BMIData
OTUID <- as.data.frame(BMIID)
colnames(OTUID) <- c("FinalID")

#Read in CSCI and land usage data by sample.
GISDataRAW <- read.csv("CSCI.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
GISData <- subset(GISDataRAW,FieldReplicate==1)
GISData$UniqueID <- paste(GISData$StationCode,GISData$SampleDate)
GISData <- subset(GISData,GISData$UniqueID %in% unique(otudata$UniqueID))
GISData <- GISData[!duplicated(GISData$UniqueID),]

#Add in aggregated land coverage at a 5km watershed buffer.
GISData$LU_2000_5K <- with(GISData,as.double(Ag_2000_5K+CODE_21_2000_5K+URBAN_2000_5K))
#Aggregate land usage into deciles for downstream analysis.
GISData$LU_Decile <- with(GISData, as.integer(LU_2000_5K/10))
#Add in sample year as a variable
GISData$Year <- year(mdy(GISData$SampleDate))

#Final step to ensure that both the factors and biological data have the same set of sample IDs.
otudata <- subset(otudata,otudata$UniqueID %in% unique(GISData$UniqueID))

#To generate an OTU table of SCCWRP data for use in Phyloseq.
for(sample in unique(otudata$UniqueID)){
  tmp1 <- subset(otudata,UniqueID==sample)
  tmp1 <- join(OTUID,tmp1,by=c("FinalID"))[,c("FinalID","Count")]
  OTUID <- cbind(OTUID,tmp1$Count)
  names(OTUID)[names(OTUID)=="tmp1$Count"]<-sample
}

#Generate diversity metrics per sample.
library("gambin")
library("sads")
siteDiversity <- as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(siteDiversity) <- c("UniqueID","gambinAlpha","Simpson","InvSimpson","nTaxa","Shannon","fisherAlpha","geomP","broken_stick_likelihood")
tmp2 <- as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(tmp2) <- c("UniqueID","gambinAlpha","Simpson","InvSimpson","nTaxa","Shannon","fisherAlpha","geomP","broken_stick_likelihood")
for(sample in unique(otudata$UniqueID)){
  tmp1 <- subset(otudata,UniqueID==sample)
  tmp2$UniqueID <- sample
  #Rarefy samples to 500 counts.  Permute and take the mean.
  RAD <- matrix(NA,nrow=1,ncol=nrow(tmp1))
  for(i in 1:20){
    tmp3 <- rrarefy(tmp1$Count,500)
    RAD <- mapply(c,as.data.frame(RAD),as.data.frame(tmp3))
  }
  RAD <- as.integer(colMeans(RAD+0.5,na.rm=TRUE))
  #To calculate the Gambin alpha parameter per sample
  tmp2$nTaxa <- nrow(tmp1)
  if(sum(RAD)>=500 & nrow(tmp1) >= 12 & min(RAD) > 0){
    gambin_fit <- fit_abundances(RAD)
    tmp2$gambinAlpha <- gambin_fit$alpha
    tmp2$Simpson <- diversity(RAD,index="simpson")
    tmp2$InvSimpson <- diversity(RAD,index="invsimpson")
    tmp2$Shannon <- diversity(RAD,index="shannon")
    tmp2$fisherAlpha <- fisher.alpha(RAD)
    geom_fit <- fitdist(RAD, distr ="geom")
    tmp2$geomP <- geom_fit$estimate
    broken_stick <- fitrad(RAD, rad="rbs")
    tmp2$broken_stick_likelihood <- -broken_stick@details$value
  }
  else{
    tmp2$gambinAlpha <- NA
    tmp2$Simpson <- NA
    tmp2$InvSimpson <- NA
    tmp2$Shannon <- NA
    tmp2$fisherAlpha <- NA
    tmp2$geomP <- NA
    tmp2$broken_stick_likelihood <- NA
  }
  siteDiversity <- rbind(siteDiversity,tmp2)
  print(sample)
}

siteDiversity <- siteDiversity[-c(1),]
GISData$gambinAlpha <- siteDiversity$gambinAlpha
GISData$Simpson <- siteDiversity$Simpson
GISData$InvSimpson <- siteDiversity$InvSimpson
GISData$nTaxa <- siteDiversity$nTaxa
GISData$Shannon <- siteDiversity$Shannon
GISData$fisherAlpha <- siteDiversity$fisherAlpha
GISData$geomP <- siteDiversity$geomP
GISData$broken_stick_likelihood <- siteDiversity$broken_stick_likelihood
write.csv(GISData,file="CSCI.csv",row.names=FALSE)

#Create Phyloseq object with the OTU table, sample factors, and taxonomic data.
otumat <- as.matrix(OTUID[,-c(1)])
otumat[is.na(otumat)] <- 0
rownames(otumat) <- OTUID$FinalID
OTU = otu_table(otumat,taxa_are_rows = TRUE)
taxmat <- as.matrix(subset(taxa,taxa$FinalID %in% OTUID$FinalID))
rownames(taxmat) <- as.data.frame(taxmat)$FinalID
taxmat <- taxmat[,-c(6)]
TAX = tax_table(taxmat)
samplemat <- as.matrix(GISData)
row.names(samplemat) <- GISData$UniqueID
sampledata <- sample_data(as.data.frame(samplemat))
sampledata$LU_2000_5K <- as.numeric(as.character(sampledata$LU_2000_5K))
sampledata$LU_Decile <- as.numeric(as.character(sampledata$LU_Decile))
physeq <- phyloseq(OTU,TAX,sampledata)

#Subset Phyloseq object by various factors and perform basic PCA and beta diversity tests.
#Unique watershed regions:
# "SMC_out"            "Ventura"            "SantaClara"         "SantaMonicaBay"    
# "SanGabriel"         "Calleguas"          "LosAngeles"         "MiddleSantaAna"    
# "UpperSantaAna"      "LowerSantaAna"      "SanJacinto"         "SanJuan"           
# "NorthernSanDiego"   "CentralSanDiego"    "MissionBaySanDiego" "SouthernSanDiego"
#physeqSubset <- subset_samples(physeq, SMCShed!="SMC_out")
physeqSubset <- subset_samples(physeq)
physeqSubset <- transform_sample_counts(physeqSubset, function(x) x/sum(x))
plot_ordination(physeqSubset,ordinate(physeqSubset,"NMDS","bray"),color="Watershed")

# Perform a PERMANOVA using a set number of permutations on a particular
# beta diversity metric and the significance of a particular design variable.
testDF = as(sample_data(physeqSubset), "data.frame")
testAdonis = adonis(distance(physeqSubset,method="bray")~Year+Watershed+LU_2000_5K+altitude,data=testDF,permutations = 1000)
testAdonis

#Plot beta diversity.
Dist = distance(physeqSubset, method = "bray")
ord = ordinate(physeqSubset, method = "PCoA", distance = Dist)
plot_scree(ord, "Scree Plot: Bray-Curtis MDS")
levelplot(as.matrix(Dist))
hist(as.vector(as.matrix(Dist)))

#Convert beta diversity matrix into a three column data frame.
distmat <- setNames(melt(as.matrix(Dist)), c("Site1","Site2","DiversityDistance"))
#Merge in the spatial coordinates for each pair of sites.
coord1 <- as.data.frame(as.matrix(sampledata[,c("UniqueID","Latitude","Longitude")]))
rownames(coord1) <- 1:nrow(coord1)
colnames(coord1) <- c("Site1","Lat1","Lon1")
distmat <- merge(distmat,coord1,by=c("Site1"))
coord2 <- as.data.frame(as.matrix(sampledata[,c("UniqueID","Latitude","Longitude")]))
rownames(coord2) <- 1:nrow(coord2)
colnames(coord2) <- c("Site2","Lat2","Lon2")
distmat <- merge(distmat,coord2,by=c("Site2"))
distmat$Lat1 <- as.numeric(as.character(distmat$Lat1))
distmat$Lon1 <- as.numeric(as.character(distmat$Lon1))
distmat$Lat2 <- as.numeric(as.character(distmat$Lat2))
distmat$Lon2 <- as.numeric(as.character(distmat$Lon2))
#Add in geodesic distances to compare separation distance to beta diversity distance.
distmat$SpatialDistance <- lapply(1:nrow(distmat),function(x) as.numeric((distm(as.matrix(distmat[x,c("Lon1","Lat1")]), as.matrix(distmat[x,c("Lon2","Lat2")]), fun = distGeo))))
distmat$SpatialDistance <- as.numeric(distmat$SpatialDistance)

#Regression between network parameters.
library(Hmisc)
library(corrplot)
library("PerformanceAnalytics")
chart.Correlation(GISData[,c("LU_2000_5K","altitude","nTaxa","Simpson","InvSimpson","gambinAlpha","fisherAlpha","Shannon","CSCI","geomP","broken_stick_likelihood")], histogram=FALSE, method="spearman")
chart.Correlation(GISData[,c("LU_2000_5K","altitude","nTaxa","CSCI")], histogram=FALSE, method="spearman")
chart.Correlation(GISData[,c("LU_2000_5K","gambinAlpha","CSCI")], histogram=FALSE, method="spearman")

#Generate map of data for a given chemical parameter in California.
library(ggmap)
library(maps)
library(mapdata)
dev.off()
MapCoordinates <- data.frame(GISData$LU_2000_5K,GISData$gambinAlpha,GISData$fisherAlpha,GISData$Simpson,GISData$InvSimpson,GISData$nTaxa,GISData$Shannon,GISData$CSCI,GISData$geomP,GISData$altitude,GISData$Longitude,GISData$Latitude,GISData$SMCShed)
colnames(MapCoordinates) = c("LU_2000_5K","gambinAlpha","fisherAlpha","Simpson","InvSimpson","nTaxa","Shannon","CSCI","geomP",'alt','lon','lat','SMCShed')
MapCoordinates <- na.omit(MapCoordinates)
mapBoundaries <- make_bbox(lon=MapCoordinates$lon,lat=MapCoordinates$lat,f=0.1)
CalMap <- get_map(location=mapBoundaries,maptype="satellite",source="google")
#CalMap <- ggmap(CalMap)+geom_point(data = MapCoordinates, mapping = aes(x = lon, y = lat, color = SMCShed),size=2)+scale_colour_gradientn(colours=rainbow(4))
CalMap <- ggmap(CalMap)+geom_point(data = MapCoordinates, mapping = aes(x = lon, y = lat, color = SMCShed),size=2)
CalMap

# Spiec-Easi network analysis.
library(devtools)
library(SpiecEasi)
library(phyloseq)
library(seqtime)
library(igraph)
library(network)
library(stringr)
spiec.out=spiec.easi(test, method="mb",icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(test)))
plot_network(spiec.graph, test, type='TAX', color="Red", label=NULL)
betaMat=as.matrix(symBeta(getOptBeta(spiec.out)))

library(maptools)
library(rgdal)
library(sf)
library(sp)
library(maps)
watersheds=readOGR("/Users/levisimons/Desktop/Data/wbdhu8_a_ca.shp")
coordinates(GISData) <- c("Latitude","Longitude")
inside.sheds <- !is.na(over(GISData, as(watersheds,"SpatialPolygons")))

