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

#Read in insect data from CEDEN sites.
insectDataCEDENRAW <- read.table("BugTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
insectDataCEDEN <- filter(insectDataCEDENRAW, CollectionReplicate==1)
insectTaxaCEDEN <- insectDataCEDEN[,c("Phylum","Class","Orders","Family","Genus","FinalID")]
names(insectTaxaCEDEN)[names(insectTaxaCEDEN)=="Orders"]<-"Order"
#Get higher order taxanomic IDs for CEDEN data.
CEDENHighTaxa <- insectTaxaCEDEN[,c("Phylum","Class","Order")]

#Read in SMC site data.
SMCDataRaw <- read.table("SMCAllSites.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
SMCData <- SMCDataRaw[,c("Order","Family","Genus","FinalID")]
SMCData <- SMCData[!duplicated(SMCData$FinalID),]
SMCTaxa <- join(SMCData,CEDENHighTaxa,by=c("Order"))
SMCTaxa <- SMCTaxa[!duplicated(SMCTaxa$FinalID),]
SMCTaxa <- SMCTaxa[,c("Phylum","Class","Order","Family","Genus","FinalID")]

#Merge to make a unified taxonomy table
taxa <- rbind(SMCTaxa,insectTaxaCEDEN)
taxa <- taxa[match(unique(taxa$FinalID), taxa$FinalID),]
taxa <- taxa[,c("Phylum","Class","Order","Family","Genus","FinalID")]
taxa[taxa==""] <- NA
taxa <- filter(taxa,taxa$FinalID!="NA")
taxaID <- unique(taxa$FinalID)

#Read in site data.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Get additional site metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
watershedTable <- as.data.frame(table(SCCWRP$Watershed))
#Merge in watershed IDs into site data.
GISBioData <- join(GISBioData,SCCWRP[,c("UniqueID","Watershed")],by=c("UniqueID"))

#Select site data which has known taxonomic data.
BMIData <- subset(GISBioData,GISBioData$FinalID %in% taxaID)

#Subset site data by watershed.
BMIData <- subset(BMIData,Watershed=="Cuyama")

#Format taxa IDs for downstream Phyloseq work and network generation.
OTUID <- as.data.frame(taxaID)
colnames(OTUID) <- c("FinalID")

#To generate an OTU table of SCCWRP data for use in Phyloseq.
for(sample in unique(BMIData$UniqueID)){
  tmp1 <- subset(BMIData,UniqueID==sample)
  tmp1 <- tmp1[!duplicated(tmp1[c("FinalID","Count")]),]
  tmp1 <- join(OTUID,tmp1,by=c("FinalID"))[,c("FinalID","Count")]
  OTUID <- cbind(OTUID,tmp1$Count)
  names(OTUID)[names(OTUID)=="tmp1$Count"]<-sample
}

#Create Phyloseq object with the OTU table, sample factors, and taxonomic data.
otumat <- as.matrix(OTUID[,-c(1)])
otumat[is.na(otumat)] <- 0
rownames(otumat) <- OTUID$FinalID
OTU = otu_table(otumat,taxa_are_rows = TRUE)
taxmat <- as.matrix(subset(taxa,taxa$FinalID %in% OTUID$FinalID))
rownames(taxmat) <- as.data.frame(taxmat)$FinalID
taxmat <- taxmat[,-c(6)]
TAX = tax_table(taxmat)
samplemat <- as.matrix(SCCWRP)
row.names(samplemat) <- SCCWRP$UniqueID
sampledata <- sample_data(as.data.frame(samplemat))
physeq <- phyloseq(OTU,TAX,sampledata)
test <- physeq

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
testAdonis = adonis(distance(physeqSubset,method="bray")~Year+Watershed+Ag_2000_5K+CODE_21_2000_5K+URBAN_2000_5K+LU_2000_5K+altitude,data=testDF,permutations = 1000)
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

