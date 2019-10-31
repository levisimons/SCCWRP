rm(list=ls())
require("plyr")
require(dplyr)
require(zetadiv)
require(sp)
require(rgdal)
require(geosphere)

#setwd("~/Desktop/SCCWRP/Metagenomics")
setwd("/home/cmb-07/sn1/alsimons/SCCWRP")

#Read in metagenomic count tables and format them as presence/absence tables.
communityInputRawPlate1 <- read.table("18SV9P1TableWithTaxonomy.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
communityInputRawPlate2 <- read.table("18SV9P2TableWithTaxonomy.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Remove spaces from column names for the count tables.
communityInputRawPlate1 <- rename(communityInputRawPlate1, OTUID = `OTU ID`)
communityInputRawPlate2 <- rename(communityInputRawPlate2, OTUID = `OTU ID`)

#Get a list of unique OTU names from count tables and convert to a data frame.
uniqueOTUs <- as.data.frame(unique(c(communityInputRawPlate1$OTUID,communityInputRawPlate2$OTUID)))
colnames(uniqueOTUs) <- c("OTUID")
uniqueOTUs$OTUID <- as.character(uniqueOTUs$OTUID)

#Create a merged metagenomic count table.
communityInput <- left_join(uniqueOTUs,communityInputRawPlate1,by="OTUID")
communityInput <- left_join(communityInput,communityInputRawPlate2,by="OTUID")

#Remove unnecessary columns.
communityInput <- communityInput[, -which(names(communityInput)  %in% c("DNAStandard","Ext-Blank1","Ext-Blank2","FB","NTC","SNAStandardII","ConsensusLineage.y","DNAstandardI","DNAstandardII","ConsensusLineage.x","202","230.x","230.y"))]

#Convert abundance to presence/absence.
rownames(communityInput) <- communityInput$OTUID
communityInput$OTUID <- NULL
communityInput[is.na(communityInput)] <- 0
communityInput[communityInput > 0] <- 1

#Read in table linking sample IDs in the metagenomic table to sample station codes.
sampleIDs <- read.table("SampleStationCodesID.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Remove dubious samples
sampleIDs <- sampleIDs[which(sampleIDs$SampleNum!=202 & sampleIDs$SampleNum!=230),]

#Read in environmental metadata by station codes table.
metadata <- read.table("MetagenomicSampleSiteMetadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Merge in sample IDs to the metadata.
metadata <- merge(sampleIDs,metadata)

#Subset metadata by which samples are present in your community data set.
metadata <- metadata[which(metadata$SampleNum %in% as.numeric(colnames(communityInput))),]
#Sort metadata by sample number.
metadata <- arrange(metadata,SampleNum)
metadata$elev_range <- as.numeric(metadata$elev_range)
metadata$max_elev <- as.numeric(metadata$max_elev)
#Create aggregate upstream land use variable.
metadata$LU <- metadata$Ag_2011_5K+metadata$URBAN_2011_5K+metadata$CODE_21_2011_5K

#Read in nitrogen and phosporus site data.
NPdata <- read.table("NandP_labdata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
NPdata[NPdata=="ND" | NPdata==""] <- NA
NPdata[,2:4] <- sapply(NPdata[,2:4],as.numeric)
NPdata[NPdata<0] <- NA

#Merge in nitrogen and phosphorus site data into metada.
metadata <- left_join(metadata,NPdata,by=c("StationCode"))

#Read in watershed by sample site location data.
SCCWRP <- read.table("MetagenomicSitesWithWatersheds.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
SCCWRP <- SCCWRP[,c("Latitude","Longitude","HUC8")]
#Remove duplicate rows
SCCWRP <- SCCWRP[!duplicated(SCCWRP),]
#Create HUC 6 watershed column.
SCCWRP$HUC6 <- substr(SCCWRP$HUC8,1,nchar(SCCWRP$HUC8)-2)
#Create HUC 4 watershed column.
SCCWRP$HUC4 <- substr(SCCWRP$HUC8,1,nchar(SCCWRP$HUC8)-4)

#Merge in watershed data into metadata
metadata <- left_join(metadata,SCCWRP,by=c("Latitude","Longitude"))
metadata <- metadata[!is.na(metadata$Latitude),]
metadata <- metadata[!is.na(metadata$Longitude),]

#Assign a geographic cluster to sample sites in the metadata
mdist <- distm(metadata[,c("Longitude","Latitude")])
km <- kmeans(mdist,centers=4)
metadata$clust <- km$cluster

#Reorder community data set columns so that their sample number order increases from left to right.
#This will match the order of the metadata data frame where sample number order increases going down.
communityInput <- communityInput[,as.character(metadata$SampleNum)]
#Create a species/site matrix for use in the zeta diversity functions.
data.spec <- as.data.frame(t(communityInput))

#Calculate how much zeta diversity of a particular order decays with distance
zetaDistance <- Zeta.ddecay(xy=metadata[,c("Latitude","Longitude")],data.spec=data.spec,order=4,distance.type="ortho",normalize=TRUE,plot=FALSE)
#How strongly correlated is zeta diversity with geographic distance?
cor.test(zetaDistance$zeta.val,zetaDistance$distance,method="spearman")
#
