rm(list=ls())
require("plyr")
require(dplyr)
require(zetadiv)
require(sp)
require(rgdal)
require(geosphere)
require(stringr)
require(tidyr)

#wd <- "/home/cmb-07/sn1/alsimons/SCCWRP"
wd <- "~/Desktop/SCCWRP/Metagenomics/"
setwd(wd)

#Read in metagenomic count tables and format them as presence/absence tables.
communityInputRawPlate1 <- read.table("18SV9P1TableWithTaxonomy.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")
communityInputRawPlate2 <- read.table("18SV9P2TableWithTaxonomy.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")
#Remove spaces from column names for the count tables.
communityInputRawPlate1 <- rename(communityInputRawPlate1, OTUID = `OTU ID`)
communityInputRawPlate2 <- rename(communityInputRawPlate2, OTUID = `OTU ID`)

#Get a list of unique OTU names from count tables and convert to a data frame.
#uniqueOTUs <- as.data.frame(unique(c(communityInputRawPlate1$OTUID,communityInputRawPlate2$OTUID)))
uniqueOTUs <- rbind(communityInputRawPlate1[,c("OTUID","ConsensusLineage")],communityInputRawPlate2[,c("OTUID","ConsensusLineage")])
#colnames(uniqueOTUs) <- c("OTUID")
uniqueOTUs$OTUID <- as.character(uniqueOTUs$OTUID)
#Split OTU names into Domain through Genus+Species.
uniqueOTUs <- suppressWarnings(separate(uniqueOTUs,'ConsensusLineage',c("Domain", "Kingdom","Phylum","Class","Order","Family","GenusSpecies"),sep=";", extra="drop"))
uniqueOTUs$GenusSpecies <- trimws(uniqueOTUs$GenusSpecies,which="left") #Remove starting blank space from genus names
uniqueOTUs <- suppressWarnings(separate(uniqueOTUs,'GenusSpecies',c("Genus","Species"),sep=" ", extra="warn"))
uniqueOTUs <- uniqueOTUs[!duplicated(uniqueOTUs$OTUID),]
#Filter out ambiguous taxonomies
uniqueOTUs <- uniqueOTUs[uniqueOTUs$Domain!="Unassigned",]

#Create a merged metagenomic count table.
communityInput <- left_join(uniqueOTUs,communityInputRawPlate1,by="OTUID")
communityInput <- left_join(communityInput,communityInputRawPlate2,by="OTUID")

#Remove unnecessary sample columns.
communityInput <- communityInput[, -which(names(communityInput)  %in% c("DNAStandard","Ext-Blank1","Ext-Blank2","FB","NTC","SNAStandardII","ConsensusLineage.y","DNAstandardI","DNAstandardII","ConsensusLineage.x","202","230.x","230.y"))]

#Choose a taxonomic level to group count data by.
#Levels are Domain, Kingdom, Phylum, Class, Order, Family, GenusSpecies, OTUID
taxonomicLevels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species", "OTUID")
taxonomicLevel <- c("Species") #Choose a taxonomic level to aggregate count data on.
taxonomicIgnore <- taxonomicLevels[taxonomicLevels != taxonomicLevel]

#Remove unnecessary sample columns.
communityInput <- communityInput[, -which(names(communityInput)  %in% taxonomicIgnore)]
#Remove ambiguous or unassigned taxa for a selected taxonomic level.
communityInput <- communityInput[!is.na(communityInput[,colnames(communityInput)==taxonomicLevel]),]
communityInput <- communityInput[!grepl("environmental",communityInput$Species),]
#Aggregate a taxonomic level to aggregate count data on.
communityInputSummarized <- aggregate(.~Species,communityInput,sum)
#Convert abundance to presence/absence.
rownames(communityInputSummarized) <- Filter(is.character, communityInputSummarized)[,1]
communityInputSummarized[,which(colnames(communityInputSummarized)==taxonomicLevel)] <- NULL
communityInputSummarized[is.na(communityInputSummarized)] <- 0
#Keep only if at least three reads are present.
communityInputSummarized[communityInputSummarized < 2] <- 0
communityInputSummarized[communityInputSummarized > 2] <- 1

#Read in table linking sample IDs in the metagenomic table to sample station codes.
sampleIDs <- read.table("SampleStationCodesID.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Remove dubious samples
sampleIDs <- sampleIDs[which(sampleIDs$SampleNum!=202 & sampleIDs$SampleNum!=230),]

#Read in environmental metadata by station codes table.
metadata <- read.table("MetagenomicSampleSiteMetadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Merge in sample IDs to the metadata.
metadata <- merge(sampleIDs,metadata)

#Subset metadata by which samples are present in your community data set.
metadata <- suppressWarnings(metadata[which(metadata$SampleNum %in% as.numeric(colnames(communityInput))),])
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

#Assign a geographic cluster to sample sites in the metadata.
#Minimize the difference in group sizes during clustering.
mdist <- distm(metadata[,c("Longitude","Latitude")])
x <- 1
#Optimize the equality between cluster group sizes.
xlist <- c()
for(i in 1:1000){
  km <- kmeans(mdist,centers=4)
  xlist <- c(xlist,sd(table(km$cluster))/mean(table(km$cluster)))
  print(min(xlist))
  print(table(km$cluster))
}
while(x > min(xlist)){
  km <- kmeans(mdist,centers=4)
  x <- sd(table(km$cluster))/mean(table(km$cluster))
  print(x)
  print(table(km$cluster))
}
metadata$clust <- km$cluster

set.seed(1)
sample_Num <- 10
zetaMax <- 5
zetaAnalysis <- data.frame()
for(j in 1:100){
  for(i in unique(km$cluster)){
    tmp <- metadata[metadata$clust==i,]
    #Randomly subsample 10 samples from each cluster by land use grouping.
    metadataSubset <- tmp[sample(nrow(tmp),sample_Num),]
    #Reorder community data set columns so that their sample number order increases from left to right.
    #This will match the order of the metadata data frame where sample number order increases going down.
    data.spec <- communityInputSummarized[,as.character(metadataSubset$SampleNum)]
    #Create a species/site matrix for use in the zeta diversity functions.
    data.spec <- as.data.frame(t(data.spec))
    zetaDecay <- Zeta.decline.ex(data.spec,orders=1:zetaMax,rescale=TRUE,plot=FALSE)
    zeta_N <- Zeta.order.ex(data.spec,order=zetaMax,rescale=TRUE)$zeta.val #Higher order zeta diversity measure.
    zeta_Nsd <- Zeta.order.ex(data.spec,order=zetaMax,rescale=TRUE)$zeta.val.sd #Higher order zeta diversity measure standard deviation.
    zeta_1 <- Zeta.order.ex(data.spec,order=1,rescale=TRUE)$zeta.val #lower order zeta diversity measure.
    zeta_1sd <- Zeta.order.ex(data.spec,order=1,rescale=TRUE)$zeta.val.sd #lower order zeta diversity measure standard deviation.
    ExpExp <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    ExpAIC <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    PLExp <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    PLAIC <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    clusterID <- i
    meanLU <- mean(metadataSubset$LU)
    sdLU <- sd(metadataSubset$LU)
    meanAL <- mean(metadataSubset$site_elev)
    sdAL<- sd(metadataSubset$site_elev)
    numWS <- length(unique(metadataSubset$HUC8))
    meanDist <- mean(distm(metadataSubset[,c("Longitude","Latitude")]))
    sdDist <- sd(distm(metadataSubset[,c("Longitude","Latitude")]))
    meanP <- mean(metadataSubset$MaxP,na.rm=T)
    sdP <- sd(metadataSubset$MaxP,na.rm=T)
    meanOrthoP <- mean(metadataSubset$MaxOrthoP,na.rm=T)
    sdOrthoP <- sd(metadataSubset$MaxOrthoP,na.rm=T)
    meanN <- mean(metadataSubset$MaxN,na.rm=T)
    sdN <- sd(metadataSubset$MaxN,na.rm=T)
    print(paste(clusterID,meanLU,sdLU,meanAL,sdAL,meanDist,sdDist,numWS,meanP,sdP,meanOrthoP,sdOrthoP,meanN,sdN,zeta_N,zeta_Nsd,zeta_1,zeta_1sd,ExpExp,ExpAIC,PLExp,PLAIC))
    dataRow <- t(as.data.frame(list(c(clusterID,meanLU,sdLU,meanAL,sdAL,meanDist,sdDist,numWS,meanP,sdP,meanOrthoP,sdOrthoP,meanN,sdN,zeta_N,zeta_Nsd,zeta_1,zeta_1sd,ExpExp,ExpAIC,PLExp,PLAIC))))
    rownames(dataRow) <- NULL
    zetaAnalysis <- rbind(zetaAnalysis,dataRow)
  }
}
colnames(zetaAnalysis) <- c("clusterID","meanLU","sdLU","meanAL","sdAL","meanDist","sdDist","numWS","meanP","sdP","meanOrthoP","sdOrthoP","meanN","sdN","zeta_N","zeta_Nsd","zeta_1","zeta_1sd","ExpExp","ExpAIC","PLExp","PLAIC")
zetaAnalysis$zeta_Nscaled <- zetaAnalysis$zeta_N/zetaAnalysis$zeta_1
#Save zeta diversity analysis for a given taxonomic level.
write.table(zetaAnalysis,paste("zetaAnalysis18SV9",taxonomicLevel,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#
zetaAnalysis <- read.table("zetaAnalysis18SV9Species.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

require(ggplot2)
require(viridis)
zetaPlot <- ggplot(zetaAnalysis, aes(x=meanLU,y=zeta_Nscaled,color=meanN))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=zeta_Nscaled))
zetaPlot+xlab("Mean LU")+ylab("Zeta_5")+scale_color_gradientn("Mean N",colours = rev(plasma(10)))
zetaPlot <- ggplot(zetaAnalysis, aes(x=meanLU,y=zeta_Nscaled,color=numWS))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=zeta_Nscaled))
zetaPlot+xlab("Mean LU")+ylab("Zeta_5")+scale_color_gradientn("Num WS",colours = rev(plasma(10)))
zetaPlot <- ggplot(zetaAnalysis, aes(x=meanDist,y=zeta_Nscaled,color=meanAL))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=zeta_Nscaled))
zetaPlot+xlab("Mean Dist")+ylab("Zeta_5")+scale_color_gradientn("Mean AL",colours = rev(plasma(10)))


#Calculate how much zeta diversity of a particular order decays with distance
zetaDistance <- Zeta.ddecay(xy=metadata[,c("Latitude","Longitude")],data.spec=data.spec,order=5,distance.type="ortho",normalize=TRUE,plot=FALSE)
#How strongly correlated is zeta diversity with geographic distance?
cor.test(zetaDistance$zeta.val,zetaDistance$distance,method="spearman")
#

require(ggplot2)
require(ggmap)
require(maps)
require(mapview)
require(mapdata)
require(munsell)
require(leaflet)
require(devtools)
require(webshot)
require(viridis)
#Map data.
CalMap = leaflet(metadata) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=plasma(10),domain=metadata$clust)
CalMap %>% addCircleMarkers(color = ~ColorScale(clust), fill = TRUE,radius=0.1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~clust,title="SCCWRP sample<br>site clusters")


data.xy <- metadata[,c("Longitude","Latitude")]
data.env <- metadata[,c("LU","MaxN","MaxOrthoP","MaxP","elev_range","max_elev")]
zetaFactors <- Zeta.msgdm(data.spec=data.spec,data.env=data.env,xy=data.xy,order=4,sam=300,reg.type="glm",distance.type="ortho",normalize=FALSE,rescale=FALSE,control=list(maxit=100))
