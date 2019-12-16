rm(list=ls())
require("plyr")
require(dplyr)
require(zetadiv)
require(sp)
require(rgdal)
require(geosphere)
require(stringr)
require(tidyr)
require(naniar)
require(taxize)

wd <- "/home/cmb-07/sn1/alsimons/SCCWRP"
#wd <- "~/Desktop/SCCWRP/Metagenomics/"
setwd(wd)

#Choose a community type.
#"" for all data.
#"algae" for algal communities.
#"BMIs" for BMI communities
communityType <- "BMIs"

#Read in metagenomic count tables and format them as presence/absence tables.
communityInputRawPlate1 <- read.table("18SV9P1TableWithTaxonomy.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")
communityInputRawPlate2 <- read.table("18SV9P2TableWithTaxonomy.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")
#Remove spaces from column names for the count tables.
communityInputRawPlate1 <- dplyr::rename(communityInputRawPlate1, OTUID = `OTU ID`)
communityInputRawPlate2 <- dplyr::rename(communityInputRawPlate2, OTUID = `OTU ID`)

#Get a list of unique OTU names from count tables and convert to a data frame.
#uniqueOTUs <- as.data.frame(unique(c(communityInputRawPlate1$OTUID,communityInputRawPlate2$OTUID)))
uniqueOTUs <- rbind(communityInputRawPlate1[,c("OTUID","ConsensusLineage")],communityInputRawPlate2[,c("OTUID","ConsensusLineage")])
#colnames(uniqueOTUs) <- c("OTUID")
uniqueOTUs$OTUID <- as.character(uniqueOTUs$OTUID)
#Remove superfluous strings from taxonomic labels.
uniqueOTUs$ConsensusLineage <- gsub("D_[0-9]+__","",uniqueOTUs$ConsensusLineage)
uniqueOTUs$ConsensusLineage <- gsub("g__","",uniqueOTUs$ConsensusLineage)
#Split OTU names into Domain through Genus+Species.
uniqueOTUs$FullTaxonomy <- uniqueOTUs$ConsensusLineage
uniqueOTUs <- suppressWarnings(separate(uniqueOTUs,'ConsensusLineage',c("Rank1", "Rank2","Rank3","Rank4","Rank5","Rank6","Rank7and8"),sep=";", extra="drop"))
uniqueOTUs$Rank7and8 <- trimws(uniqueOTUs$Rank7and8,which="left") #Remove starting blank space from genus names
uniqueOTUs <- suppressWarnings(separate(uniqueOTUs,'Rank7and8',c("Rank7","Rank8"),sep=" ", extra="warn"))
uniqueOTUs <- uniqueOTUs[!duplicated(uniqueOTUs$OTUID),]
#Filter out ambiguous taxonomies
uniqueOTUs <- uniqueOTUs[uniqueOTUs$Rank1!="Unassigned" & uniqueOTUs$Rank1!="Ambiguous_taxa",]
ambiguousList <- c("Incertae Sedis","metagenome","sp.","environmental","eukaryote","uncultured","soil","Ambiguous_taxa","group","cf.","aff.","gen.","marine","cf","unidentified","Uncultured","invertebrate")
ambiguousList <- as.list(ambiguousList)
uniqueOTUs <- data.frame(lapply(uniqueOTUs, trimws), stringsAsFactors = FALSE)
uniqueOTUs <- replace_with_na_all(data=uniqueOTUs,condition=~.x %in% as.list(ambiguousList))

#Subset 18Sv9 OTU table to only contain taxonomic BMI data from SCCWRP
#Generated here: https://github.com/levisimons/SCCWRP/blob/master/SCCWRPTaxonomyGenerator.R
uniqueBMIs <- read.table("BMITaxonomies.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")

#Subset 18Sv9 OTU table to only contain taxonomic algal data from SCCWRP
#Generated here: https://github.com/levisimons/SCCWRP/blob/master/SCCWRPTaxonomyGenerator.R
uniqueAlgae <- read.table("AlgaeTaxonomies.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")

#Create a merged metagenomic count table.
if(communityType==""){
  communityInput <- left_join(uniqueOTUs,communityInputRawPlate1,by="OTUID")
}
if(communityType=="algae"){
  communityInput <- left_join(uniqueAlgae,communityInputRawPlate1,by="OTUID")
}
if(communityType=="BMIs"){
  communityInput <- left_join(uniqueBMIs,communityInputRawPlate1,by="OTUID")
}
communityInput <- left_join(communityInput,communityInputRawPlate2,by="OTUID")
#Remove unnecessary sample columns.
communityInput <- communityInput[, -which(names(communityInput)  %in% c("DNAStandard","Ext-Blank1","Ext-Blank2","FB","NTC","SNAStandardII","ConsensusLineage.y","DNAstandardI","DNAstandardII","ConsensusLineage.x","202","230.x","230.y","FullTaxonomy"))]

#Choose a taxonomic level to group count data by.
#Levels are Domain, Kingdom, Phylum, Class, Order, Family, GenusSpecies, OTUID
taxonomicLevels <- c("kingdom", "phylum", "class", "order", "family", "genus","species", "OTUID")
taxonomicLevel <- c("OTUID") #Choose a taxonomic level to aggregate count data on.
taxonomicIgnore <- taxonomicLevels[taxonomicLevels != taxonomicLevel]

#Remove unnecessary sample columns.
communityInput <- communityInput[, -which(names(communityInput)  %in% taxonomicIgnore)]
communityInput <- communityInput[!is.na(communityInput[,colnames(communityInput)==taxonomicLevel]),]
#Aggregate a taxonomic level to aggregate count data on.
communityInput[is.na(communityInput)] <- 0
if(taxonomicLevel=="OTUID"){
  communityInputSummarized <- as.data.frame(communityInput)
}
if(taxonomicLevel!="OTUID"){
  communityInputSummarized <- as.data.frame(aggregate(formula(paste0(". ~ ",taxonomicLevel)),communityInput,sum,na.action = na.omit))
}
#Convert abundance to presence/absence.
rownames(communityInputSummarized) <- Filter(is.character, communityInputSummarized)[,1]
communityInputSummarized[,which(colnames(communityInputSummarized)==taxonomicLevel)] <- NULL
communityInputSummarized[is.na(communityInputSummarized)] <- 0
#Keep only if at least three reads are present.
communityInputSummarized[communityInputSummarized <= 2] <- 0
communityInputSummarized[communityInputSummarized > 2] <- 1

#Calculte taxonomic richness by sample.
communityRichness <- as.data.frame(colSums(communityInputSummarized))
communityRichness$SampleNum <- as.numeric(rownames(communityRichness))
colnames(communityRichness) <- c("Richness","SampleNum")

#Read in table linking sample IDs in the metagenomic table to sample station codes.
sampleIDs <- read.table("SampleStationCodesID.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Remove dubious samples
sampleIDs <- sampleIDs[which(sampleIDs$SampleNum!=202 & sampleIDs$SampleNum!=230),]

#Read in environmental metadata by station codes table.
metadata <- read.table("MetagenomicSampleSiteMetadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Merge in sample IDs to the metadata.
metadata <- merge(sampleIDs,metadata,by=c("StationCode","Date"))

#Subset metadata by which samples are present in your community data set.
metadata <- suppressWarnings(metadata[which(metadata$SampleNum %in% as.numeric(colnames(communityInput))),])
#Sort metadata by sample number.
metadata <- arrange(metadata,SampleNum)
metadata$elev_range <- as.numeric(metadata$elev_range)
metadata$max_elev <- as.numeric(metadata$max_elev)
#Create aggregate upstream land use variable.
metadata$LU <- metadata$Ag_2011_5K+metadata$URBAN_2011_5K+metadata$CODE_21_2011_5K
#Merge in community richness.
metadata <- left_join(metadata,communityRichness,by="SampleNum")

#Read in nitrogen and phosporus site data.
NPdata <- read.table("NandP_labdata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
NPdata[NPdata=="ND" | NPdata==""] <- NA
NPdata[,2:4] <- sapply(NPdata[,2:4],as.numeric)
NPdata[NPdata<0] <- NA

#Merge in nitrogen and phosphorus site data into metada.
metadata <- left_join(metadata,NPdata,by=c("StationCode"))

#Read in watershed by sample site location data.
SCCWRP <- read.table("MetagenomicSitesWithWatershedsEcoregions.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
SCCWRP <- SCCWRP[,c("Latitude","Longitude","HUC8","PSA6","PSA8")]
#Remove duplicate rows
SCCWRP <- SCCWRP[!duplicated(SCCWRP),]
#Create HUC 6 watershed column.
SCCWRP$HUC6 <- substr(SCCWRP$HUC8,1,nchar(SCCWRP$HUC8)-2)
#Create HUC 4 watershed column.
SCCWRP$HUC4 <- substr(SCCWRP$HUC8,1,nchar(SCCWRP$HUC8)-4)
#Create PSA4 level ecoregion
SCCWRP$PSA4 <- SCCWRP$PSA6
SCCWRP$PSA4 <- gsub("Deserts Modoc","HighInterior",SCCWRP$PSA4)
SCCWRP$PSA4 <- gsub("Sierra Nevada","HighInterior",SCCWRP$PSA4)
SCCWRP$PSA4 <- gsub("Central Valley","CentralCalifornia",SCCWRP$PSA4)
SCCWRP$PSA4 <- gsub("Chaparral","CentralCalifornia",SCCWRP$PSA4)
SCCWRP$Ecoregion <- as.factor(SCCWRP$PSA4)
levels(SCCWRP$Ecoregion) <- 1:length(unique(SCCWRP$Ecoregion))

#Merge in watershed data into metadata
metadata <- left_join(metadata,SCCWRP,by=c("Latitude","Longitude"))
metadata <- metadata[!is.na(metadata$Latitude),]
metadata <- metadata[!is.na(metadata$Longitude),]

set.seed(1)
sample_Num <- 10
zetaMax <- 5
zetaAnalysis <- data.frame()
for(j in 1:100){
  for(i in unique(metadata$Ecoregion)){
    #Subset by ecoregion.
    tmp <- metadata[metadata$Ecoregion==i,]
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
    Ecoregion <- i
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
    print(j)
    print(paste(Ecoregion,meanLU,sdLU,meanAL,sdAL,meanDist,sdDist,numWS,meanP,sdP,meanOrthoP,sdOrthoP,meanN,sdN,zeta_N,zeta_Nsd,zeta_1,zeta_1sd,ExpExp,ExpAIC,PLExp,PLAIC))
    dataRow <- t(as.data.frame(list(c(Ecoregion,meanLU,sdLU,meanAL,sdAL,meanDist,sdDist,numWS,meanP,sdP,meanOrthoP,sdOrthoP,meanN,sdN,zeta_N,zeta_Nsd,zeta_1,zeta_1sd,ExpExp,ExpAIC,PLExp,PLAIC))))
    rownames(dataRow) <- NULL
    zetaAnalysis <- rbind(zetaAnalysis,dataRow)
  }
}
colnames(zetaAnalysis) <- c("Ecoregion","meanLU","sdLU","meanAL","sdAL","meanDist","sdDist","numWS","meanP","sdP","meanOrthoP","sdOrthoP","meanN","sdN","zeta_N","zeta_Nsd","zeta_1","zeta_1sd","ExpExp","ExpAIC","PLExp","PLAIC")
indx <- sapply(zetaAnalysis, is.factor)
zetaAnalysis[indx] <- lapply(zetaAnalysis[indx], function(x) as.numeric(as.character(x)))
zetaAnalysis$zeta_Nscaled <- zetaAnalysis$zeta_N/zetaAnalysis$zeta_1
#Save zeta diversity analysis for a given taxonomic level.
write.table(zetaAnalysis,paste("zetaAnalysis18SV9",communityType,taxonomicLevel,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#

##To run locally.
zetaAnalysis <- read.table(paste("zetaAnalysis18SV9",communityType,taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

require(ggplot2)
require(viridis)
zetaPlot <- ggplot(zetaAnalysis, aes(x=meanLU,y=zeta_Nscaled,color=meanN))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=zeta_Nscaled))
zetaPlot+xlab("Mean LU")+ylab("Zeta_5")+scale_color_gradientn("Mean N",colours = rev(plasma(10)))
zetaPlot <- ggplot(zetaAnalysis, aes(x=meanLU,y=zeta_Nscaled,color=meanAL))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=zeta_Nscaled))
zetaPlot+xlab("Mean LU")+ylab("Zeta_5")+scale_color_gradientn("Mean Al",colours = rev(plasma(10)))
zetaPlot <- ggplot(zetaAnalysis, aes(x=meanLU,y=PLExp,color=meanN))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=zeta_Nscaled))
zetaPlot+xlab("Mean LU")+ylab("Zeta_1")+scale_color_gradientn("Mean N",colours = rev(plasma(10)))

#Check for correlation patterns between zeta diversity and environmental parameters.
require("PerformanceAnalytics")
chart.Correlation(zetaAnalysis[,c("zeta_Nscaled","meanLU","meanAL","meanN","meanP","meanOrthoP")], histogram=TRUE, method="spearman")

#Calculate how much zeta diversity of a particular order decays with distance
data.spec <- communityInputSummarized[,as.character(metadata$SampleNum)]
#Create a species/site matrix for use in the zeta diversity functions.
data.spec <- as.data.frame(t(data.spec))
zetaDistance <- Zeta.ddecay(xy=metadata[,c("Latitude","Longitude")],data.spec=data.spec,order=5,distance.type="ortho",normalize="Jaccard",plot=FALSE)
#How strongly correlated is zeta diversity with geographic distance?
cor.test(zetaDistance$zeta.val,zetaDistance$distance,method="spearman")

#Compare community assembly profiles.
require(IDPmisc)
mean(NaRV.omit(zetaAnalysis$PLAIC))
mean(NaRV.omit(zetaAnalysis$ExpAIC))
t.test(NaRV.omit(zetaAnalysis$PLAIC),NaRV.omit(zetaAnalysis$ExpAIC),alternative="two.sided")

#Zeta.varpart returns a data frame with one column containing the variation explained by each component 
#a (the variation explained by distance alone),
#b (the variation explained by either distance or the environment),
#c (the variation explained by the environment alone) and 
#d (the unexplained variation).
data.env <- metadata[,c("LU","MaxN","MaxOrthoP","MaxP","elev_range","max_elev")]
rownames(data.env) <- metadata[,c("SampleNum")]
data.env <- na.omit(data.env)
data.xy <- metadata[,c("Longitude","Latitude")]
rownames(data.xy) <- metadata[,c("SampleNum")]
data.xy <- data.xy[which(rownames(data.xy) %in% rownames(data.env)),]
tmp <- data.spec[which(rownames(data.spec) %in% rownames(data.xy)),]
zetaFactors <- Zeta.msgdm(data.spec=tmp,data.env=data.env,xy=data.xy,order=5,reg.type="glm",distance.type="ortho",normalize=FALSE,rescale=FALSE,control=list(maxit=100))
zetaVar <- Zeta.varpart(zetaFactors)

#Mapping metadata
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
metadata$Ecoregion <- as.numeric(metadata$Ecoregion)
CalMap = leaflet(metadata) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=plasma(10),domain=metadata$Ecoregion)
CalMap %>% addCircleMarkers(color = ~ColorScale(Ecoregion), fill = TRUE,radius=0.1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~Ecoregion,title="SCCWRP sample<br>site ecoregions")
