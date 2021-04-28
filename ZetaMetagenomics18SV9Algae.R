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
require(picante)
require(relaimpo)

wd <- "/home/cmb-07/sn1/alsimons/SCCWRP"
wd <- "~/Desktop/SCCWRP/Metagenomics/"
setwd(wd)

#Read in metagenomic count tables and format them as presence/absence tables.
communityInputRawPlate1 <- read.table("18SV9P1TableWithTaxonomy.txt", header=T, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")
communityInputRawPlate2 <- read.table("18SV9P2TableWithTaxonomy.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=F, encoding = "UTF-8")
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

#Subset 18Sv9 OTU table to only contain taxonomic algal data from SCCWRP
#Generated here: https://github.com/levisimons/SCCWRP/blob/master/SCCWRPTaxonomyGenerator.R
uniqueAlgae <- read.table("AlgaeTaxonomies18SV9.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")

#Merge algal community data.
communityInput <- dplyr::left_join(uniqueAlgae,communityInputRawPlate1,by="OTUID")
communityInput <- dplyr::left_join(communityInput,communityInputRawPlate2,by="OTUID")
#Remove unnecessary sample columns.
communityInput <- communityInput[, -which(names(communityInput)  %in% c("DNAStandard","Ext-Blank1","Ext-Blank2","FB","NTC","SNAStandardII","ConsensusLineage.y","DNAstandardI","DNAstandardII","ConsensusLineage.x","202","230.x","230.y","FullTaxonomy"))]

#Choose a taxonomic level to group count data by.
#Levels are domain, kingdom, phylum, class, order, family, genus, species, OTUID
taxonomicLevels <- colnames(communityInput[,grep("^[A-Za-z]", colnames(communityInput))])
taxonomicLevel <- c("family") #Choose a taxonomic level to aggregate count data on.
taxonomicIgnore <- taxonomicLevels[taxonomicLevels != taxonomicLevel]

#Read in table linking sample IDs in the metagenomic table to sample station codes.
sampleIDs <- read.table("SampleStationCodesID.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Standarize date format.
sampleIDs$Date <- as.Date(sampleIDs$Date,format="%m/%d/%y")
sampleIDs$Date <- format(sampleIDs$Date,format="%m/%d/%y")
#Remove dubious samples
sampleIDs <- sampleIDs[which(sampleIDs$SampleNum!=202 & sampleIDs$SampleNum!=230),]
#Create unique identifier column.
sampleIDs$UniqueID <- paste(sampleIDs$StationCode,sampleIDs$Date)
#Only keep sample numbers and unique IDs.
sampleIDs <- sampleIDs[,c("SampleNum","UniqueID")]

#Read in environmental metadata by station codes table.
metadata <- read.table("MetagenomicSampleSiteMetadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Standardize date format and create a unique sample ID.
metadata$Date <- as.Date(metadata$Date,format="%m/%d/%y")
metadata$Date <- format(metadata$Date,format="%m/%d/%y")
metadata$UniqueID <- paste(metadata$StationCode,metadata$Date)

#Merge in sample IDs to the metadata.
metadata <- dplyr::left_join(sampleIDs,metadata,by=c("UniqueID"))

#Standardize elevation and land use data.
metadata$elev_range <- as.numeric(metadata$elev_range)
metadata$max_elev <- as.numeric(metadata$max_elev)
#Create aggregate upstream land use variable.
metadata$LU <- metadata$Ag_2011_5K+metadata$URBAN_2011_5K+metadata$CODE_21_2011_5K
#Extract sampling year
metadata$Year <- as.numeric(format(as.Date(metadata$Date,"%m/%d/%y"),"%Y"))

#Extract CSCI values.
CSCI <- read.table("ALS_DataRequest_CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Standardize date format and create a unique sample ID.
CSCI$sampledate <- as.Date(CSCI$sampledate,format="%Y-%m-%d")
CSCI$sampledate <- format(CSCI$sampledate,format="%m/%d/%y")
CSCI$UniqueID <- paste(CSCI$stationcode,CSCI$sampledate)
#Aggregate CSCI scores by sample over replicates.
CSCI <- CSCI[,c("UniqueID","csci")]
colnames(CSCI) <- c("UniqueID","csci")
CSCI <- as.data.frame(aggregate(csci ~ .,data=CSCI,mean))

#Merge in CSCI data.
metadata <- dplyr::left_join(metadata,CSCI,by="UniqueID")

#Read in algal streams condition index data.
ASCI <- read.table("ASCI.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
#Standardize the date format and create unique identifier column.
ASCI$sampledate <- as.Date(ASCI$sampledate,format="%Y-%m-%d")
ASCI$sampledate <- format(ASCI$sampledate,format="%m/%d/%y")
ASCI$UniqueID <- paste(ASCI$stationcode,ASCI$sampledate)
#Only keep ASCI results for all algae.
ASCI <- ASCI[ASCI$assemblage=="Hybrid" & ASCI$metric=="ASCI",]

#Aggregate CSCI scores by sample over replicates.
ASCI <- ASCI[,c("UniqueID","result")]
colnames(ASCI) <- c("UniqueID","asci")
ASCI <- as.data.frame(aggregate(asci ~ .,data=ASCI,mean))

#Merge in ASCI data.
metadata <- dplyr::left_join(metadata,ASCI,by="UniqueID")

#Add land use bands.  Band1 (Low): 0-3%, Band2 (Intermediate): 3-15%, Band3 (High): 15-100%.
metadata$LUBand <- as.numeric(case_when(metadata$LU <= 3 ~ "1", metadata$LU > 3 & metadata$LU <= 15 ~ "2", metadata$LU > 15 ~ "3", TRUE ~ as.character(metadata$LU)))
#Remove samples with missing metadata.
metadata <- metadata[!is.na(metadata$LUBand) & !is.na(metadata$csci) & !is.na(metadata$asci),]
#Remove duplicate metadata rows.
metadata$SampleNum <- NULL
metadata <- metadata[!duplicated(metadata),]

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
#Move taxonomic names onto data frame row names.
rownames(communityInputSummarized) <- Filter(is.character, communityInputSummarized)[,1]
communityInputSummarized[,which(colnames(communityInputSummarized)==taxonomicLevel)] <- NULL
#Determine the average number of reads per sample for a given unique site/location.
communityInputSummarized <- as.data.frame(t(communityInputSummarized))
communityInputSummarized$SampleNum <- as.numeric(rownames(communityInputSummarized))
communityInputSummarized <- dplyr::left_join(communityInputSummarized,sampleIDs,by=c("SampleNum"))
communityInputSummarized$SampleNum <- NULL
communityInputSummarized <- as.data.frame(aggregate(formula(paste0(". ~ UniqueID")),communityInputSummarized,mean,na.action = na.omit))
rownames(communityInputSummarized) <- communityInputSummarized$UniqueID
communityInputSummarized$UniqueID <- NULL

#Convert abundance to presence/absence.
communityInputSummarized[is.na(communityInputSummarized)] <- 0
#Keep only if at least two reads are present.
communityInputSummarized[communityInputSummarized <= 2] <- 0
communityInputSummarized[communityInputSummarized > 2] <- 1

#Remove community entries with missing metadata.
communityInputSummarized <- communityInputSummarized[rownames(communityInputSummarized) %in% metadata$UniqueID,]

set.seed(1)
sample_Num <- 15
zetaMax <- 10
zetaAnalysis <- data.frame()
for(j in 1:100){
  for(i in unique(metadata$LUBand)){
    #Subset by land use band.
    tmp <- metadata[metadata$LUBand==i,]
    #Randomly subsample sample_Num samples from each cluster by land use grouping.
    metadataSubset <- tmp[sample(nrow(tmp),sample_Num),]
    #Subset the sample by taxa presence/absence data frame to only contain the same samples represented by metadata sites.
    #This is the input for zeta diversity calculations.
    data.spec <- communityInputSummarized[rownames(communityInputSummarized) %in% metadataSubset$UniqueID,]
    #Calculate zeta diversity decay statistics.
    zetaDecay <- Zeta.decline.ex(data.spec,orders=1:zetaMax,rescale=TRUE,plot=FALSE)
    zeta_1 <- Zeta.order.ex(data.spec,order=1,rescale=TRUE)$zeta.val #lower order zeta diversity measure.
    zeta_1sd <- Zeta.order.ex(data.spec,order=1,rescale=TRUE)$zeta.val.sd #lower order zeta diversity measure standard deviation.
    zeta_2 <- zetaDecay$zeta.val[2] #lower order zeta diversity measure.
    zeta_2sd <- zetaDecay$zeta.val.sd[2] #lower order zeta diversity measure standard deviation.
    zeta_N <- zetaDecay$zeta.val[zetaMax] #Higher order zeta diversity measure.
    zeta_Nsd <- zetaDecay$zeta.val.sd[zetaMax] #Higher order zeta diversity measure standard deviation.
    ExpExp <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    ExpAIC <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    PLExp <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    PLAIC <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    LUBand <- i
    #Calculate mean environmental variables per sample group.
    meanLU <- mean(metadataSubset$LU)
    sdLU <- sd(metadataSubset$LU)
    meanAL <- mean(metadataSubset$site_elev)
    sdAL<- sd(metadataSubset$site_elev)
    meanDist <- mean(distm(metadataSubset[,c("new_long","new_lat")]))
    sdDist <- sd(distm(metadataSubset[,c("new_long","new_lat")]))
    #Get the mean and standard deviation on the ASCI scores per sample group
    meanASCI <- mean(metadataSubset$asci,na.rm=T)
    sdASCI <- sd(metadataSubset$asci,na.rm=T)
    #Aggregate zeta analysis output
    print(paste(LUBand,meanLU,sdLU,meanAL,sdAL,meanDist,sdDist,meanASCI,sdASCI,zeta_N,zeta_Nsd,zeta_1,zeta_1sd,zeta_2,zeta_2sd,ExpExp,ExpAIC,PLExp,PLAIC))
    dataRow <- t(as.data.frame(list(c(LUBand,meanLU,sdLU,meanAL,sdAL,meanDist,sdDist,meanASCI,sdASCI,zeta_N,zeta_Nsd,zeta_1,zeta_1sd,zeta_2,zeta_2sd,ExpExp,ExpAIC,PLExp,PLAIC))))
    rownames(dataRow) <- NULL
    zetaAnalysis <- rbind(zetaAnalysis,dataRow)
  }
}
colnames(zetaAnalysis) <- c("LUBand","meanLU","sdLU","meanAL","sdAL","meanDist","sdDist","meanASCI","sdASCI","zeta_N","zeta_Nsd","zeta_1","zeta_1sd","zeta_2","zeta_2sd","ExpExp","ExpAIC","PLExp","PLAIC")
indx <- sapply(zetaAnalysis, is.factor)
zetaAnalysis[indx] <- lapply(zetaAnalysis[indx], function(x) as.numeric(as.character(x)))
zetaAnalysis <- do.call(data.frame,lapply(zetaAnalysis, function(x) replace(x, is.infinite(x),NA)))
#Save zeta diversity analysis for a given taxonomic level.
write.table(zetaAnalysis,paste("zetaAnalysis18SV9Algae",taxonomicLevel,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#

##To run locally.
zetaAnalysis <- read.table(paste("zetaAnalysis18SV9Algae",taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

zetaModel <- lm(meanASCI~zeta_1+zeta_2+zeta_N,data=zetaAnalysis)
plot(zetaModel$model$meanASCI,zetaModel$fitted.values)
cor.test(zetaModel$model$meanASCI,zetaModel$fitted.values)
zetaAnalysis$modeledASCI <- zetaModel$fitted.values
calc.relimp(zetaModel)
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(zetaModel)

# Assessing R2 shrinkage using 10-Fold Cross-Validation 
require(bootstrap)
require(caret)
set.seed(1)
train.control <- trainControl(method="repeatedcv",number=10,repeats=10)
ASCImodel <- train(meanASCI~zeta_1+zeta_2+zeta_N,data=zetaAnalysis,method="lm",trControl=train.control)
print(ASCImodel)

#Comparing mean and modeled CSCI versus environmental parameters
zetaVarModel1 <- lm(modeledASCI~meanAL+meanLU+meanDist,data=zetaAnalysis)
summary(zetaVarModel1)
calc.relimp(zetaVarModel1)
anova(zetaVarModel1)
zetaVarModel2 <- lm(meanASCI~meanAL+meanLU+meanDist,data=zetaAnalysis)
summary(zetaVarModel2)
calc.relimp(zetaVarModel2)
anova(zetaVarModel2)


require(ggplot2)
require(viridis)
zetaPlot <- ggplot(zetaAnalysis, aes(x=meanCSCI,y=modeledCSCI,color=meanLU))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=modeledCSCI))
zetaPlot+xlab("Mean CSCI")+ylab("Modeled CSCI")+scale_color_gradientn("Mean LU",colours = rev(plasma(10)))
zetaPlot <- ggplot(zetaAnalysis, aes(x=meanCSCI,y=modeledCSCI,color=meanAL))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=modeledCSCI))
zetaPlot+xlab("Mean CSCI")+ylab("Modeled CSCI")+scale_color_gradientn("Mean AL",colours = rev(plasma(10)))
zetaPlot <- ggplot(zetaAnalysis, aes(x=meanLU,y=zeta_N,color=meanLU))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=zeta_N))
zetaPlot+xlab("Mean LU")+ylab("Zeta_5")+scale_color_gradientn("Mean LU",colours = rev(plasma(10)))
zetaPlot <- ggplot(zetaAnalysis, aes(x=meanLU,y=PLExp,color=meanN))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=PLExp))
zetaPlot+xlab("Mean LU")+ylab("Zeta_1")+scale_color_gradientn("Mean N",colours = rev(plasma(10)))

#Check for correlation patterns between zeta diversity and environmental parameters.
require("PerformanceAnalytics")
chart.Correlation(zetaAnalysis[,c("modeledCSCI","meanCSCI","meanLU","meanAL","meanDist","zeta_1","zeta_2","zeta_3","zeta_4","zeta_N")], histogram=TRUE, method="pearson")

#Correlation plots of mean and modeled CSCI scores, zeta diversity measures, and environmental parameters.
require(Hmisc)
require(corrplot)
communityType <- "Algae"
taxonomicLevels <- c("species","genus","family","order")
for(taxonomicLevel in taxonomicLevels){
  zetaAnalysis <- read.table(paste("zetaAnalysis18SV9",communityType,taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
  zetaModel <- lm(meanASCI~zeta_1+zeta_2+zeta_N,data=zetaAnalysis)
  zetaAnalysis$modeledASCI <- zetaModel$fitted.values
  zetaCor <- zetaAnalysis[,c("meanLU","meanAL","meanDist","zeta_1","zeta_2","zeta_N","meanASCI","modeledASCI")]
  zetaCor <- rcorr(as.matrix(zetaCor),type="pearson")
  corr <- zetaCor$r
  p.mat <- zetaCor$P
  colnames(corr) <- c("Land Use","Altitude","Distance",":zeta[1]",":zeta[2]",":zeta[10]","Mean H_ASCI","Modeled H_ASCI")
  rownames(corr) <- c("Land Use","Altitude","Distance",":zeta[1]",":zeta[2]",":zeta[10]","Mean H_ASCI","Modeled H_ASCI")
  png(paste("zetaAnalysis18SV9Algae",taxonomicLevel,".png",sep=""),width=7,height=7,units="in",res=600)
  corrplot(corr = corr, p.mat = p.mat, diag = FALSE, type="lower", sig.level = 0.0001, tl.col="black", tl.srt=45, tl.cex=1.3, order="original",mar=c(0,0,3,0))
  mtext(paste("18S V9 soft-bodied algae aggregated to",taxonomicLevel), at=2.5, line=3, cex=1.3)
  dev.off()
}

#Calculate how much zeta diversity of a particular order decays with distance
data.spec <- communityInputSummarized[,as.character(metadata$SampleNum)]
#Create a species/site matrix for use in the zeta diversity functions.
data.spec <- as.data.frame(t(data.spec))
zetaDistance <- Zeta.ddecay(xy=metadata[,c("Latitude","Longitude")],data.spec=data.spec,order=5,distance.type="ortho",normalize="Jaccard",plot=FALSE)
#How strongly correlated is zeta diversity with geographic distance?
cor.test(zetaDistance$zeta.val,zetaDistance$distance,method="spearman")
zetaDecay <- Zeta.decline.ex(data.spec,orders=1:zetaMax,rescale=TRUE,plot=TRUE)

#Compare community assembly profiles.
require(IDPmisc)
mean(NaRV.omit(zetaAnalysis$ExpAIC))
mean(NaRV.omit(zetaAnalysis$PLAIC))
cor.test(zetaAnalysis$meanLU,zetaAnalysis$PLAIC-zetaAnalysis$ExpAIC)
t.test(NaRV.omit(zetaAnalysis$PLAIC),NaRV.omit(zetaAnalysis$ExpAIC),alternative="two.sided")

#Zeta.varpart returns a data frame with one column containing the variation explained by each component 
#a (the variation explained by distance alone),
#b (the variation explained by either distance or the environment),
#c (the variation explained by the environment alone) and 
#d (the unexplained variation).
data.env <- metadata[,c("LU","site_elev")]
data.env[] <- lapply(data.env,as.numeric)
rownames(data.env) <- metadata[,c("SampleNum")]
#data.env <- na.omit(data.env)
data.xy <- metadata[,c("Longitude","Latitude")]
rownames(data.xy) <- metadata[,c("SampleNum")]
data.xy <- data.xy[which(rownames(data.xy) %in% rownames(data.env)),]
data.spec <- as.data.frame(t(communityInputSummarized))
tmp <- data.spec[which(rownames(data.spec) %in% rownames(data.xy)),]
zetaFactors <- Zeta.msgdm(data.spec=tmp,data.env=data.env,xy=data.xy,order=2,reg.type="glm",distance.type="ortho",normalize=FALSE,empty.row=0,rescale=FALSE,control=list(maxit=1000))
zetaVar <- Zeta.varpart(zetaFactors)
zetaVar

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
#metadata$Ecoregion <- as.numeric(metadata$LU)
CalMap = leaflet(metadata) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=plasma(10),domain=metadata$LU)
CalMap %>% addCircleMarkers(color = ~ColorScale(LU), fill = TRUE,radius=5,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~Ecoregion,title="SCCWRP sample<br>% Land Use")
