rm(list=ls())
require(dplyr)
require(zetadiv)
require(geosphere)
require(relaimpo)

wd <- "/project/noujdine_61/alsimons/SCCWRP"
wd <- "~/Desktop/SCCWRP/MorphologicalBMIsAlgae/"
setwd(wd)

#Read in morphological stream community data describing all algae (soft algae and diatoms)
#and format them as presence/absence tables.
AlgalInput <- read.table("AlgTaxa.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
#Standarize date format.
AlgalInput$sampledate <- as.Date(AlgalInput$sampledate,format="%Y-%m-%d")
AlgalInput$sampledate <- format(AlgalInput$sampledate,format="%m/%d/%y")
#Create unique sample identifier.
AlgalInput$UniqueID <- paste(AlgalInput$stationcode,AlgalInput$sampledate)
#Subset columns of interest
AlgalInput <- AlgalInput[,c("stationcode","sampledate","replicate","UniqueID","finalid")]

#Add in additional algal morphology data.
AlgalInput2 <- read.table("AlgTaxa2.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
#Standarize date format.
AlgalInput2$sampledate <- as.Date(AlgalInput2$sampledate,format="%m/%d/%y")
AlgalInput2$sampledate <- format(AlgalInput2$sampledate,format="%m/%d/%y")
#Create unique sample identifier.
AlgalInput2$UniqueID <- paste(AlgalInput2$stationcode,AlgalInput2$sampledate)
#Subset columns of interest
AlgalInput2 <- AlgalInput2[,c("stationcode","sampledate","replicate","UniqueID","finalid")]

#Create merged algal morphology data set.
AlgalInput <- rbind(AlgalInput,AlgalInput2)

#Generated here: https://github.com/nuzhdinlab/SCCWRP/blob/master/MorphologicalTaxonomyGenerator.R
uniqueAlgae <- read.table("AlgalTaxonomiesMorphological.txt", header=T, sep="\t",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")

#Create merged data frame with sample data and full taxonomic information.
communityInput <- dplyr::left_join(uniqueAlgae,AlgalInput[,c("stationcode","sampledate","finalid")],by=c("LeafTaxa"="finalid"))

#Choose a taxonomic level to group count data by.
#Levels are domain, kingdom, phylum, class, order, family, genus, species, OTUID
taxonomicLevels <- colnames(uniqueAlgae[,grep("^[A-Za-z]", colnames(uniqueAlgae))])
taxonomicLevel <- c("order") #Choose a taxonomic level to aggregate count data on.
taxonomicIgnore <- taxonomicLevels[taxonomicLevels != taxonomicLevel]

#Remove unnecessary sample columns.
communityInput$UniqueID <- paste(communityInput$stationcode,communityInput$sampledate)
communityInput <- communityInput[,c("UniqueID",taxonomicLevel)]
#Remove empty taxonomic rows.
communityInput <- communityInput[!is.na(communityInput[,taxonomicLevel]),]

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
ASCI <- ASCI[ASCI$assemblage=="Diatom" & ASCI$metric=="ASCI",]

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
metadata$SampleNum <- NULL
metadata <- metadata[!duplicated(metadata),]
rownames(metadata) <- metadata$UniqueID

#Only keep community data for samples which have both morphological and metagenomic data.
communityInput <- communityInput[communityInput$UniqueID %in% metadata$UniqueID,]

#Collapse duplicate rows so each row is a unique pairing of sample by taxa.
communityInput <- communityInput[!duplicated(communityInput),]

#Convert community data to a presence/absence data frame.
#Rows for sample IDs and columns for taxa aggregated to a particular level.
communityInputSummarized <- data.frame(matrix(nrow=length(unique(communityInput[,"UniqueID"])),ncol=length(unique(communityInput[,taxonomicLevel]))))
colnames(communityInputSummarized) <- unique(communityInput[,taxonomicLevel])
rownames(communityInputSummarized) <- unique(communityInput[,"UniqueID"])
for(sample in unique(communityInput[,"UniqueID"])){
  communityInputSummarized[sample,] <- as.numeric(unique(communityInput[,taxonomicLevel]) %in% communityInput[communityInput$UniqueID==sample,taxonomicLevel])
}
#Ensure the rows in the site/species data frame include all the sites found with metadata.
un1 <- union(rownames(communityInputSummarized),rownames(metadata))
tmp1 <- as.data.frame(matrix(0,ncol=ncol(communityInputSummarized),nrow=length(un1),dimnames=list(un1,names(communityInputSummarized))))
tmp2 <- tmp1
tmp1[rownames(tmp1) %in% rownames(communityInputSummarized),] <- communityInputSummarized
communityInputSummarized <- tmp1

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
write.table(zetaAnalysis,paste("zetaAnalysisMorphologicalAlgae",taxonomicLevel,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#

##To run locally.
taxonomicLevel <- "order"
zetaAnalysis <- read.table(paste("zetaAnalysisMorphologicalAlgae",taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

zetaModel <- lm(meanASCI~zeta_1+zeta_2+zeta_N,data=zetaAnalysis)
plot(zetaModel$model$meanASCI,zetaModel$fitted.values)
cor.test(zetaModel$model$meanASCI,zetaModel$fitted.values)
zetaAnalysis$modeledASCI <- zetaModel$fitted.values
calc.relimp(zetaModel)
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(zetaModel)
dev.off()

# Assessing R2 shrinkage using 10-Fold Cross-Validation 
#require(bootstrap)
require(caret)
set.seed(1)
train.control <- trainControl(method="repeatedcv",number=10,repeats=10)
ASCImodel <- train(meanASCI~zeta_1+zeta_2+zeta_N,data=zetaAnalysis,method="lm",trControl=train.control)
print(cor.test(zetaModel$model$meanASCI,zetaModel$fitted.values))
print(ASCImodel)

#Comparing mean and modeled ASCI versus environmental parameters
zetaVarModel1 <- lm(modeledASCI~meanAL+meanLU+meanDist,data=zetaAnalysis)
summary(zetaVarModel1)
calc.relimp(zetaVarModel1)
anova(zetaVarModel1)
zetaVarModel2 <- lm(meanASCI~meanAL+meanLU+meanDist,data=zetaAnalysis)
summary(zetaVarModel2)
calc.relimp(zetaVarModel2)
anova(zetaVarModel2)

#Correlation plots of mean and modeled ASCI scores, zeta diversity measures, and environmental parameters.
require(Hmisc)
require(corrplot)
communityType <- "Algae"
taxonomicLevels <- c("species","genus","family","order")
for(taxonomicLevel in taxonomicLevels){
  zetaAnalysis <- read.table(paste("zetaAnalysisMorphologicalAlgae",taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
  zetaModel <- lm(meanASCI~zeta_1+zeta_2+zeta_N,data=zetaAnalysis)
  zetaAnalysis$modeledASCI <- zetaModel$fitted.values
  zetaCor <- zetaAnalysis[,c("meanLU","meanAL","meanDist","zeta_1","zeta_2","zeta_N","meanASCI","modeledASCI")]
  zetaCor <- rcorr(as.matrix(zetaCor),type="pearson")
  corr <- zetaCor$r
  p.mat <- zetaCor$P
  colnames(corr) <- c("Land Use","Altitude","Distance",":zeta[1]",":zeta[2]",":zeta[10]","Mean ASCI","Modeled ASCI")
  rownames(corr) <- c("Land Use","Altitude","Distance",":zeta[1]",":zeta[2]",":zeta[10]","Mean ASCI","Modeled ASCI")
  par(xpd=TRUE)
  png(paste("zetaAnalysisMorphologicalAlgae",taxonomicLevel,".png",sep=""),width=7,height=7,units="in",res=600)
  corrplot(corr = corr, p.mat = p.mat, diag = FALSE, type="lower", sig.level = 0.0001, tl.col="black", tl.srt=45, tl.cex=1.3, order="original",mar=c(0,0,3,0), cl.align.text = "r")
  mtext(paste("Morphologically sorted soft-bodied algae aggregated to",taxonomicLevel), at=2.5, line=3, cex=1.3)
  dev.off() 
}
