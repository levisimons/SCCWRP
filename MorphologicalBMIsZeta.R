rm(list=ls())
require(dplyr)
require(zetadiv)
require(geosphere)
require(relaimpo)

#wd <- "~/Desktop/SCCWRP/MorphologicalBMIsAlgae/"
wd <- "/project/noujdine_61/alsimons/SCCWRP"
setwd(wd)

#Read in morphological stream community data describing all BMIs
#and format them as presence/absence tables.
BMIInput <- read.table("BugTaxa.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
#Standarize date format.
BMIInput$sampledate <- as.Date(BMIInput$sampledate,format="%Y-%m-%d")
BMIInput$sampledate <- format(BMIInput$sampledate,format="%m/%d/%y")
#Create unique sample identifier.
BMIInput$UniqueID <- paste(BMIInput$stationcode,BMIInput$sampledate)
#Subset columns of interest
BMIInput <- BMIInput[,c("stationcode","sampledate","replicate","UniqueID","finalid")]

#Add in additional BMI morphology data.
BMIInput2 <- read.table("BugTaxa2.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
#Standarize date format.
BMIInput2$sampledate <- as.Date(BMIInput2$sampledate,format="%m/%d/%y")
BMIInput2$sampledate <- format(BMIInput2$sampledate,format="%m/%d/%y")
#Create unique sample identifier.
BMIInput2$UniqueID <- paste(BMIInput2$stationcode,BMIInput2$sampledate)
#Subset columns of interest
BMIInput2 <- BMIInput2[,c("stationcode","sampledate","replicate","UniqueID","finalid")]

#Create merged algal morphology data set.
BMIInput <- rbind(BMIInput,BMIInput2)

#Generated here: https://github.com/nuzhdinlab/SCCWRP/blob/master/MorphologicalTaxonomyGenerator.R
uniqueBMIs <- read.table("BMITaxonomiesMorphological.txt", header=T, sep="\t",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")

#Create merged data frame with sample data and full taxonomic information.
communityInput <- dplyr::left_join(uniqueBMIs,BMIInput[,c("stationcode","sampledate","finalid")],by=c("LeafTaxa"="finalid"))

#Standarize date format.
#communityInput$sampledate <- as.Date(communityInput$sampledate,format="%Y-%m-%d")
#communityInput$sampledate <- format(communityInput$sampledate,format="%m/%d/%y")

#Choose a taxonomic level to group count data by.
#Levels are domain, kingdom, phylum, class, order, family, genus, species, OTUID
taxonomicLevels <- colnames(uniqueBMIs[,grep("^[A-Za-z]", colnames(uniqueBMIs))])
taxonomicLevel <- c("order") #Choose a taxonomic level to aggregate count data on.
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

#Remove duplicate entries.
sampleIDs <- sampleIDs[,c("StationCode","Date","UniqueID")]
sampleIDs <- sampleIDs[!duplicated(sampleIDs),]

#Remove unnecessary sample columns.
communityInput <- communityInput[,c("stationcode","sampledate",taxonomicLevel)]

#Remove empty taxonomic rows.
communityInput <- communityInput[!is.na(communityInput[,taxonomicLevel]),]

#Create unique identifier column.
communityInput$UniqueID <- paste(communityInput$stationcode,communityInput$sampledate)

#Only keep community data for samples which have both morphological and metagenomic data.
communityInput <- communityInput[communityInput$UniqueID %in% sampleIDs$UniqueID,]

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

#Read in california stream condition index data.
CSCI <- read.table("CSCI.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")

#Standardize the date format.
CSCI$sampledate <- as.Date(CSCI$sampledate,format="%Y-%m-%d")
CSCI$sampledate <- format(CSCI$sampledate,format="%m/%d/%y")
#Create unique identifier column.
CSCI$UniqueID <- paste(CSCI$stationcode,CSCI$sampledate)

#Filter CSCI data frame to contain only samples shared between morphological and metagenomic methods.
CSCI <- CSCI[CSCI$UniqueID %in% communityInput$UniqueID,]

#Read in sample station metadata.
metadata <- read.table("GIS.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")

#Filter metadata to only contain sites shared between morphological and metagenomic methods.
metadata <- metadata[metadata$stationcode %in% communityInput$stationcode,]

#Create upstream land use column
metadata$LU <- metadata$ag_2011_5k+metadata$urban_2011_5k+metadata$code_21_2011_5k

#Add land use bands.  Band1 (Low): 0-3%, Band2 (Intermediate): 3-15%, Band3 (High): 15-100%.
metadata$LUBand <- case_when(metadata$LU <= 3 ~ "1", metadata$LU > 3 & metadata$LU <= 15 ~ "2", metadata$LU > 15 ~ "3", TRUE ~ as.character(metadata$LU))

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
    data.spec <- communityInputSummarized[grepl(paste(metadataSubset$stationcode, collapse="|"),rownames(communityInputSummarized)),]
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
    meanCSCI <- mean(CSCI[CSCI$UniqueID %in% rownames(data.spec),"csci"])
    sdCSCI <- sd(CSCI[CSCI$UniqueID %in% rownames(data.spec),"csci"])
    #Aggregate zeta analysis output
    print(paste(LUBand,meanLU,sdLU,meanAL,sdAL,meanDist,sdDist,meanCSCI,sdCSCI,zeta_N,zeta_Nsd,zeta_1,zeta_1sd,zeta_2,zeta_2sd,ExpExp,ExpAIC,PLExp,PLAIC))
    dataRow <- t(as.data.frame(list(c(LUBand,meanLU,sdLU,meanAL,sdAL,meanDist,sdDist,meanCSCI,sdCSCI,zeta_N,zeta_Nsd,zeta_1,zeta_1sd,zeta_2,zeta_2sd,ExpExp,ExpAIC,PLExp,PLAIC))))
    rownames(dataRow) <- NULL
    zetaAnalysis <- rbind(zetaAnalysis,dataRow)
  }
}
colnames(zetaAnalysis) <- c("LUBand","meanLU","sdLU","meanAL","sdAL","meanDist","sdDist","meanCSCI","sdCSCI","zeta_N","zeta_Nsd","zeta_1","zeta_1sd","zeta_2","zeta_2sd","ExpExp","ExpAIC","PLExp","PLAIC")
indx <- sapply(zetaAnalysis, is.factor)
zetaAnalysis[indx] <- lapply(zetaAnalysis[indx], function(x) as.numeric(as.character(x)))
zetaAnalysis <- do.call(data.frame,lapply(zetaAnalysis, function(x) replace(x, is.infinite(x),NA)))
#Save zeta diversity analysis for a given taxonomic level.
write.table(zetaAnalysis,paste("zetaAnalysisMorphologicalBMIs",taxonomicLevel,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#

##To run locally.
taxonomicLevel <- "order"
zetaAnalysis <- read.table(paste("zetaAnalysisMorphologicalBMIs",taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

zetaModel <- lm(meanCSCI~zeta_1+zeta_2+zeta_N,data=zetaAnalysis)
plot(zetaModel$model$meanCSCI,zetaModel$fitted.values)
cor.test(zetaModel$model$meanCSCI,zetaModel$fitted.values)
zetaAnalysis$modeledCSCI <- zetaModel$fitted.values
calc.relimp(zetaModel)
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(zetaModel)
dev.off()

# Assessing R2 shrinkage using 10-Fold Cross-Validation 
#require(bootstrap)
require(caret)
set.seed(1)
train.control <- trainControl(method="repeatedcv",number=10,repeats=10)
CSCImodel <- train(meanCSCI~zeta_1+zeta_2+zeta_N,data=zetaAnalysis,method="lm",trControl=train.control)
print(cor.test(zetaModel$model$meanCSCI,zetaModel$fitted.values))
print(CSCImodel)

#Comparing mean and modeled CSCI versus environmental parameters
zetaVarModel1 <- lm(modeledCSCI~meanAL+meanLU+meanDist,data=zetaAnalysis)
summary(zetaVarModel1)
calc.relimp(zetaVarModel1)
anova(zetaVarModel1)
zetaVarModel2 <- lm(meanCSCI~meanAL+meanLU+meanDist,data=zetaAnalysis)
summary(zetaVarModel2)
calc.relimp(zetaVarModel2)
anova(zetaVarModel2)

#Correlation plots of mean and modeled CSCI scores, zeta diversity measures, and environmental parameters.
require(Hmisc)
require(corrplot)
communityType <- "BMIs"
taxonomicLevel <- "species"
zetaAnalysis <- read.table(paste("zetaAnalysisMorphologicalBMIs",taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
zetaModel <- lm(meanCSCI~zeta_1+zeta_2+zeta_N,data=zetaAnalysis)
zetaAnalysis$modeledCSCI <- zetaModel$fitted.values
zetaCor <- zetaAnalysis[,c("meanLU","meanAL","meanDist","zeta_1","zeta_2","zeta_N","meanCSCI","modeledCSCI")]
zetaCor <- rcorr(as.matrix(zetaCor),type="pearson")
corr <- zetaCor$r
p.mat <- zetaCor$P
colnames(corr) <- c("Land Use","Altitude","Distance",":zeta[1]",":zeta[2]",":zeta[10]","Mean CSCI","Modeled CSCI")
rownames(corr) <- c("Land Use","Altitude","Distance",":zeta[1]",":zeta[2]",":zeta[10]","Mean CSCI","Modeled CSCI")
par(xpd=TRUE)
corrplot(corr = corr, p.mat = p.mat, diag = FALSE, type="lower", sig.level = 0.0001, tl.col="black", tl.srt=45, tl.cex=1.3, order="original",mar=c(0,0,3,0), cl.align.text = "r")
mtext(paste("Morphologically sorted",communityType,"aggregated to",taxonomicLevel), at=2.5, line=3, cex=1.3)
