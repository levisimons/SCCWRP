library("plyr")
library(dplyr)
library("ggplot2")
library(ggpubr)
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)
library(MASS)
library(zetadiv)
library(magrittr)
library(relaimpo)
library(geosphere)

setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage values.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Remove duplicate measures.
GISBioData <- GISBioData[!duplicated(GISBioData[,c("UniqueID","FinalID","Count")]),]
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Merge in altitude.
GISBioData <- join(GISBioData,SCCWRP[,c("UniqueID","altitude")],by=c("UniqueID"))
#Add the number of genera per sample
taxaBySample <- count(GISBioData,UniqueID)
colnames(taxaBySample) <- c("UniqueID","nTaxa")
GISBioData <- join(GISBioData,taxaBySample,by=c("UniqueID"))
#Add the number of samples per watershed
WSBySample <- as.data.frame(table(SCCWRP$Watershed))
colnames(WSBySample) <- c("Watershed","NSamples")
GISBioData <- join(GISBioData,WSBySample,by=c("Watershed"))
#Read in watershed metadata
#HUC2 = State level, HUC4 = super-regional level, HUC6 = regional level, HUC8 = local level
Watersheds <- read.table("SCCWRPWatersheds.tsv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
Watersheds$Watershed <- gsub('\\s+', '', Watersheds$Watershed)
#Merge in watershed area
GISBioData <- join(GISBioData,Watersheds,by=c("Watershed"))
#Read in functional feeding group for each taxon.
#Abbreviations used in denoting functional feeding groups are as follows ( http://www.safit.org/Docs/CABW_std_taxonomic_effort.pdf ):
#P= predator MH= macrophyte herbivore OM= omnivore
#PA= parasite PH= piercer herbivore XY= xylophage (wood eater)
#CG= collector-gatherer SC= scraper
#CF= collector filterer SH= shredder 
FFG <- read.table("metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
# Filter data so only known functional feeding groups are kept.
FFG <- subset(FFG, FunctionalFeedingGroup != "")
# Generate functional feeding group data frame.
FFG <- FFG[,c("FinalID","LifeStageCode","FunctionalFeedingGroup")]
FFG <- subset(FFG,LifeStageCode=="L" | LifeStageCode=="X" | FinalID=="Hydrophilidae" | FinalID=="Hydraenidae")
FFG <- FFG[!duplicated(FFG$FinalID),]
#Merge in functional feeding groups into sample data.
GISBioData <- join(GISBioData,FFG[,c("FinalID","FunctionalFeedingGroup")],by=c("FinalID"))
FFGCounts <- na.omit(as.data.frame(unique(GISBioData$FunctionalFeedingGroup)))
colnames(FFGCounts) <- c("FunctionalFeedingGroup")
FFGCounts$FunctionalFeedingGroup <- as.character(as.factor(FFGCounts$FunctionalFeedingGroup))
FFGCounts <- arrange(FFGCounts,FunctionalFeedingGroup)
#Add column containing the sum of taxa, by functional feeding groups, within each sample.
tmp <- GISBioData[,c("UniqueID","FunctionalFeedingGroup","Count")]
colnames(tmp) <- c("UniqueID","FunctionalFeedingGroup","FFGCount")
tmp <- aggregate(tmp$FFGCount, by=list(Category=tmp$UniqueID,tmp$FunctionalFeedingGroup), FUN=sum)
colnames(tmp) <- c("UniqueID","FunctionalFeedingGroup","FFGCount")
tmp <- arrange(tmp,UniqueID,FunctionalFeedingGroup)
GISBioData <- join(GISBioData,tmp,by=c("UniqueID","FunctionalFeedingGroup"))

#Find watersheds with a larger enough set of samples for downstream analysis.
sampleMin <- 25 #Minimum of 25 samples per watershed
samplingNum <- 10 #Number of samples to select per sampling group within watershed.
GISBioDataLWS <- subset(GISBioData,NSamples>=sampleMin)
#Determine the number of samples a particular taxon is present in.
taxaSamples <- as.data.frame(table(GISBioDataLWS$FinalID))
colnames(taxaSamples) <- c("FinalID","NumSamples")
GISBioDataLWS <- join(GISBioDataLWS ,taxaSamples,by=c("FinalID"))
#Determine the number of watersheds a particular taxon is present in.
taxaWS <- as.data.frame(count(GISBioDataLWS ,FinalID,HUC8))
colnames(taxaWS) <- c("FinalID","HUC8","NumWatersheds")
taxaWS$NumWatersheds[taxaWS$NumWatersheds>0] <- 1
taxaWS <- as.data.table(table(taxaWS$FinalID))
colnames(taxaWS) <- c("FinalID","NumWatersheds")
GISBioDataLWS <- join(GISBioDataLWS,taxaWS,by=c("FinalID"))

zetaAnalysis <- data.frame()
#Determine zeta diversity per subsample group per watershed.  Run repeated subsampling
#to help with determinining significance of correlations between environmental factors
#and zeta diversity decay parameters.
set.seed(1)
for(j in 1:100){
  for(WS in unique(GISBioDataLWS$HUC8)){
    GISBioDataLocal <- subset(GISBioDataLWS,HUC8==WS) #Subsample by watershed.
    GISBioDataLocal <- GISBioDataLocal[GISBioDataLocal$UniqueID %in% sample(unique(GISBioDataLocal$UniqueID),samplingNum),] #Subsample to a uniform sample number.
    GISBioDataLocal <- GISBioDataLocal[!is.na(GISBioDataLocal$FunctionalFeedingGroup),]
    selected <- GISBioDataLocal
    metadata <- GISBioDataLocal[,c("UniqueID","LU_2000_5K","altitude","Longitude","Latitude","CSCI")]
    metadata <- metadata[!duplicated(metadata),] #Get unique environmental parameters per watershed set of samples.
    #Get geographic distances between samples
    MeanDist <- mean(distm(metadata[,c('Longitude','Latitude')], metadata[,c('Longitude','Latitude')], fun=distGeo))
    SDDist <- sd(distm(metadata[,c('Longitude','Latitude')], metadata[,c('Longitude','Latitude')], fun=distGeo))
    #Get all unique taxa in statewide data set.
    uniqueTaxa <- as.data.frame(unique(selected$FinalID))
    colnames(uniqueTaxa) <- c("FinalID")
    uniqueTaxa <- arrange(uniqueTaxa,FinalID)
    #Get all unique taxa in statewide data set.
    uniqueFFG <- as.data.frame(unique(selected$FunctionalFeedingGroup))
    colnames(uniqueFFG) <- c("FunctionalFeedingGroup")
    uniqueFFG <- arrange(uniqueFFG,FunctionalFeedingGroup)
    #Create presence/absence matrix of taxa in samples.
    #Rows for sample ID and columns 
    PresenceAbsence <- uniqueTaxa
    #Create presence/absence matrix of FFGs in samples.
    #Rows for sample ID and columns
    PAFFG <- uniqueFFG
    for(ID in unique(selected$UniqueID)){
      #Presence/Absence matrix for taxa.
      sampleDF <- subset(selected,UniqueID == ID)
      sampleDF <- sampleDF[,c("FinalID","Count")]
      tmp <- merge(sampleDF,uniqueTaxa,by=c("FinalID"),all=TRUE)
      colnames(tmp) <- c("FinalID",ID)
      PresenceAbsence <- cbind(PresenceAbsence,tmp[,c(2)])
      colnames(PresenceAbsence)[ncol(PresenceAbsence)] <- ID
      #Presence/Absence matrix for FFGs.
      sampleDF <- subset(selected,UniqueID == ID)
      sampleDF <- sampleDF[,c("FunctionalFeedingGroup","FFGCount")]
      sampleDF <- sampleDF[!duplicated(sampleDF),]
      tmp <- merge(sampleDF,uniqueFFG,by=c("FunctionalFeedingGroup"),all=TRUE)
      colnames(tmp) <- c("FunctionalFeedingGroup",ID)
      PAFFG <- cbind(PAFFG,tmp[,c(2)])
      colnames(PAFFG)[ncol(PAFFG)] <- ID
    }
    #Generate a presence/absence dataframe for zeta diversity analysis of taxa.
    #Rows for samples, columns for taxa IDs.
    PresenceAbsence[is.na(PresenceAbsence)] <- 0
    PresenceAbsence[PresenceAbsence > 0] <- 1
    data.SCCWRP <- as.data.frame(t(PresenceAbsence[,-c(1)]))
    colnames(data.SCCWRP) <- uniqueTaxa$FinalID
    #Generate a presence/absence dataframe for zeta diversity analysis of randomly assigned genera.
    #Rows for samples, columns for genera IDs.
    data.rand.SCCWRP <- as.data.frame(matrix(sample(unlist(data.SCCWRP),ncol(data.SCCWRP)*nrow(data.SCCWRP)),nrow=nrow(data.SCCWRP), byrow=T))
    colnames(data.rand.SCCWRP) <- colnames(data.SCCWRP)
    rownames(data.rand.SCCWRP) <- rownames(data.SCCWRP)
    #Generate a presence/absence dataframe for zeta diversity analysis of FFGs.
    #Rows for samples, columns for FFGs IDs.
    PAFFG[is.na(PAFFG)] <- 0
    PAFFG[PAFFG > 0] <- 1
    ffg.SCCWRP <- as.data.frame(t(PAFFG[,-c(1)]))
    colnames(ffg.SCCWRP) <- uniqueFFG$FunctionalFeedingGroup
    #Generate a presence/absence dataframe for zeta diversity analysis of randomly assigned FFGs.
    #Rows for samples, columns for FFGs IDs.
    ffg.rand.SCCWRP <- as.data.frame(matrix(sample(unlist(ffg.SCCWRP),ncol(ffg.SCCWRP)*nrow(ffg.SCCWRP)),nrow=nrow(ffg.SCCWRP), byrow=T))
    colnames(ffg.rand.SCCWRP) <- colnames(ffg.SCCWRP)
    rownames(ffg.rand.SCCWRP) <- rownames(ffg.SCCWRP)
    #Calculate zeta diversity decay for taxa within each watershed set of samples.
    zetaDecay <- Zeta.decline.ex(data.SCCWRP,orders=1:samplingNum,plot=FALSE,rescale=TRUE)
    ExpIntercept <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    ExpExp <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    ExpAIC <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    zeta_N <- Zeta.order.ex(data.SCCWRP,order=samplingNum,rescale=TRUE)$zeta.val #Higher order zeta diversity measure.
    PLExp <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    PLAIC <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    zeta_1 <- Zeta.order.ex(data.SCCWRP,order=1,rescale=TRUE) #Zeta order 1 for taxa.
    MeanAlpha <- zeta_1$zeta.val #Mean alpha diversity for taxa.
    SDAlpha <- zeta_1$zeta.val.sd #Standard deviation on alpha diversity for taxa.
    Beta <- mean(vegdist(data.SCCWRP,method="jaccard")) #Beta diversity for taxa based on Jaccard dissimilarity.
    #Calculate zeta diversity decay for randomly reassigned gener within each watershed set of samples.
    zetaDecay <- Zeta.decline.ex(data.rand.SCCWRP,orders=1:samplingNum,plot=FALSE,rescale=TRUE)
    ExpInterceptRand <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    ExpExpRand <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    ExpAICRand <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    zeta_NRand <- Zeta.order.ex(data.rand.SCCWRP,order=samplingNum,rescale=TRUE)$zeta.val #Higher order zeta diversity measure.
    PLExpRand <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    PLAICRand <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    zeta_1 <- Zeta.order.ex(data.rand.SCCWRP,order=1,rescale=TRUE) #Zeta order 1 for FFGs.
    MeanAlphaRand <- zeta_1$zeta.val #Mean alpha diversity for FFGs.
    SDAlphaRand <- zeta_1$zeta.val.sd #Standard deviation on alpha diversity for FFGs.
    BetaRand <- mean(vegdist(data.rand.SCCWRP,method="jaccard")) #Beta diversity for genera based on Jaccard dissimilarity.
    #Calculate zeta diversity decay for FFGs within each watershed set of samples.
    zetaDecay <- Zeta.decline.ex(ffg.SCCWRP,orders=1:samplingNum,plot=FALSE,rescale=TRUE)
    ExpInterceptFFG <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    ExpExpFFG <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    ExpAICFFG <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    zeta_NFFG <- Zeta.order.ex(ffg.SCCWRP,order=samplingNum,rescale=TRUE)$zeta.val #Higher order zeta diversity measure.
    PLExpFFG <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    PLAICFFG <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    zeta_1 <- Zeta.order.ex(ffg.SCCWRP,order=1,rescale=TRUE) #Zeta order 1 for FFGs.
    MeanAlphaFFG <- zeta_1$zeta.val #Mean alpha diversity for FFGs.
    SDAlphaFFG <- zeta_1$zeta.val.sd #Standard deviation on alpha diversity for FFGs.
    BetaFFG <- mean(vegdist(ffg.SCCWRP,method="jaccard")) #Beta diversity for FFGs based on Jaccard dissimilarity.
    #Calculate zeta diversity decay for randomly reassigned FFGs within each watershed set of samples.
    zetaDecay <- Zeta.decline.ex(ffg.rand.SCCWRP,orders=1:samplingNum,plot=FALSE,rescale=TRUE)
    ExpInterceptFFGRand <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
    ExpExpFFGRand <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
    ExpAICFFGRand <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
    zeta_NFFGRand <- Zeta.order.ex(ffg.rand.SCCWRP,order=samplingNum,rescale=TRUE)$zeta.val #Higher order zeta diversity measure.
    PLExpFFGRand <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
    PLAICFFGRand <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
    zeta_1 <- Zeta.order.ex(ffg.rand.SCCWRP,order=1,rescale=TRUE) #Zeta order 1 for FFGs.
    MeanAlphaFFGRand <- zeta_1$zeta.val #Mean alpha diversity for FFGs.
    SDAlphaFFGRand <- zeta_1$zeta.val.sd #Standard deviation on alpha diversity for FFGs.
    BetaFFGRand <- mean(vegdist(ffg.rand.SCCWRP,method="jaccard")) #Beta diversity for FFGs based on Jaccard dissimilarity.
    row <- t(as.data.frame(c(WS,mean(metadata$CSCI),mean(metadata$LU_2000_5K),sd(metadata$LU_2000_5K),mean(metadata$altitude),sd(metadata$altitude),MeanDist,SDDist,MeanAlpha,SDAlpha,Beta,ExpIntercept,ExpExp,ExpAIC,zeta_N,PLExp,PLAIC,MeanAlphaFFG,SDAlphaFFG,BetaFFG,ExpInterceptFFG,ExpExpFFG,ExpAICFFG,zeta_NFFG,PLExpFFG,PLAICFFG,MeanAlphaFFGRand,SDAlphaFFGRand,BetaFFGRand,ExpInterceptFFGRand,ExpExpFFGRand,ExpAICFFGRand,zeta_NFFGRand,PLExpFFGRand,PLAICFFGRand,MeanAlphaRand,SDAlphaRand,BetaRand,ExpInterceptRand,ExpExpRand,ExpAICRand,zeta_NRand,PLExpRand,PLAICRand)))
    zetaAnalysis <- rbind(zetaAnalysis,row)
    print(paste(j,WS))
  }
}
colnames(zetaAnalysis) <- c("Watershed","MeanCSCI","MeanLU","SDLU","MeanAltitude","SDAltitude","MeanDist","SDDist","MeanAlpha","SDAlpha","Beta","ExpIntercept","ExpExp","ExpAIC","zeta_N","PLExp","PLAIC","MeanAlphaFFG","SDAlphaFFG","BetaFFG","ExpInterceptFFG","ExpExpFFG","ExpAICFFG","zeta_NFFG","PLExpFFG","PLAICFFG","MeanAlphaFFGRand","SDAlphaFFGRand","BetaFFGRand","ExpInterceptFFGRand","ExpExpFFGRand","ExpAICFFGRand","zeta_NFFGRand","PLExpFFGRand","PLAICFFGRand","MeanAlphaRand","SDAlphaRand","BetaRand","ExpInterceptRand","ExpExpRand","ExpAICRand","zeta_NRand","PLExpRand","PLAICRand")
zetaAnalysis[,1:ncol(zetaAnalysis)] <- sapply(zetaAnalysis[,1:ncol(zetaAnalysis)],as.character)
zetaAnalysis[,2:ncol(zetaAnalysis)] <- sapply(zetaAnalysis[,2:ncol(zetaAnalysis)],as.numeric)
rownames(zetaAnalysis) <- 1:nrow(zetaAnalysis)
write.table(zetaAnalysis,"zetaAnalysis100PermutationsHUC8Dist.txt",quote=FALSE,sep="\t",row.names = FALSE)

##If the zeta diversity tables are already generated, read data in here for further analysis.
zetaAnalysis <- read.table("zetaAnalysis100PermutationsHUC8Dist.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

zetaModel <- lm(MeanCSCI ~ MeanAlpha+MeanAlphaFFG+Beta+BetaFFG+PLExp+PLExpFFG+zeta_N+zeta_NFFG, data=zetaAnalysis)
printCoefmat(coef(summary(step(zetaModel)))) #Determine parameters to drop from model via AIC scores.
zetaModel2 <- lm(MeanCSCI ~ MeanAlpha+MeanAlphaFFG+Beta+BetaFFG+zeta_N+zeta_NFFG, data=zetaAnalysis)# Drop one parameter by p.
printCoefmat(coef(summary(step(zetaModel2))))

#Model summaries
summary(zetaModel)
zetaModel$coefficients
anova(zetaModel)
calc.relimp(zetaModel,type="lmg",rela=TRUE)
# Model 2 is for publication
summary(zetaModel2)
zetaModel2$coefficients
anova(zetaModel2)
calc.relimp(zetaModel2,type="lmg",rela=TRUE)


#calculate a linear model fit.
zetaAnalysis$fit1 <- zetaModel$coefficients[1]+zetaAnalysis$MeanAlpha*zetaModel$coefficients[2]+zetaAnalysis$MeanAlphaFFG*zetaModel$coefficients[3]+zetaAnalysis$Beta*zetaModel$coefficients[4]+zetaAnalysis$BetaFFG*zetaModel$coefficients[5]+zetaAnalysis$PLExp*zetaModel$coefficients[6]+zetaAnalysis$PLExpFFG*zetaModel$coefficients[7]+zetaAnalysis$zeta_N*zetaModel$coefficients[8]+zetaAnalysis$zeta_NFFG*zetaModel$coefficients[9]
zetaAnalysis$fit2 <- zetaModel2$coefficients[1]+zetaAnalysis$MeanAlpha*zetaModel2$coefficients[2]+zetaAnalysis$MeanAlphaFFG*zetaModel2$coefficients[3]+zetaAnalysis$Beta*zetaModel2$coefficients[4]+zetaAnalysis$BetaFFG*zetaModel2$coefficients[5]+zetaAnalysis$zeta_N*zetaModel2$coefficients[6]+zetaAnalysis$zeta_NFFG*zetaModel2$coefficients[7]

#To generate map of data for a given environmental parameter in California.
library(ggmap)
library(maps)
library(mapview)
library(mapdata)
library(munsell)
library(leaflet)
library(devtools)
library(webshot)
dev.off()
#Aggregate mean health indices by watershed for mapping.
zetaMean <- zetaAnalysis[,c("Watershed","MeanCSCI","fit2","zeta_N","zeta_NFFG","Beta","BetaFFG","MeanAlpha","MeanAlphaFFG")]
zetaMean <- aggregate(zetaMean[,2:9],list(zetaMean$Watershed),mean)
colnames(zetaMean) <- c("HUC8","MeanCSCI","ModelIndex","zeta_10_taxa","zeta_10_functional","beta_taxa","beta_functional","alpha_taxa","alpha_functional")
MapCoordinates <- join(GISBioDataLWS,zetaMean,by=c("HUC8"))
MapCoordinates <- MapCoordinates[,c("LU_2000_5K","altitude","MeanCSCI","ModelIndex","zeta_10_taxa","zeta_10_functional","beta_taxa","beta_functional","alpha_taxa","alpha_functional","Longitude","Latitude")]
colnames(MapCoordinates) = c('LandUse','Altitude','MeanCSCI','ModelIndex','zeta_10_taxa','zeta_10_functional',"beta_taxa","beta_functional","alpha_taxa","alpha_functional",'longitude','latitude')
MapCoordinates <- MapCoordinates[!duplicated(MapCoordinates),]
MapCoordinates <- na.omit(MapCoordinates)
#Map data.
CalMap = leaflet(MapCoordinates) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=rainbow(10),domain=MapCoordinates$ModelIndex)
CalMap %>% addCircleMarkers(color = ~ColorScale(ModelIndex), fill = TRUE,radius=0.1,fillOpacity = 0.1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  leaflet::addLegend(position="topright", pal=ColorScale,values=~ModelIndex,title="Mean Modeled Index")

#Plotting indices against each other.
plot(zetaAnalysis$MeanCSCI,zetaAnalysis$fit1,xlab = "Mean CSCI score per sample group", ylab="Fitted stream conditions index")
abline(lm(fit1~MeanCSCI,data=zetaAnalysis),col="red")
legend("topleft", bty="n", legend=paste("R2 =", format(modelCor$estimate, digits=4),"\np < 2.2e-16"))

#Plotting taxonomic alpha, beta, and zeta diversity
dev.off()
par(mfrow=c(3,3))
plot(zetaAnalysis$MeanAltitude,zetaAnalysis$MeanAlpha,xlab="AL (m)",ylab=expression("Mean" ~ alpha[taxa]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(a)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(MeanAlpha~MeanAltitude,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanAltitude,zetaAnalysis$Beta,xlab="AL (m)",ylab=expression("Mean" ~ beta[taxa]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(b)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(Beta~MeanAltitude,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanAltitude,zetaAnalysis$zeta_N,xlab="AL (m)",ylab=expression("Mean" ~ zeta["10,taxa"]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(c)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(zeta_N~MeanAltitude,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanLU,zetaAnalysis$MeanAlpha,xlab="LU (%cover)",ylab=expression("Mean" ~ alpha[taxa]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(d)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(MeanAlpha~MeanLU,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanLU,zetaAnalysis$Beta,xlab="LU (%cover)",ylab=expression("Mean" ~ beta[taxa]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(e)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(Beta~MeanLU,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanLU,zetaAnalysis$zeta_N,xlab="LU (%cover)",ylab=expression("Mean" ~ zeta["10,taxa"]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(f)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(zeta_N~MeanLU,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanDist,zetaAnalysis$MeanAlpha,xlab="MD (m)",ylab=expression("Mean" ~ alpha[taxa]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(g)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(MeanAlpha~MeanDist,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanDist,zetaAnalysis$Beta,xlab="MD (m)",ylab=expression("Mean" ~ beta[taxa]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(h)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(Beta~MeanDist,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanDist,zetaAnalysis$zeta_N,xlab="MD (m)",ylab=expression("Mean" ~ zeta["10,taxa"]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(i)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(zeta_N~MeanDist,zetaAnalysis),col="red")

#Plotting functional alpha, beta, and zeta diversity
dev.off()
par(mfrow=c(3,3))
plot(zetaAnalysis$MeanAltitude,zetaAnalysis$MeanAlphaFFG,xlab="AL (m)",ylab=expression("Mean" ~ alpha[functional]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(a)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(MeanAlphaFFG~MeanAltitude,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanAltitude,zetaAnalysis$BetaFFG,xlab="AL (m)",ylab=expression("Mean" ~ beta[functional]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(b)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(BetaFFG~MeanAltitude,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanAltitude,zetaAnalysis$zeta_NFFG,xlab="AL (m)",ylab=expression("Mean" ~ zeta["10,functional"]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(c)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(zeta_NFFG~MeanAltitude,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanLU,zetaAnalysis$MeanAlphaFFG,xlab="LU (%cover)",ylab=expression("Mean" ~ alpha[functional]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(d)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(MeanAlphaFFG~MeanLU,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanLU,zetaAnalysis$BetaFFG,xlab="LU (%cover)",ylab=expression("Mean" ~ beta[functional]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(e)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(BetaFFG~MeanLU,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanLU,zetaAnalysis$zeta_NFFG,xlab="LU (%cover)",ylab=expression("Mean" ~ zeta["10,functional"]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(f)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(zeta_NFFG~MeanLU,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanDist,zetaAnalysis$MeanAlphaFFG,xlab="MD (m)",ylab=expression("Mean" ~ alpha[functional]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(g)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(MeanAlphaFFG~MeanDist,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanDist,zetaAnalysis$BetaFFG,xlab="MD (m)",ylab=expression("Mean" ~ beta[functional]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(h)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(BetaFFG~MeanDist,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanDist,zetaAnalysis$zeta_NFFG,xlab="MD (m)",ylab=expression("Mean" ~ zeta["10,functional"]))
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(i)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(zeta_NFFG~MeanDist,zetaAnalysis),col="red")

#Plotting CSCI and modeled index
dev.off()
par(mfrow=c(3,2))
plot(zetaAnalysis$MeanAltitude,zetaAnalysis$MeanCSCI,xlab="AL (m)",ylab="Mean CSCI")
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(a)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(MeanCSCI~MeanAltitude,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanAltitude,zetaAnalysis$fit2,xlab="AL (m)",ylab="Mean Modeled Index")
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(b)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(fit2~MeanAltitude,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanLU,zetaAnalysis$MeanCSCI,xlab="LU (%cover)",ylab="Mean CSCI")
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(c)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(MeanCSCI~MeanLU,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanLU,zetaAnalysis$fit2,xlab="LU (%cover)",ylab="Mean Modeled Index")
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(d)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(fit2~MeanLU,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanDist,zetaAnalysis$MeanCSCI,xlab="MD (m)",ylab="Mean CSCI")
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(e)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(MeanCSCI~MeanDist,zetaAnalysis),col="red")
plot(zetaAnalysis$MeanDist,zetaAnalysis$fit2,xlab="MD (m)",ylab="Mean Modeled Index")
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], "(f)",    adj = c( 1, 1 ), col = "blue" )
abline(lm(fit2~MeanDist,zetaAnalysis),col="red")

#Regression between network parameters.
library(Hmisc)
library(corrplot)
library("PerformanceAnalytics")
#Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
chart.Correlation(zetaAnalysis[,c("MeanCSCI","MeanAltitude","MeanLU","MeanAlpha","MeanAlphaFFG","zeta_N","zeta_NFFG","Beta","BetaFFG","PLExp","PLExpFFG","fit2")], histogram=FALSE, method="pearson")

#Probabilities of detecting a member of a particular taxon or functional feeding group.
histDataFFG <- as.data.frame(table(GISBioDataLWS$FunctionalFeedingGroup))
histDataFFG$Freq <- histDataFFG$Freq/nrow(GISBioDataLWS)
histData <- as.data.frame(table(GISBioDataLWS$FinalID))
histData$Freq <- histData$Freq/length(unique(GISBioDataLWS$UniqueID))

#Investigate trends in the abundances of taxa by FFGs.
GISBioDataLWSFFG <- GISBioDataLWS[,c("UniqueID","LU_2000_5K","altitude","FunctionalFeedingGroup","FFGCount","ActualOrganismCount")]
GISBioDataLWSFFG <- GISBioDataLWSFFG[!duplicated(GISBioDataLWSFFG),]
GISBioDataLWSFFG$FFGRA <- GISBioDataLWSFFG$FFGCount/GISBioDataLWSFFG$ActualOrganismCount
cor.test(GISBioDataLWSFFG[GISBioDataLWSFFG$FunctionalFeedingGroup=="CG",]$LU_2000_5K,GISBioDataLWSFFG[GISBioDataLWSFFG$FunctionalFeedingGroup=="CG",]$FFGRA)
