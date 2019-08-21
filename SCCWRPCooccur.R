#This script creates co-occurrence networks, and calculates some of their topological measures
#from presence/absence taxa by site data curated from the SCCWRP stream database.
#The approach used is described here: https://www.benjaminblonder.org/papers/2016_ECOG.pdf
require(netassoc)
require(igraph)
require(network)
require(plyr)
require(dplyr)
require(data.table)
require(geosphere)

setwd("~/Desktop/SCCWRP")
#setwd("/home/cmb-07/sn1/alsimons/SCCWRP")

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
taxaBySample <- as.data.frame(table(GISBioData$UniqueID))
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
FFGUnique <- as.data.frame(unique(FFG$FunctionalFeedingGroup))
colnames(FFGUnique) <- c("FunctionalFeedingGroup")
FFGUnique$FunctionalFeedingGroup <- as.character(FFGUnique$FunctionalFeedingGroup)
#Merge in functional feeding groups into sample data.
GISBioData <- join(GISBioData,FFG[,c("FinalID","FunctionalFeedingGroup")],by=c("FinalID"))
#Get the number of instances of taxa by functional feeding group by sample site.
FFGCount <- as.data.frame(table(GISBioData[,c("UniqueID","FunctionalFeedingGroup")]))
colnames(FFGCount) <- c("UniqueID","FunctionalFeedingGroup","FFGFreq")
FFGCount <- merge(FFGCount,FFGUnique,all=TRUE)
#Get the number of unique functional feeding groups by sample site.
FFGTypeCount <- as.data.frame(table(FFGCount[FFGCount$FFGFreq!=0,c("UniqueID")]))
colnames(FFGTypeCount) <- c("UniqueID","FFGTypeCount")
#Merge these counts back into the data frame.
GISBioData <- join(GISBioData,FFGCount,by=c("UniqueID","FunctionalFeedingGroup"))
#Merge in the functional feeding group type count back into the data frame.
GISBioData <- join(GISBioData,FFGTypeCount,by=c("UniqueID"))
#Find watersheds with a larger enough set of samples for downstream analysis.
sampleMin <- 15 #Minimum of samples per watershed
samplingNum <- 10 #Number of samples to select per sampling group within watershed.
GISBioDataLWS <- subset(GISBioData,NSamples>=sampleMin)
#Determine land use quantiles
GISBioDataLWS$LUquantile <- as.numeric(with(GISBioDataLWS, cut(LU_2000_5K, breaks=quantile(LU_2000_5K, probs=seq(0,1,0.2), na.rm=TRUE),include.lowest=TRUE)))
#Determine altitude quantiles
GISBioDataLWS$ALquantile <- as.numeric(with(GISBioDataLWS, cut(altitude, breaks=quantile(altitude, probs=seq(0,1,0.2), na.rm=TRUE),include.lowest=TRUE)))

networkAnalysis <- data.frame()

set.seed(1)
for(j in 1:100){
  for(WS in unique(GISBioDataLWS$HUC8)){
    for(LU in unique(GISBioDataLWS$LUquantile)){
      for(AL in unique(GISBioDataLWS$ALquantile)){
        GISBioDataLocal <- subset(GISBioDataLWS,HUC8==WS & LUquantile==LU & ALquantile==AL) #Subsample by watershed.
        if(length(unique(GISBioDataLocal$UniqueID))>=sampleMin){
          GISBioDataLocal <- GISBioDataLocal[GISBioDataLocal$UniqueID %in% sample(unique(GISBioDataLocal$UniqueID),samplingNum),] #Subsample to a uniform sample number.
          Watershed <- unique(GISBioDataLocal$Watershed)
          selected <- GISBioDataLocal
          metadata <- GISBioDataLocal[,c("UniqueID","LU_2000_5K","altitude","Longitude","Latitude","CSCI")]
          metadata <- metadata[!duplicated(metadata),] #Get unique environmental parameters per watershed set of samples.
          #Get geographic distances between samples
          MeanDist <- mean(distm(metadata[,c('Longitude','Latitude')], metadata[,c('Longitude','Latitude')], fun=distGeo))
          SDDist <- sd(distm(metadata[,c('Longitude','Latitude')], metadata[,c('Longitude','Latitude')], fun=distGeo))
          #Mean land use for samples
          MeanLU <- mean(metadata$LU_2000_5K)
          SDLU <- sd(metadata$LU_2000_5K)
          #Mean altitude for samples
          MeanAltitude <- mean(metadata$altitude)
          SDAltitude <- sd(metadata$altitude)
          #Mean CSCI for samples
          MeanCSCI <- mean(metadata$CSCI)
          SDCSCI <- sd(metadata$CSCI)
          #Get all unique taxa in statewide data set.
          uniqueTaxa <- as.data.frame(unique(selected$FinalID))
          colnames(uniqueTaxa) <- c("FinalID")
          uniqueTaxa <- arrange(uniqueTaxa,FinalID)
          #Create presence/absence matrix of taxa in samples.
          #Rows for sample ID and columns 
          PresenceAbsence <- uniqueTaxa
          for(ID in unique(selected$UniqueID)){
            #Presence/Absence matrix for taxa.
            sampleDF <- subset(selected,UniqueID == ID)
            sampleDF <- sampleDF[,c("FinalID","Count")]
            tmp <- merge(sampleDF,uniqueTaxa,by=c("FinalID"),all=TRUE)
            colnames(tmp) <- c("FinalID",ID)
            PresenceAbsence <- cbind(PresenceAbsence,tmp[,c(2)])
            colnames(PresenceAbsence)[ncol(PresenceAbsence)] <- ID
          }
          #Rows for samples, columns for taxa IDs.
          m_obs <- PresenceAbsence[,-c(1)]
          m_obs[is.na(m_obs)] <- 0
          m_obs[m_obs > 0] <- 1
          #Calculate co-occurrence network using presence/absence data.
          n <- make_netassoc_network(m_obs, vegan::permatfull(m_obs,fixedmar="both",mtype="prab",times=100)$perm[[1]],method="partial_correlation",args=list(method="shrinkage"),p.method='fdr', numnulls=1000, plot=FALSE,alpha=1e-4,verbose=FALSE)
          #Analyze the network's topology.
          networkgraph <- as.undirected(n$network_all,mode="collapse") #Generate full graph from co-occurrence patterns.
          networkgraph_pos <- as.undirected(n$network_pos,mode="collapse") #Generate positive edge graph from co-occurrence patterns.
          networkgraph_neg <- as.undirected(n$network_neg,mode="collapse") #Generate negative edge graph from co-occurrence patterns.
          if(gsize(networkgraph)>0 & gsize(networkgraph_pos)>0 & gsize(networkgraph_neg)>0){
            networkmatrix <- as.matrix(get.adjacency(networkgraph,attr='weight')) #Generate matrix from co-occurrence graph.
            lambda_1 <- Re(eigen(networkmatrix)$values[1]) #The real component of the most positive eigenvalue
            S <- mean(networkmatrix[networkmatrix!=0]) #Mean edge weight as classified by standardized effect size values
            M <- modularity(cluster_edge_betweenness(networkgraph, weights=NULL,directed=FALSE)) #Modularity
            L <- mean_distance(networkgraph,directed=FALSE) #Characteristic path length
            N <- gorder(networkgraph) #Number of nodes
            C <- edge_density(networkgraph) #Connectance
            E <- gsize(networkgraph) #Number of edges
            K <- (2*E)/N#Mean node degree
            M_rand <- (1-(2/sqrt(N)))*(2/(K))^(2/3) #Modularity of a corresponding random network
            L_rand <- 0.5+((log(N)-0.5772156649)/log(K))
            #Calculate the degree heterogeneity.
            networkmatrix[upper.tri(networkmatrix)] <- 0
            networkmatrix <- ifelse(networkmatrix!=0,1,networkmatrix)
            zeta <- mean(colSums(networkmatrix)^2)/mean(colSums(networkmatrix))^2
            #Calculate the degree heterogeneity of a corresponding random network.
            networkmatrix_rand <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
            zeta_rand <- mean(colSums(networkmatrix_rand)^2)/mean(colSums(networkmatrix_rand))^2
            #Analyze the positive network's topology.
            networkgraph <- as.undirected(n$network_pos,mode="collapse") #Generate graph from co-occurrence patterns.
            #networkgraph <- as.undirected(graph_from_adjacency_matrix(as.matrix(n$matrix_spsp_ses_thresholded)),mode="collapse")
            networkmatrix <- as.matrix(get.adjacency(networkgraph,attr='weight')) #Generate matrix from co-occurrence graph.
            S_pos <- mean(networkmatrix[networkmatrix!=0]) #Mean edge weight as classified by standardized effect size values
            M_pos <- modularity(cluster_edge_betweenness(networkgraph, weights=NULL,directed=FALSE)) #Modularity
            L_pos <- mean_distance(networkgraph,directed=FALSE) #Characteristic path length
            N_pos <- gorder(networkgraph) #Number of nodes
            C_pos <- edge_density(networkgraph) #Connectance
            E_pos <- gsize(networkgraph) #Number of edges
            K_pos <- (2*E_pos)/N_pos#Mean node degree
            M_rand_pos <- (1-(2/sqrt(N_pos)))*(2/(K_pos))^(2/3) #Modularity of a corresponding random network
            L_rand_pos <- 0.5+((log(N_pos)-0.5772156649)/log(K_pos))
            #Calculate the degree heterogeneity.
            networkmatrix[upper.tri(networkmatrix)] <- 0
            networkmatrix <- ifelse(networkmatrix!=0,1,networkmatrix)
            zeta_pos <- mean(colSums(networkmatrix)^2)/mean(colSums(networkmatrix))^2
            #Calculate the degree heterogeneity of a corresponding random network.
            networkmatrix_rand <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
            zeta_rand_pos <- mean(colSums(networkmatrix_rand)^2)/mean(colSums(networkmatrix_rand))^2
            #Analyze the negative network's topology.
            networkgraph <- as.undirected(n$network_neg,mode="collapse") #Generate graph from co-occurrence patterns.
            #networkgraph <- as.undirected(graph_from_adjacency_matrix(as.matrix(n$matrix_spsp_ses_thresholded)),mode="collapse")
            networkmatrix <- as.matrix(get.adjacency(networkgraph,attr='weight')) #Generate matrix from co-occurrence graph.
            S_neg <- mean(networkmatrix[networkmatrix!=0]) #Mean edge weight as classified by standardized effect size values
            M_neg <- modularity(cluster_edge_betweenness(networkgraph, weights=NULL,directed=FALSE)) #Modularity
            L_neg <- mean_distance(networkgraph,directed=FALSE) #Characteristic path length
            N_neg <- gorder(networkgraph) #Number of nodes
            C_neg <- edge_density(networkgraph) #Connectance
            E_neg <- gsize(networkgraph) #Number of edges
            K_neg <- (2*E_neg)/N_neg#Mean node degree
            M_rand_neg <- (1-(2/sqrt(N_neg)))*(2/(K_neg))^(2/3) #Modularity of a corresponding random network
            L_rand_neg <- 0.5+((log(N_neg)-0.5772156649)/log(K_neg))
            #Calculate the degree heterogeneity.
            networkmatrix[upper.tri(networkmatrix)] <- 0
            networkmatrix <- ifelse(networkmatrix!=0,1,networkmatrix)
            zeta_neg <- mean(colSums(networkmatrix)^2)/mean(colSums(networkmatrix))^2
            #Calculate the degree heterogeneity of a corresponding random network.
            networkmatrix_rand <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
            zeta_rand_neg <- mean(colSums(networkmatrix_rand)^2)/mean(colSums(networkmatrix_rand))^2
            row <- t(as.data.frame(c(Watershed,MeanLU,SDLU,MeanAltitude,SDAltitude,MeanDist,SDDist,MeanCSCI,SDCSCI,N,C,E,M,L,zeta,zeta_rand,K,M_rand,L_rand,lambda_1,S,N_pos,C_pos,E_pos,M_pos,L_pos,zeta_pos,zeta_rand_pos,K_pos,M_rand_pos,L_rand_pos,S,N_neg,C_neg,E_neg,M_neg,L_neg,zeta_neg,zeta_rand_neg,K_neg,M_rand_neg,L_rand_neg,S_neg)))
            networkAnalysis <- rbind(networkAnalysis,row)
            print(paste(j,Watershed))
            print(row)
          }
        }
      }
    }
  }
}
colnames(networkAnalysis) <- c("Watershed","MeanLU","SDLU","MeanAltitude","SDAltitude","MeanDist","SDDist","MeanCSCI","SDCSCI","N","C","E","M","L","zeta","zeta_rand","K","M_rand","L_rand","lambda_1","S","N_pos","C_pos","E_pos","M_pos","L_pos","zeta_pos","zeta_rand_pos","K_pos","M_rand_pos","L_rand_pos","S_pos","N_neg","C_neg","E_neg","M_neg","L_neg","zeta_neg","zeta_rand_neg","K_neg","M_rand_neg","L_rand_neg","S_neg")
networkAnalysis[,1:ncol(networkAnalysis)] <- sapply(networkAnalysis[,1:ncol(networkAnalysis)],as.character)
networkAnalysis[,2:ncol(networkAnalysis)] <- sapply(networkAnalysis[,2:ncol(networkAnalysis)],as.numeric)
rownames(networkAnalysis) <- 1:nrow(networkAnalysis)
write.table(networkAnalysis,"CooccurrenceAnalysis.txt",quote=FALSE,sep="\t",row.names = FALSE)
#################################


##This part is used to run locally.
require(relaimpo)
require(vegan)
#First simple linear model of network parameters.
networkAnalysis <- read.table("CooccurrenceAnalysis.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
networkAnalysis <- networkAnalysis[!is.na(networkAnalysis$MeanCSCI),]
networkAnalysis$CoVLU <- networkAnalysis$SDLU/networkAnalysis$MeanLU
networkAnalysis$CoVAltitude <-  networkAnalysis$SDAltitude/networkAnalysis$MeanAltitude
#networkModel <- lm(MeanCSCI ~ N+C+S+M+zeta,data=networkAnalysis)
networkModel <- lm(MeanCSCI ~ C+M+S+zeta,data=networkAnalysis)
summary(networkModel)
anova(networkModel)
#adonis(MeanCSCI ~ N+C+S+M+zeta,data=networkAnalysis[!is.na(networkAnalysis$MeanCSCI),],strata = networkAnalysis[!is.na(networkAnalysis$MeanCSCI),"Watershed"],permutations=10,method="manhattan")
calc.relimp(networkModel,type="lmg",rela=FALSE)
networkAnalysis$ModeledCSCI <- networkModel$coefficients[1]+networkModel$coefficients[2]*networkAnalysis$N+networkModel$coefficients[3]*networkAnalysis$C+networkModel$coefficients[4]*networkAnalysis$S+networkModel$coefficients[5]*networkAnalysis$M+networkModel$coefficients[6]*networkAnalysis$zeta
printCoefmat(coef(summary(step(networkModel))))
#networkAnalysis$ModeledCSCI <- networkModel$coefficients[1]+networkModel$coefficients[2]*networkAnalysis$N+networkModel$coefficients[3]*networkAnalysis$C_pos+networkModel$coefficients[4]*networkAnalysis$M_pos+networkModel$coefficients[5]*networkAnalysis$S_pos+networkModel$coefficients[6]*networkAnalysis$zeta_pos

# Assessing R2 shrinkage using 10-Fold Cross-Validation 
require(bootstrap)
# define functions 
theta.fit <- function(x,y){lsfit(x,y)}
theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef}
#Run for zero horizon images.
# matrix of predictors
#X <- as.matrix(networkAnalysis[,c("N","C","M","S","zeta")])
X <- as.matrix(networkAnalysis[,c("C","M","S","zeta")])
# vector of predicted values
y <- as.matrix(networkAnalysis[,c("MeanCSCI")])
#Run cross-validation
results <- crossval(X,y,theta.fit,theta.predict,ngroup=10)
cor(y, networkModel$fitted.values)**2 # raw R2 
cor(y,results$cv.fit)**2 # cross-validated R2

#Evaluating the mean and modeled CSCI against environmental parameters.
calc.relimp(lm(MeanCSCI ~ MeanLU+SDLU+MeanAltitude+SDAltitude+MeanDist,data=networkAnalysis))
anova(lm(MeanCSCI ~ MeanLU+SDLU+MeanAltitude+SDAltitude+MeanDist,data=networkAnalysis))
calc.relimp(lm(ModeledCSCI ~ MeanLU+SDLU+MeanAltitude+SDAltitude+MeanDist,data=networkAnalysis))
anova(lm(ModeledCSCI ~ MeanLU+SDLU+MeanAltitude+SDAltitude+MeanDist,data=networkAnalysis))

#Investigate trends between number of genera per functional feeding groups per site and other factors.
FFGTrends <- GISBioDataLWS[,c("UniqueID","Watershed","LU_2000_5K","altitude","FunctionalFeedingGroup","FFGFreq","nTaxa","FFGTypeCount")]
FFGTrends <- FFGTrends[!duplicated(FFGTrends),]
FFGTrends <- merge(FFGUnique,FFGTrends,by=c("FunctionalFeedingGroup"),all.x=TRUE)
FFGTrends[is.na(FFGTrends)] <- 0
#Normalize the FFG frequency by the number of overall taxonomic groups per site.
FFGTrends$FFGFreq_normalized <- FFGTrends$FFGFreq/FFGTrends$nTaxa
#Compare FFG frequency by the number of overall taxonomic groups per site with altitude and land use.
summary(lm(FFGTrends[FFGTrends$FunctionalFeedingGroup=="OM","FFGFreq_normalized"]~FFGTrends[FFGTrends$FunctionalFeedingGroup=="OM","LU_2000_5K"]+FFGTrends[FFGTrends$FunctionalFeedingGroup=="OM","altitude"]))
#Decline in functional diversity with land use.
cor.test(FFGTrends$LU_2000_5K,FFGTrends$FFGTypeCount)

#Relations between functional diversity and co-occurrence network topology.
networkFFG <- join(networkAnalysis,FFGTrends,by=c("Watershed"))
cor.test(networkFFG[networkFFG$FunctionalFeedingGroup=="OM","FFGFreq_normalized"],networkFFG[networkFFG$FunctionalFeedingGroup=="OM","C"])
cor.test(networkFFG[networkFFG$FunctionalFeedingGroup=="CG","FFGFreq_normalized"],networkFFG[networkFFG$FunctionalFeedingGroup=="CG","S"])
anova(lm(networkFFG[networkFFG$FunctionalFeedingGroup=="OM","C"]~networkFFG[networkFFG$FunctionalFeedingGroup=="OM","FFGFreq_normalized"]+networkFFG[networkFFG$FunctionalFeedingGroup=="OM","Watershed"]))

#Plot network patterns
require(ggplot2)
require(viridis)
networkPlot <- ggplot(networkAnalysis, aes(x=MeanCSCI,y=ModeledCSCI,color=MeanLU))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=ModeledCSCI))
networkPlot+xlab("Mean CSCI")+ylab("Modeled CSCI")+scale_color_gradientn("Land use\n(% cover)",colours = rev(plasma(10)),limits=c(0,100))
#
networkPlot <- ggplot(networkAnalysis, aes(x=MeanCSCI,y=ModeledCSCI,color=log10(MeanAltitude)))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=ModeledCSCI))
networkPlot+xlab("Mean CSCI")+ylab("Modeled CSCI")+scale_color_gradientn("Altitude (m)",colours = rev(plasma(10)))
#
networkPlot <- ggplot(networkAnalysis, aes(x=MeanCSCI,y=ModeledCSCI,color=log10(MeanDist)))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=ModeledCSCI))
networkPlot+xlab("Mean CSCI")+ylab("Modeled CSCI")+scale_color_gradientn("Distance (m)",colours = rev(plasma(10)))
#
