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

#Find watersheds with a larger enough set of samples for downstream analysis.
sampleMin <- 15 #Minimum of 25 samples per watershed
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

networkAnalysis <- data.frame()

set.seed(1)
for(j in 1:1){
  for(WS in unique(GISBioDataLWS$HUC8)){
    GISBioDataLocal <- subset(GISBioDataLWS,HUC8==WS) #Subsample by watershed.
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
    networkgraph <- as.undirected(n$network_all,mode="collapse") #Generate graph from co-occurrence patterns.
    networkmatrix <- as.matrix(get.adjacency(networkgraph,attr='weight')) #Generate matrix from co-occurrence graph.
    M <- modularity(cluster_edge_betweenness(networkgraph, weights=NULL,directed=FALSE)) #Modularity
    L <- mean_distance(networkgraph,directed=FALSE) #Characteristic path length
    Network <- as.network(get.adjacency(networkgraph)) #Generate network from graph.
    N <- network.size(Network) #Number of nodes
    C <- network.density(Network) #Connectance
    E <- network.edgecount(Network) #Number of edges
    #Calculate the degree heterogeneity.
    networkmatrix[upper.tri(networkmatrix)] <- 0
    networkmatrix <- ifelse(networkmatrix!=0,1,networkmatrix)
    zeta <- mean(colSums(networkmatrix)^2)/mean(colSums(networkmatrix))^2
    row <- t(as.data.frame(c(Watershed,MeanLU,SDLU,MeanAltitude,SDAltitude,MeanDist,SDDist,MeanCSCI,SDCSCI,N,C,E,M,L,zeta)))
    networkAnalysis <- rbind(networkAnalysis,row)
    print(paste(j,Watershed))
  }
}
colnames(networkAnalysis) <- c("Watershed","MeanLU","SDLU","MeanAltitude","SDAltitude","MeanDist","SDDist","MeanCSCI","SDCSCI","N","C","E","M","L","zeta")
networkAnalysis[,1:ncol(networkAnalysis)] <- sapply(networkAnalysis[,1:ncol(networkAnalysis)],as.character)
networkAnalysis[,2:ncol(networkAnalysis)] <- sapply(networkAnalysis[,2:ncol(networkAnalysis)],as.numeric)
rownames(networkAnalysis) <- 1:nrow(networkAnalysis)
write.table(networkAnalysis,"CooccurrenceAnalysis.txt",quote=FALSE,sep="\t",row.names = FALSE)

networkAnalysis <- read.table("CooccurrenceAnalysis.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
networkModel <- lm(MeanCSCI ~ N+C+E+M+L+zeta,data=networkAnalysis)
summary(networkModel)
calc.relimp(networkModel,type="lmg",rela=FALSE)
networkAnalysis$ModeledCSCI <- networkModel$coefficients[1]+networkModel$coefficients[2]*networkAnalysis$N+networkModel$coefficients[3]*networkAnalysis$C+networkModel$coefficients[4]*networkAnalysis$E+networkModel$coefficients[5]*networkAnalysis$M+networkModel$coefficients[6]*networkAnalysis$L+networkModel$coefficients[7]*networkAnalysis$zeta
