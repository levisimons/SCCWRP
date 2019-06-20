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
networkAnalysis <- read.table("CooccurrenceAnalysis.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
networkModel <- lm(MeanCSCI ~ N+C+M+S+zeta,data=networkAnalysis)
summary(networkModel)
calc.relimp(networkModel,type="lmg",rela=FALSE)
networkAnalysis$ModeledCSCI <- networkModel$coefficients[1]+networkModel$coefficients[2]*networkAnalysis$N+networkModel$coefficients[3]*networkAnalysis$C+networkModel$coefficients[4]*networkAnalysis$M+networkModel$coefficients[5]*networkAnalysis$S+networkModel$coefficients[6]*networkAnalysis$zeta
networkAnalysis$R_pos <- networkAnalysis$E_pos/networkAnalysis$E
networkAnalysis$K_normalized <- networkAnalysis$K/networkAnalysis$E
networkModel2 <- lm(MeanCSCI ~ N+C+M+zeta,data=networkAnalysis)
summary(networkModel2)
calc.relimp(networkModel2,type="lmg",rela=FALSE)
networkAnalysis$ModeledCSCI2 <- networkModel2$coefficients[1]+networkModel2$coefficients[2]*networkAnalysis$N+networkModel2$coefficients[3]*networkAnalysis$C+networkModel2$coefficients[4]*networkAnalysis$M+networkModel2$coefficients[5]*networkAnalysis$zeta

calc.relimp(lm(MeanCSCI ~ MeanLU+MeanAltitude+MeanDist,data=networkAnalysis))
calc.relimp(lm(ModeledCSCI2 ~ MeanLU+MeanAltitude+MeanDist,data=networkAnalysis))
