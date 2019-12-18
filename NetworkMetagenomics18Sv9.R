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
require(netassoc)
require(igraph)
require(network)

wd <- "/home/cmb-07/sn1/alsimons/SCCWRP"
wd <- "~/Desktop/SCCWRP/Metagenomics/"
setwd(wd)

#Choose a community type.
#"" for all data.
#"algae" for algal communities.
#"BMIs" for BMI communities
communityType <- "algae"

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
taxonomicLevel <- c("family") #Choose a taxonomic level to aggregate count data on.
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
networkAnalysis <- data.frame()
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
    m_obs <- as.data.frame(t(data.spec[rowSums(data.spec!=0) > 0,colSums(data.spec!=0) > 0]))
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
      row <- t(as.data.frame(c(i,meanLU,sdLU,meanAL,sdAL,meanDist,sdDist,meanP,sdP,meanOrthoP,sdOrthoP,meanN,sdN,N,C,E,M,L,zeta,zeta_rand,K,M_rand,L_rand,lambda_1,S,N_pos,C_pos,E_pos,M_pos,L_pos,zeta_pos,zeta_rand_pos,K_pos,M_rand_pos,L_rand_pos,S,N_neg,C_neg,E_neg,M_neg,L_neg,zeta_neg,zeta_rand_neg,K_neg,M_rand_neg,L_rand_neg,S_neg)))
      networkAnalysis <- rbind(networkAnalysis,row)
      print(paste(j,i))
    }
  }
}
colnames(networkAnalysis) <- c("Ecoregion","meanLU","sdLU","meanAL","sdAL","meanDist","sdDist","meanP","sdP","meanOrthoP","sdOrthoP","meanN","sdN","N","C","E","M","L","zeta","zeta_rand","K","M_rand","L_rand","lambda_1","S","N_pos","C_pos","E_pos","M_pos","L_pos","zeta_pos","zeta_rand_pos","K_pos","M_rand_pos","L_rand_pos","S_pos","N_neg","C_neg","E_neg","M_neg","L_neg","zeta_neg","zeta_rand_neg","K_neg","M_rand_neg","L_rand_neg","S_neg")
networkAnalysis[,1:ncol(networkAnalysis)] <- sapply(networkAnalysis[,1:ncol(networkAnalysis)],as.character)
networkAnalysis[,2:ncol(networkAnalysis)] <- sapply(networkAnalysis[,2:ncol(networkAnalysis)],as.numeric)
rownames(networkAnalysis) <- 1:nrow(networkAnalysis)
write.table(networkAnalysis,paste("networkAnalysis18SV9",communityType,taxonomicLevel,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)

##To run locally.
networkAnalysis <- read.table(paste("networkAnalysis18SV9",communityType,taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Check for correlation patterns between co-occurrence networks topology and environmental parameters.
require("PerformanceAnalytics")
chart.Correlation(networkAnalysis[,c("meanLU","meanAL","meanDist","meanN","meanP","meanOrthoP","N","C","E","S","M","L","zeta")], histogram=TRUE, method="pearson")
chart.Correlation(networkAnalysis[,c("meanLU","meanAL","meanDist","meanN","meanP","meanOrthoP","N_pos","C_pos","E_pos","S_pos","M_pos","L_pos","zeta_pos")], histogram=TRUE, method="pearson")
chart.Correlation(networkAnalysis[,c("meanLU","meanAL","meanDist","meanN","meanP","meanOrthoP","N_neg","C_neg","E_neg","S_neg","M_neg","L_neg","zeta_neg")], histogram=TRUE, method="spearman")

require(Hmisc)
require(xtable)
# x is a matrix containing the data
# method : correlation method. "pearson"" or "spearman"" is supported
# removeTriangle : remove upper or lower triangle
# results :  if "html" or "latex"
# the results will be displayed in html or latex format
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex") 
  }
} 

communityType <- "BMIs"
taxonomicLevels <- c("class", "order", "family", "genus","species", "OTUID")
netMatTotal <- data.frame()
for(taxonomicLevel in taxonomicLevels){
  networkAnalysis <- read.table(paste("networkAnalysis18SV9",communityType,taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
  x <- networkAnalysis[,c("meanLU","meanAL","meanN","meanP","meanOrthoP","N","C","S","M","L","zeta")]
  netMat <- corstars(x,method="pearson",removeTriangle="upper",result="none")
  netMat <- as.data.frame(netMat[rownames(netMat) %in% c("N","C","S","M","L","zeta"),colnames(netMat) %in% c("meanLU","meanAL","meanN","meanP","meanOrthoP")])
  netMat$TaxonomicLevel <- taxonomicLevel
  netMatTotal <- rbind(netMatTotal,netMat)
}
