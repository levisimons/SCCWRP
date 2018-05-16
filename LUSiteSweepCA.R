library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)

setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage
#values.  If this file is not yet generated then proceed with the following commands
#to generate it in the first place.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Ensure that all sites have a LU_2000_5K value.
GISBioData <- subset(GISBioData, LU_2000_5K != "NA")
#Order data by LU_2000_5K.
GISBioData <- arrange(GISBioData,LU_2000_5K)

#Get number of unique LU_2000_5K values.
sitesNum <- length(unique(GISBioData$UniqueID))
#Enter number of divisions for subsampling.
divisionNum = 60
#Obtain subsampling number.
sampleNum <- as.integer(sitesNum/divisionNum)
uniqueSamples <- as.data.frame(unique(GISBioData$UniqueID))
colnames(uniqueSamples) <- c("UniqueID")

for(i in 1:divisionNum){
  lowNum=(i-1)*sampleNum+1
  highNum=i*sampleNum
  GISBioData <- arrange(GISBioData,LU_2000_5K)
  uniqueSampleSubset <- as.data.frame(uniqueSamples[lowNum:highNum,1])
  colnames(uniqueSampleSubset) <- c("UniqueID")
  GISBioDataSubset <- GISBioData[GISBioData$UniqueID %in% as.vector(uniqueSampleSubset$UniqueID),]
  #Determine the average LU_2000_5K per subsample of sites.
  meanLU_2000_5K = mean(na.omit(GISBioDataSubset$LU_2000_5K))
  #Initialize a data frame where the rows are all of the unique measurements for a given
  #subset of the data.
  #Order the data frame by measurement name.
  selected <- arrange(GISBioDataSubset,Year,UniqueID)
  eLSAInput <- as.data.frame(unique(selected$FinalID))
  colnames(eLSAInput)<-c("FinalID")
  eLSAInput <- as.data.frame(eLSAInput[order(as.character(eLSAInput$FinalID)),])
  colnames(eLSAInput)<-c("FinalID")
  
  #Add the relative taxa abundances by column to a new dataframe.
  #The rows are the unique taxa in a given subset of data.
  selected <- selected[order(selected$Year,selected$UniqueID,selected$FinalID),]
  for(ID in unique(selected$UniqueID)){
    tmp <- filter(selected, UniqueID == ID)[,c("FinalID","Measurement","UniqueID")]
    tmp <- as.data.frame(tmp[order(tmp$FinalID),])
    tmp <- tmp[-c(3)]
    colnames(tmp)<-c("FinalID",paste("Measurement",ID,sep=" "))
    eLSAInput <- join(eLSAInput,tmp,by="FinalID")
    eLSAInput$FinalID=as.character(eLSAInput$FinalID)
    eLSAInput <- eLSAInput %>% group_by(FinalID) %>% summarise_if(is.numeric,mean,na.rm=TRUE)
    #print(ID)
  }
  
  eLSAInput[is.na(eLSAInput)] <- "NA"
  
  #Determine the number of time points in the eLSA input file.
  spotNum = length(unique(selected$Year))
  #Determine the number of replicates per time point in the eLSA input file.
  #In order to ensure a uniform number of replicates per year this needs to
  #be the maximum number of replicates for all of the years available.
  repMax = 0
  for(year in unique(selected$Year)){
    tmp <- filter(selected, Year == year)[,c("UniqueID","Year")]
    repNum = length(unique(tmp$UniqueID))
    if(repNum >= repMax){repMax = repNum}
    #print (paste(repMax,repNum,year,sep=" "))
  }
  repNum = repMax
  
  #Now insert the replicates with actual data in between the "NA" dummy columns
  #which ensure that the final eLSA input file has an even number of replicates
  #per year regardless of the variations in the actual number of sites (replicates)
  #sampled per year.
  eLSAtmp <- eLSAInput[,1]
  j=1
  k=1
  nulCol <- data.frame(matrix(ncol=repNum*spotNum,nrow=length(unique(selected$FinalID))))
  nulCol[,1] <- eLSAInput[,1]
  for(year in unique(selected$Year)){
    tmp <- filter(selected, Year == year)
    rep = length(unique(tmp$UniqueID))
    for(i in 1:repNum){
      if(i <= rep){
        repLabel = paste(year,"DoneRep",i,sep="")
        j=j+1
        k=k+1
        eLSAtmp[,k] <- eLSAInput[,j]
      }
      else{
        repLabel = as.character(paste(year,"Rep",i,sep=""))
        k=k+1
        eLSAtmp[,k] <- "NA"
        #print(paste(k,repLabel,sep=" "))
      }
    }
  }
  
  eLSAInput <- eLSAtmp
  
  #Designate a unique filename.
  #N is the number of samples in the subsample group.
  #S is the number of spots, or years represented in the subsample group.
  #R is the number of replicates per year.  Many of the years will have null replicates, but a uniform number is needed for eLSA.
  #M is the mean LU_2000_5K score per subsample group.
  filename = paste("LUSweepCA",sampleNum,"Samples",lowNum,"to",highNum,"S",spotNum,"R",repNum,"M",meanLU_2000_5K,sep="")
  
  #Output file for use in eLSA.
  write.table(eLSAInput,paste(filename,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  eLSACommand = paste("lsa_compute ",filename,".txt ","-r ",repNum," -s ",spotNum," ",filename,"Network.txt;",sep="")
  print(eLSACommand)
}

#Read in eLSA output.
#Compute network statistics of the likeliest association networks between taxa.
library(igraph)
library(network)
library(stringr)
#Read in site data.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Ensure that all sites have a CSCI value.
GISBioData <- subset(GISBioData, LU_2000_5K != "NA")
networkfiles <- Sys.glob("LUSweepCA83Samples*Network.txt")
networkAnalysis <- data.frame()
networkConTaxa <- data.frame()
networkCovTaxa <- data.frame()
#Define a 'not in' function.
'%!in%' <- function(x,y)!('%in%'(x,y))
for(networkFile in networkfiles){
  networkdata <- read.delim(networkFile,header=TRUE, sep="\t",as.is=T,check.names=FALSE)
  #Filter out association network data based on P and Q scores, for the local similarity
  #between two factors, with values less than 0.01.
  networkdata <- filter(networkdata, P <= 1e-4)
  networkdata <- filter(networkdata, Q <= 1e-4)
  names(networkdata)[names(networkdata)=="LS"]<-"weight"
  #networkdata <- filter(networkdata,abs(weight)>=0.3)
  meanLU <- as.numeric(str_match(networkFile,"M(.*?)Network")[2])
  #Generate network graph and begin calculating network parameters.
  networkgraph=graph.data.frame(networkdata,directed=FALSE)
  Network_size<-network.size(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
  if(ecount(networkgraph)>0){
    #Get the full weighted adjacency matrix.
    networkmatrix <- as.matrix(get.adjacency(networkgraph,attr='weight'))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_network <- eigen(networkmatrix)
    #Get the real component first eigenvalue.
    lambda_network_m <- Re(lambda_network$values[1])
    #Generate randomized version of full weighted adjacency matrix.
    set.seed(1)
    randnetworkmatrix <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_rand <- eigen(randnetworkmatrix)
    #Get the real component of the first eigenvalue.
    lambda_rand_m <- Re(lambda_rand$values[1])
    #Calculate stability parameter.
    gamma <- lambda_network_m/lambda_rand_m
    #Calculate the degree heterogeneity.
    networkmatrix[upper.tri(networkmatrix)] <- 0
    networkmatrix <- ifelse(networkmatrix!=0,1,networkmatrix)
    zeta <- mean(colSums(networkmatrix)^2)/mean(colSums(networkmatrix))^2
  }
  #Filter contravariant network data based on local similarity scores.
  networkdataCon <- subset(networkdata,networkdata$weight<0)
  #Aggregate significantly contravarying taxa.
  networkdataConTemp <- networkdataCon[,c("X","Y","weight")]
  networkdataConTemp <- as.data.frame(table(append(networkdataConTemp$X,networkdataConTemp$Y,after=length(networkdataConTemp$X))))
  networkdataConTemp$meanLU <- meanLU
  networkConTaxa <- rbind(networkConTaxa,networkdataConTemp)
  #Generate network graph and begin calculating network parameters.
  networkgraphCon=graph.data.frame(networkdataCon,directed=FALSE)
  if(ecount(networkgraphCon)>0){
    #Get the full weighted adjacency matrix.
    networkmatrix <- as.matrix(get.adjacency(networkgraphCon,attr='weight'))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_network <- eigen(networkmatrix)
    #Get the real component first eigenvalue.
    lambda_network_m_Con <- Re(lambda_network$values[1])
    #Generate randomized version of full weighted adjacency matrix.
    set.seed(1)
    randnetworkmatrix <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_rand_Con <- eigen(randnetworkmatrix)
    #Get the real component of the first eigenvalue.
    lambda_rand_Con <- Re(lambda_rand_Con$values[1])
    #Calculate stability parameter.
    gamma_Con <- lambda_network_m_Con/lambda_rand_Con
    #Calculate the degree heterogeneity.
    networkmatrixCon <- networkmatrix
    networkmatrixCon[upper.tri(networkmatrixCon)] <- 0
    networkmatrixCon <- ifelse(networkmatrixCon!=0,1,networkmatrixCon)
    zeta_Con <- mean(colSums(networkmatrixCon)^2)/mean(colSums(networkmatrixCon))^2
    # Generate adjacency matrix of relative taxa abundance correlations
    adj= as.network(get.adjacency(networkgraphCon,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
    # Get the number of unique network edges
    networkEdgecount <- network.edgecount(adj)
    networkEdgecountCon <- networkEdgecount
    # Get the number of nodes
    networkNodecount <- network.size(adj)
    # Get the average degree per node.
    k <- (2*networkEdgecount)/networkNodecount
    # Get the random characteristic path length.
    networkRandLength <- 0.5+((log(networkNodecount)-0.5772156649)/log(k))
    # Get the random clustering coefficient.
    networkRandClustering <- k/networkNodecount
    # Get the network density.
    networkDensity <- network.density(adj)
    con_C <- networkDensity
    # Calculate the modularity of the network.
    networkModularity <- modularity(cluster_edge_betweenness(networkgraphCon, weights=NULL,directed=FALSE))
    con_M <- networkModularity
    # Calculate the number of groups related to the modularity value.
    networkModGroups <- length(cluster_edge_betweenness(networkgraphCon, weights=NULL,directed=FALSE))
    # Calculate the average network path length
    networkLength <- mean_distance(networkgraphCon,directed=FALSE)
    con_L <- networkLength
    # Calculate the clustering coefficient
    networkClustering <- transitivity(networkgraphCon,type="globalundirected",isolate="zero")
    con_Cl <- networkClustering
    # Calcuate the log ratio of clustering coefficients.
    l_con_rCl <- log(networkClustering/networkRandClustering)
    # Calculate the modularity of the random network.
    networkRandModularity <- (1-(2/sqrt(networkNodecount)))*(2/k)^(2/3)
    # Calculate the log ratio of the modularities.
    l_con_rM <- log(networkModularity/networkRandModularity)
    # Get log ratio of characteristic path lengths.
    l_con_rL <- log(networkLength/networkRandLength)
  }
  #Filter covariant network data based on local similarity scores.
  networkdataCov <- subset(networkdata,networkdata$weight>0)
  #Aggregate significantly contravarying taxa.
  networkdataCovTemp <- networkdataCov[,c("X","Y","weight")]
  networkdataCovTemp <- as.data.frame(table(append(networkdataCovTemp$X,networkdataCovTemp$Y,after=length(networkdataCovTemp$X))))
  networkdataCovTemp$meanLU <- meanLU
  networkCovTaxa <- rbind(networkCovTaxa,networkdataCovTemp)
  #Generate network graph and begin calculating network parameters.
  networkgraphCov=graph.data.frame(networkdataCov,directed=FALSE)
  if(ecount(networkgraph)>0){
    #Get the full weighted adjacency matrix.
    networkmatrix <- as.matrix(get.adjacency(networkgraphCov,attr='weight'))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_network <- eigen(networkmatrix)
    #Get the real component first eigenvalue.
    lambda_network_m_Cov <- Re(lambda_network$values[1])
    #Generate randomized version of full weighted adjacency matrix.
    set.seed(1)
    randnetworkmatrix <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_rand_Cov <- eigen(randnetworkmatrix)
    #Get the real component of the first eigenvalue.
    lambda_rand_Cov <- Re(lambda_rand_Cov$values[1])
    #Calculate stability parameter.
    gamma_Cov <- lambda_network_m_Cov/lambda_rand_Cov
    #Calculate the degree heterogeneity.
    networkmatrixCov <- networkmatrix
    networkmatrixCov[upper.tri(networkmatrixCov)] <- 0
    networkmatrixCov <- ifelse(networkmatrixCov!=0,1,networkmatrixCov)
    zeta_Cov <- mean(colSums(networkmatrixCov)^2)/mean(colSums(networkmatrixCov))^2
    # Generate adjacency matrix of relative taxa abundance correlations
    adj= as.network(get.adjacency(networkgraphCov,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
    # Get the number of unique network edges
    networkEdgecount <- network.edgecount(adj)
    networkEdgecountCov <- networkEdgecount
    # Get the number of nodes
    networkNodecount <- network.size(adj)
    # Get the average degree per node.
    k <- (2*networkEdgecount)/networkNodecount
    # Get the random characteristic path length.
    networkRandLength <- 0.5+((log(networkNodecount)-0.5772156649)/log(k))
    # Get the random clustering coefficient.
    networkRandClustering <- k/networkNodecount
    # Get the network density.
    networkDensity <- network.density(adj)
    cov_C <- networkDensity
    # Calculate the modularity of the network.
    networkModularity <- modularity(cluster_edge_betweenness(networkgraphCov, weights=NULL,directed=FALSE))
    cov_M <- networkModularity
    # Calculate the number of groups related to the modularity value.
    networkModGroups <- length(cluster_edge_betweenness(networkgraphCov, weights=NULL,directed=FALSE))
    # Calculate the average network path length
    networkLength <- mean_distance(networkgraphCov,directed=FALSE)
    cov_L <- networkLength
    # Calculate the clustering coefficient
    networkClustering <- transitivity(networkgraphCov,type="globalundirected",isolate="zero")
    cov_Cl <- networkClustering
    # Calcuate the log ratio of clustering coefficients.
    l_cov_rCl <- log(networkClustering/networkRandClustering)
    # Calculate the modularity of the random network.
    networkRandModularity <- (1-(2/sqrt(networkNodecount)))*(2/k)^(2/3)
    # Calculate the log ratio of the modularities.
    l_cov_rM <- log(networkModularity/networkRandModularity)
    # Get log ratio of characteristic path lengths.
    l_cov_rL <- log(networkLength/networkRandLength)
  }
  dat <- data.frame()
  dat[1,1] <- networkFile
  dat[1,2] <- meanLU
  dat[1,3] <- l_con_rL
  dat[1,4] <- l_con_rCl
  dat[1,5] <- l_con_rM
  dat[1,6] <- l_cov_rL
  dat[1,7] <- l_cov_rCl
  dat[1,8] <- l_cov_rM
  dat[1,9] <- lambda_network_m
  dat[1,10] <- con_L
  dat[1,11] <- con_Cl
  dat[1,12] <- con_M
  dat[1,13] <- cov_L
  dat[1,14] <- cov_Cl
  dat[1,15] <- cov_M
  dat[1,16] <- zeta
  dat[1,17] <- con_C
  dat[1,18] <- cov_C
  dat[1,19] <- Network_size
  dat[1,20] <- Pm <- networkEdgecountCov/(networkEdgecountCov+networkEdgecountCon)
  dat[1,21] <- lambda_network_m_Con
  dat[1,22] <- lambda_network_m_Cov
  dat[1,23] <- zeta_Con
  dat[1,24] <- zeta_Cov
  dat[1,25] <- gamma_Con
  dat[1,26] <- gamma_Cov
  dat[1,27] <- gamma
  dat[1,28] <- lambda_rand_m
  dat[1,29] <- lambda_rand_Con
  dat[1,30] <- lambda_rand_Cov
  networkAnalysis <- rbind(networkAnalysis,dat)
  print(paste(networkFile,meanLU,l_con_rL,l_con_rCl,l_con_rM,l_cov_rL,l_cov_rCl,l_cov_rM,lambda_network_m,con_L,con_Cl,con_M,cov_L,cov_Cl,cov_M,zeta,con_C,cov_C,Network_size,Pm,lambda_network_m_Con,lambda_network_m_Cov,zeta_Con,zeta_Cov,gamma_Con,gamma_Cov,gamma,lambda_rand_m,lambda_rand_Con,lambda_rand_Cov))
}
colnames(networkAnalysis) <- c("filename","meanLU","l_con_rL","l_con_rCl","l_con_rM","l_cov_rL","l_cov_rCl","l_cov_rM","lambda_network_m","con_L","con_Cl","con_M","cov_L","cov_Cl","cov_M","zeta","con_C","cov_C","Network_size","Pm","lambda_network_m_Con","lambda_network_m_Cov","zeta_Con","zeta_Cov","gamma_Con","gamma_Cov","gamma","lambda_rand_m","lambda_rand_Con","lambda_rand_Cov")
networkAnalysis[networkAnalysis=="-Inf"] <- NA
networkAnalysis[networkAnalysis=="Inf"] <- NA
networkAnalysis <- arrange(networkAnalysis,meanLU)
#write.table(networkAnalysis,"CSCISiteSweepCA.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Regression between network parameters.
library(Hmisc)
res2<-rcorr(as.matrix(networkAnalysis[,-c(1)]),type="spearman")
library(corrplot)
#Correlation plot.  Insignificant correlations are leaved blank.
corrplot(res2$r, type="upper", order="hclust", p.mat = res2$P, sig.level = 0.01, insig = "blank")
library("PerformanceAnalytics")
#Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
chart.Correlation(networkAnalysis[,-c(1)], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("zeta","lambda_network_m")], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("l_cov_rM","lambda_network_m_Cov")], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("l_con_rM","lambda_network_m_Con")], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("cov_C","lambda_network_m_Cov")], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("con_C","lambda_network_m_Con")], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("Pm","lambda_network_m")], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","l_con_rL")], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","l_cov_rL")], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","l_con_rM")], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","l_cov_rM")], histogram=TRUE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta","l_cov_rM","l_con_rM")], histogram=TRUE, method="spearman")

#To generate map of data for a given environmental parameter in California.
library(ggmap)
library(maps)
library(mapdata)
dev.off()
MapCoordinates <- data.frame(GISBioData$LU_2000_5K,GISBioData$Longitude,GISBioData$Latitude)
colnames(MapCoordinates) = c('LU','lon','lat')
MapCoordinates <- na.omit(MapCoordinates)
mapBoundaries <- make_bbox(lon=MapCoordinates$lon,lat=MapCoordinates$lat,f=0.1)
CalMap <- get_map(location=mapBoundaries,maptype="satellite",source="google")
CalMap <- ggmap(CalMap)+geom_point(data = MapCoordinates, mapping = aes(x = lon, y = lat, color = LU))
CalMap
