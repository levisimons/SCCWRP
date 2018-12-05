require(plyr)
require(dplyr)
require(ape)
require(vegan)
require(data.table)
require(igraph)
require(network)
require(stringr)

#Read in eLSA output.
#Compute network statistics of the likeliest association networks between taxa.
setwd("~/Desktop/SCCWRP")
#setwd("~/panfs/SCCWRP")

#Read in site data.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Ensure that all sites have a land use value.
GISBioData <- subset(GISBioData, LU_2000_5K != "NA")
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
networkfiles <- Sys.glob("LUSweepCAWatershedPerm100*Years*Reps*MeanLU*nTaxa*Network.txt")
networkAnalysis <- data.frame()
networkConTaxa <- data.frame()
networkCovTaxa <- data.frame()

for(networkFile in networkfiles){
  networkdata <- read.delim(networkFile,header=TRUE, sep="\t",as.is=T,check.names=FALSE)
  #Filter out association network data based on P and Q scores, for the local similarity
  #between two factors, with values less than a particuar threshold.
  networkdata <- filter(networkdata, P <= 1e-4)
  networkdata <- filter(networkdata, Q <= 1e-4)
  names(networkdata)[names(networkdata)=="LS"]<-"weight"
  meanLU <- as.numeric(str_match(networkFile, "LUSweepCAWatershedPerm100(.*?)Years(.*?)Reps(.*?)MeanLU(.*?)nTaxa(.*?)Network.txt")[5])
  #Generate network graph and begin calculating network parameters.
  if(nrow(networkdata) > 0){
    networkgraph=graph.data.frame(networkdata,directed=FALSE)
    Network_size<-network.size(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
    if(ecount(networkgraph)>0){
      #Get the full weighted adjacency matrix.
      networkmatrix <- as.matrix(get.adjacency(networkgraph,attr='weight'))
      #Mean interaction strength
      meanStrength <- mean(abs(networkmatrix))
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
      #Calculate modularity
      networkModularity <- modularity(cluster_edge_betweenness(networkgraph, weights=NULL,directed=FALSE))
      M <- networkModularity
      networkNodecount <-network.size(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
      # Get the number of unique network edges
      networkEdgecount <- network.edgecount(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
      # Get the number of nodes
      networkNodecount <- network.size(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
      # Get the average degree per node.
      k <- (2*networkEdgecount)/networkNodecount
      # Calculate the modularity of the random network.
      networkRandModularity <- (1-(2/sqrt(networkNodecount)))*(2/k)^(2/3)
      # Calculate the log ratio of the modularities.
      l_rM <- log(networkModularity/networkRandModularity)
    }
    #Filter contravariant network data based on local similarity scores.
    networkdataCon <- subset(networkdata,networkdata$weight<0)
    if(nrow(networkdataCon)>0){
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
        #Mean interaction strength
        meanStrength_Con <- mean(abs(networkmatrix))
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
        #Calculate the degree heterogeneity of the corresponding random network.
        randnetworkmatrixCon <- randnetworkmatrix
        randnetworkmatrixCon[upper.tri(randnetworkmatrixCon)] <- 0
        randnetworkmatrixCon <- ifelse(randnetworkmatrixCon!=0,1,randnetworkmatrixCon)
        zeta_rand_Con <- mean(colSums(randnetworkmatrixCon)^2)/mean(colSums(randnetworkmatrixCon))^2
        # Log response ratio of degree heterogeneity.
        l_con_rzeta <- log(zeta_Con/zeta_rand_Con)
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
    }
    
    #Filter covariant network data based on local similarity scores.
    networkdataCov <- subset(networkdata,networkdata$weight>0)
    if(nrow(networkdataCov)>0){
      #Aggregate significantly contravarying taxa.
      networkdataCovTemp <- networkdataCov[,c("X","Y","weight")]
      networkdataCovTemp <- as.data.frame(table(append(networkdataCovTemp$X,networkdataCovTemp$Y,after=length(networkdataCovTemp$X))))
      networkdataCovTemp$meanLU <- meanLU
      networkCovTaxa <- rbind(networkCovTaxa,networkdataCovTemp)
      #Generate network graph and begin calculating network parameters.
      networkgraphCov=graph.data.frame(networkdataCov,directed=FALSE)
      if(ecount(networkgraphCov)>0){
        #Get the full weighted adjacency matrix.
        networkmatrix <- as.matrix(get.adjacency(networkgraphCov,attr='weight'))
        #Mean interaction strength
        meanStrength_Cov <- mean(abs(networkmatrix))
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
        #Calculate the degree heterogeneity of the corresponding random network.
        randnetworkmatrixCov <- randnetworkmatrix
        randnetworkmatrixCov[upper.tri(randnetworkmatrixCov)] <- 0
        randnetworkmatrixCov <- ifelse(randnetworkmatrixCov!=0,1,randnetworkmatrixCov)
        zeta_rand_Cov <- mean(colSums(randnetworkmatrixCov)^2)/mean(colSums(randnetworkmatrixCov))^2
        # Log response ratio of degree heterogeneity.
        l_cov_rzeta <- log(zeta_Cov/zeta_rand_Cov)
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
    dat[1,31] <- M
    dat[1,32] <- l_rM
    dat[1,33] <- meanStrength
    dat[1,34] <- meanStrength_Cov
    dat[1,35] <- meanStrength_Con
    dat[1,36] <- zeta_rand_Con
    dat[1,37] <- l_con_rzeta
    dat[1,38] <- zeta_rand_Cov
    dat[1,39] <- l_cov_rzeta
    networkAnalysis <- rbind(networkAnalysis,dat)
    print(paste(networkFile,meanLU,l_con_rL,l_con_rCl,l_con_rM,l_cov_rL,l_cov_rCl,l_cov_rM,lambda_network_m,con_L,con_Cl,con_M,cov_L,cov_Cl,cov_M,zeta,con_C,cov_C,Network_size,Pm,lambda_network_m_Con,lambda_network_m_Cov,zeta_Con,zeta_Cov,gamma_Con,gamma_Cov,gamma,lambda_rand_m,lambda_rand_Con,lambda_rand_Cov,M,l_rM,meanStrength,meanStrength_Cov,meanStrength_Con,zeta_rand_Con,l_con_rzeta,zeta_rand_Cov,l_cov_rzeta))
  }
}
colnames(networkAnalysis) <- c("filename","meanLU","l_con_rL","l_con_rCl","l_con_rM","l_cov_rL","l_cov_rCl","l_cov_rM","lambda_network_m","con_L","con_Cl","con_M","cov_L","cov_Cl","cov_M","zeta","con_C","cov_C","Network_size","Pm","lambda_network_m_Con","lambda_network_m_Cov","zeta_Con","zeta_Cov","gamma_Con","gamma_Cov","gamma","lambda_rand_m","lambda_rand_Con","lambda_rand_Cov","M","l_rM","meanStrength","meanStrength_Cov","meanStrength_Con","zeta_rand_Con","l_con_rzeta","zeta_rand_Cov","l_cov_rzeta")
networkAnalysis[networkAnalysis=="-Inf"] <- NA
networkAnalysis[networkAnalysis=="Inf"] <- NA
networkAnalysis <- arrange(networkAnalysis,meanLU)
write.table(networkAnalysis,"LUSweepCAWatershedPerm100Networks.txt",quote=FALSE,sep="\t",row.names = FALSE)
#Send off output file to be analyzed by the script NetworkTopologyAnalysis.R
