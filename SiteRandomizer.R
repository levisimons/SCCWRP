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
siteNum=100
for(i in 1:siteNum){
  #Read in site data containing biological counts, water chemistry, and land usage
  #values.  If this file is not yet generated then proceed with the following commands
  #to generate it in the first place.
  GISBiochemData <- read.table("GISBiochemData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
  #Ensure that all sites have a CSCI value.
  GISBiochemData <- subset(GISBiochemData, CSCI != "NA")
  
  #Randomly subsample the site data by individual sample location-date (UniqueID).
  SampleSize=50
  subsample <- unique(GISBiochemData$UniqueID)
  subsample <- sample(subsample,SampleSize)
  
  #Filter data by a subsample of UniqueID.
  GISBiochemData <- GISBiochemData[GISBiochemData$UniqueID %in% subsample,]
  
  #Determine the average CSCI per subsample of sites.
  meanCSCI = mean(na.omit(GISBiochemData$CSCI))
  
  #Initialize a data frame where the rows are all of the unique measurements for a given
  #subset of the data.
  #Order the data frame by measurement name.
  selected <- arrange(GISBiochemData,Year,UniqueID)
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
    tmp <- filter(selected, Year == year)[,c(6,7)]
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
  #M is the mean CSCI score per subsample group.
  filename = paste("N",SampleSize,"S",spotNum,"R",repNum,"M",meanCSCI,sep="")
  
  #Output file for use in eLSA.
  write.table(eLSAInput,paste(filename,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  print(filename)
}


#Read in eLSA output.
#Compute network statistics of the likeliest association networks between taxa.
library(igraph)
library(network)
library(stringr)
#Read in site data.
GISBiochemData <- read.table("GISBiochemData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Ensure that all sites have a CSCI value.
GISBiochemData <- subset(GISBiochemData, CSCI != "NA")
#Get unique identifiers for algal, invertebrate, and chemical measurement types.
algae <- subset(GISBiochemData,GISBiochemData$MeasurementType=="Algal relative abundance")
algaeID <- unique(algae$FinalID)
insect <-subset(GISBiochemData,GISBiochemData$MeasurementType=="Invertebrate relative abundances" | GISBiochemData$MeasurementType=="Invertebrate relative abundance")
insectID <- unique(insect$FinalID)
chem <- subset(GISBiochemData,GISBiochemData$MeasurementType!="Algal relative abundance" & GISBiochemData$MeasurementType!="Invertebrate relative abundances" & GISBiochemData$MeasurementType!="Invertebrate relative abundance")
chemID <- unique(chem$FinalID)
networkfiles <- Sys.glob("N50*Network.txt")
networkAnalysis <- data.frame()
#Define a 'not in' function.
'%!in%' <- function(x,y)!('%in%'(x,y))
for(networkFile in networkfiles){
  networkdata <- read.delim(networkFile,header=TRUE, sep="\t",as.is=T,check.names=FALSE)
  #Filter out association network data based on P scores, for the local similarity
  #between two factors, with values less than 0.05.
  networkdata <- filter(networkdata, P <= 0.01)
  names(networkdata)[names(networkdata)=="LS"]<-"weight"
  meanCSCI <- as.numeric(str_match(networkFile,"M(.*?)Network")[2])
  #Remove some subset of chemical and biological factors as nodes from the network.
  networkdata <- subset(networkdata,networkdata$X %!in% chemID)
  networkdata <- subset(networkdata,networkdata$Y %!in% chemID)
  #Generate network graph and begin calculating network parameters.
  networkgraph=graph.data.frame(networkdata,directed=FALSE)
  if(ecount(networkgraph)>0){
    #Get the full weighted adjacency matrix.
    networkmatrix <- as.matrix(get.adjacency(networkgraph,attr='weight'))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_network <- eigen(networkmatrix)
    #Get the real component first eigenvalue.
    lambda_network_1 <- Re(lambda_network$values[1])
    #Generate randomized version of full weighted adjacency matrix.
    set.seed(1)
    randnetworkmatrix <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_rand <- eigen(randnetworkmatrix)
    #Get the real component of the first eigenvalue.
    lambda_rand_1 <- Re(lambda_rand$values[1])
    #Calculate stability parameter.
    gamma <- lambda_network_1/lambda_rand_1
  }
  #Filter contravariant network data based on local similarity scores.
  networkdataCon <- subset(networkdata,networkdata$weight<0)
  #Generate network graph and begin calculating network parameters.
  networkgraphCon=graph.data.frame(networkdataCon,directed=FALSE)
  if(ecount(networkgraphCon)>0){
    # Generate adjacency matrix of relative taxa abundance correlations
    adj= as.network(get.adjacency(networkgraphCon,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
    # Get the number of unique network edges
    networkEdgecount <- network.edgecount(adj)
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
  #Generate network graph and begin calculating network parameters.
  networkgraphCov=graph.data.frame(networkdataCov,directed=FALSE)
  if(ecount(networkgraph)>0){
    # Generate adjacency matrix of relative taxa abundance correlations
    adj= as.network(get.adjacency(networkgraphCov,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
    # Get the number of unique network edges
    networkEdgecount <- network.edgecount(adj)
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
  dat[1,2] <- meanCSCI
  dat[1,3] <- l_con_rL
  dat[1,4] <- l_con_rCl
  dat[1,5] <- l_con_rM
  dat[1,6] <- l_cov_rL
  dat[1,7] <- l_cov_rCl
  dat[1,8] <- l_cov_rM
  dat[1,9] <- gamma
  dat[1,10] <- con_L
  dat[1,11] <- con_Cl
  dat[1,12] <- con_M
  dat[1,13] <- cov_L
  dat[1,14] <- cov_Cl
  dat[1,15] <- cov_M
  networkAnalysis <- rbind(networkAnalysis,dat)
  print(paste(networkFile,meanCSCI,l_con_rL,l_con_rCl,l_con_rM,l_cov_rL,l_cov_rCl,l_cov_rM,gamma,con_L,con_Cl,con_M,cov_L,cov_Cl,cov_M))
}
colnames(networkAnalysis) <- c("filename","meanCSCI","l_con_rL","l_con_rCl","l_con_rM","l_cov_rL","l_cov_rCl","l_cov_rM","gamma","con_L","con_Cl","con_M","cov_L","cov_Cl","cov_M")
networkAnalysis[networkAnalysis=="-Inf"] <- NA
networkAnalysis[networkAnalysis=="Inf"] <- NA

#Logistic regression between network parameters and CSCI
library(PerformanceAnalytics)
library(aod)
library(glmm)
library(rcompanion)
model.vars <- names(networkAnalysis)[-c(1,2)]
model.list <- lapply(model.vars, function(x){
  glm(substitute(meanCSCI ~ i, list(i=as.name(x))),data=networkAnalysis)
})
lapply(model.list,summary)
lapply(model.list,nagelkerke)
