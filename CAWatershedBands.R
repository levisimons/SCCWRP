library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)

#This script looks at co-occurrence networks generated on a HUC-8 watershed scale from the SCCWRP
#California data set.  Each watershed is divided into land use bands, and then the core networks
#which persist across all land use bands by watershed are compared to infer the relationship between
#co-occurrence network structure and stress.
setwd("~/Desktop/SCCWRP")

#Read in site data containing biological counts, water chemistry, and land usage values.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#How many samples per group for generating a co-occurrence network?
groupNum=10
#Run through analysis on SCCWRP archive on a watershed-level scale.
#Select watersheds with a large enough set of samples for analysis.
watersheds <- subset(as.data.frame(table(SCCWRP$Watershed)),Freq>=3*groupNum)
colnames(watersheds) <- c("Watershed","Samples")

#Watersheds with a large enough set of samples to generate co-occurrence networks
#for each land use band.
LargeWS <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(LargeWS) <- c("Watershed","Samples")

#This portion of code determines which watersheds have enough samples in all land use bands
#to generate networks.
#Set land use thresholds.
PD <- 5 #Partially degraded land use threshold.
HD <- 15 #Highly degraded land use threshold.
#Subset data by watershed, and then by land use bands.
for(WS in watersheds$Watershed){
  #Get samples per watershed.
  GISBioDataSubset <- subset(GISBioData,Watershed==WS)
  j=0
  for(i in 1:3){
    if(i==1){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K<=PD)
      if(length(unique(GISBioDataSubsetLU$UniqueID))>=groupNum){
        j=j+1
      }
    }
    if(i==2){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K > PD & LU_2000_5K<=HD)
      if(length(unique(GISBioDataSubsetLU$UniqueID))>=groupNum){
        j=j+1
      }
    }
    if(i==3){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K > HD)
      if(length(unique(GISBioDataSubsetLU$UniqueID))>=groupNum){
        j=j+1
      }
    }
    #Get watersheds with a minimum number of samples per watershed and land use band.
    if(j==3){
      LargeWS <- rbind(LargeWS,subset(watersheds,watersheds$Watershed==WS))
      }
  }
}

#Generate LSA input files for each watershed which contains a sufficient number of samples
#for each land use band.
for(WS in LargeWS$Watershed){
  GISBioDataSubset <- subset(GISBioData,Watershed==WS)
  for(i in 1:3){
    if(i==1){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K<=PD)
    }
    if(i==2){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K > PD & LU_2000_5K<=HD)
    }
    if(i==3){
      GISBioDataSubsetLU <- subset(GISBioDataSubset,LU_2000_5K > HD)
    }
    selected <- GISBioDataSubsetLU
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
    #S is the number of spots, or years represented in the subsample group.
    #R is the number of replicates per year.  Many of the years will have null replicates, but a uniform number is needed for eLSA.
    #M is the mean LU_2000_5K score per subsample group.
    sampleNum <- length(unique(selected$UniqueID))
    meanLU = mean(na.omit(selected$LU_2000_5K))
    filename <- paste("Watershed",WS,"NSamples",sampleNum,"S",spotNum,"R",repNum,"M",meanLU,sep="")
    #Output file for use in eLSA.
    write.table(eLSAInput,paste(filename,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    eLSACommand = paste("lsa_compute ",filename,".txt ","-r ",repNum," -s ",spotNum," ",filename,"Network.txt;",sep="")
    print(eLSACommand)
  }
}


#Read in eLSA output.
#Compute network statistics of the likeliest association networks between taxa.
library(igraph)
library(network)
library(stringr)
#Read in site data.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Ensure that all sites have a land use value.
GISBioData <- subset(GISBioData, LU_2000_5K != "NA")
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Get list of taxa of interest.
taxaList <- as.data.frame(table(GISBioData$FinalID))
colnames(taxaList) <- c("FinalID","Freq")
#Get network files to analyze.
networkfiles <- as.data.frame(Sys.glob("Watershed*NSamples*Network.txt"))
colnames(networkfiles) <- c("filename")
#Add watersheds to the filenames
networkfiles$Watershed <- str_split_fixed(str_split_fixed(networkfiles$filename,"Watershed",2)[,2],"NSamples",2)[,1]
networkfiles <- arrange(networkfiles,Watershed)

#Analyze network topologies for co-occurrence networks generated from taxa which occur
#across all land use bands on a per watershed basis.
networkAnalysis <- data.frame()
for(WS in unique(networkfiles$Watershed)){
  WSSubset <- subset(networkfiles,Watershed==WS)
  coreTaxa <- list()
  #Get the taxa which are common to a watershed across all land use bands.
  for(networkFile in unique(WSSubset$filename)){
    networkdata <- read.delim(networkFile,header=TRUE, sep="\t",as.is=T,check.names=FALSE)
    #Filter out taxa not being focused on in a particular study.
    networkdata <- networkdata[networkdata$X %in% as.vector(taxaList$FinalID) | networkdata$Y %in% as.vector(taxaList$FinalID),]
    #Filter out association network data based on P and Q scores, for the local similarity
    #between two factors, with values less than a particuar threshold.
    networkdata <- filter(networkdata, P <= 1e-2)
    networkdata <- filter(networkdata, Q <= 1e-2)
    names(networkdata)[names(networkdata)=="LS"]<-"weight"
    #Get taxa present for each co-occurrence network.
    uniqueTaxa <- append(unique(networkdata$X),unique(networkdata$Y))
    if(length(coreTaxa)==0){
      coreTaxa <- uniqueTaxa
    }
    if(length(coreTaxa)>0){
      coreTaxa <- Reduce(intersect,list(coreTaxa,uniqueTaxa))
    }
    #print(paste(WS,length(coreTaxa)))
  }
  for(networkFile in unique(WSSubset$filename)){
    networkdata <- read.delim(networkFile,header=TRUE, sep="\t",as.is=T,check.names=FALSE)
    #Focus only on taxa common to all land use bands in a watershed.
    networkdata <- networkdata[networkdata$X %in% as.vector(coreTaxa) | networkdata$Y %in% as.vector(coreTaxa),]
    #Filter out association network data based on P and Q scores, for the local similarity
    #between two factors, with values less than a particuar threshold.
    networkdata <- filter(networkdata, P <= 1e-2)
    networkdata <- filter(networkdata, Q <= 1e-2)
    names(networkdata)[names(networkdata)=="LS"]<-"weight"
    #Get mean land use sample group.
    options(digits=15)
    meanLU <- as.numeric(str_match(networkFile,"M(.*?)Network")[2])
    #Generate network graph and begin calculating network parameters.
    networkgraph=graph.data.frame(networkdata,directed=FALSE)
    #Covariant co-occurrence network
    networkgraphCov <- graph.data.frame(subset(networkdata,weight>0),directed=FALSE)
    if(ecount(networkgraphCov)>0){
      # Generate adjacency matrix of relative taxa abundance correlations
      adj= as.network(get.adjacency(networkgraphCov,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
      cov_C = network.density(adj)
      #Calculate the degree heterogeneity.
      networkmatrix <- as.matrix(get.adjacency(networkgraphCov,attr='weight'))
      #Mean interaction strength
      cov_meanStrength <- mean(abs(networkmatrix))
      networkmatrix[upper.tri(networkmatrix)] <- 0
      networkmatrix <- ifelse(networkmatrix!=0,1,networkmatrix)
      cov_zeta <- mean(colSums(networkmatrix)^2)/mean(colSums(networkmatrix))^2
      #Generate randomized version of full weighted adjacency matrix.
      #Calculate the degree heterogeneity of the corresponding random network.
      set.seed(1)
      randnetworkmatrixCov <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
      randnetworkmatrixCov[upper.tri(randnetworkmatrixCov)] <- 0
      randnetworkmatrixCov <- ifelse(randnetworkmatrixCov!=0,1,randnetworkmatrixCov)
      cov_rand_zeta <- mean(colSums(randnetworkmatrixCov)^2)/mean(colSums(randnetworkmatrixCov))^2
      # Log response ratio of degree heterogeneity.
      cov_lr_zeta <- log(cov_zeta/cov_rand_zeta)
      #Calculate modularity
      cov_M <- modularity(cluster_edge_betweenness(networkgraphCov, weights=NULL,directed=FALSE))
      #Edge counts
      cov_Edge <- ecount(networkgraphCov)
      #Node counts
      cov_nodes <- vcount(networkgraphCov)
      #Get the average degree per node.
      cov_k <- (2*cov_Edge)/cov_nodes
      #Characteristic path length.
      cov_L <- mean_distance(networkgraphCov,directed=FALSE)
      #Log ratio of characteristic path length to its random counterpart.
      cov_lr_L <- log(mean_distance(networkgraphCov,directed=FALSE)/(0.5+((log(cov_nodes)-0.5772156649)/log(cov_k))))
      #Clustering coefficient.
      cov_Cl <- transitivity(networkgraphCon,type="globalundirected",isolate="zero")
      #Log ratio of clustering coefficient to its random counterpart.
      cov_lr_Cl <- log(transitivity(networkgraphCov,type="globalundirected",isolate="zero")/(cov_k/cov_nodes))
    }
    #Contravariant co-occurrence network
    networkgraphCon <- graph.data.frame(subset(networkdata,weight<0),directed=FALSE)
    if(ecount(networkgraphCon)>0){
      # Generate adjacency matrix of relative taxa abundance correlations
      adj= as.network(get.adjacency(networkgraphCon,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
      con_C = network.density(adj)
      #Calculate the degree heterogeneity.
      networkmatrix <- as.matrix(get.adjacency(networkgraphCon,attr='weight'))
      #Mean interaction strength
      con_meanStrength <- mean(abs(networkmatrix))
      networkmatrix[upper.tri(networkmatrix)] <- 0
      networkmatrix <- ifelse(networkmatrix!=0,1,networkmatrix)
      con_zeta <- mean(colSums(networkmatrix)^2)/mean(colSums(networkmatrix))^2
      #Generate randomized version of full weighted adjacency matrix.
      #Calculate the degree heterogeneity of the corresponding random network.
      set.seed(1)
      randnetworkmatrixCon <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
      randnetworkmatrixCon[upper.tri(randnetworkmatrixCon)] <- 0
      randnetworkmatrixCon <- ifelse(randnetworkmatrixCon!=0,1,randnetworkmatrixCon)
      con_rand_zeta <- mean(colSums(randnetworkmatrixCon)^2)/mean(colSums(randnetworkmatrixCon))^2
      # Log response ratio of degree heterogeneity.
      con_lr_zeta <- log(con_zeta/con_rand_zeta)
      #Calculate modularity
      con_M <- modularity(cluster_edge_betweenness(networkgraphCon, weights=NULL,directed=FALSE))
      #Edge counts
      con_Edge <- ecount(networkgraphCon)
      #Node counts
      con_nodes <- vcount(networkgraphCon)
      #Get the average degree per node.
      con_k <- (2*con_Edge)/con_nodes
      #Characteristic path length.
      con_L <- mean_distance(networkgraphCon,directed=FALSE)
      #Log ratio of characteristic path length to its random counterpart.
      con_lr_L <- log(mean_distance(networkgraphCon,directed=FALSE)/(0.5+((log(con_nodes)-0.5772156649)/log(con_k))))
      #Clustering coefficient.
      con_Cl <- transitivity(networkgraphCon,type="globalundirected",isolate="zero")
      #Log ratio of clustering coefficient to its random counterpart.
      con_lr_Cl <- log(transitivity(networkgraphCon,type="globalundirected",isolate="zero")/(con_k/con_nodes))
    }
    #Add all network parameters to a dataframe.
    networkAnalysis <- rbindlist(list(networkAnalysis,list(networkFile,WS,meanLU,cov_nodes,cov_C,cov_Edge,cov_meanStrength,cov_zeta,cov_lr_zeta,cov_M,cov_L,cov_Cl,con_nodes,con_C,con_Edge,con_meanStrength,con_zeta,con_lr_zeta,con_M,con_L,con_Cl)))
    print(paste(networkFile,WS,meanLU,cov_nodes,cov_C,cov_Edge,cov_meanStrength,cov_zeta,cov_lr_zeta,cov_M,cov_lr_L,cov_lr_Cl,con_nodes,con_C,con_Edge,con_meanStrength,con_zeta,con_lr_zeta,con_M,con_lr_L,con_lr_Cl))
  }
}
colnames(networkAnalysis) <- c("filename","Watershed","meanLU","cov_nodes","cov_C","cov_Edge","cov_meanStrength","cov_zeta","cov_lr_zeta","cov_M","cov_L","cov_Cl","con_nodes","con_C","con_Edge","con_meanStrength","con_zeta","con_lr_zeta","con_M","con_L","con_Cl")
networkAnalysis <- dplyr::arrange(networkAnalysis,Watershed,meanLU)
#Remove any watersheds with spurious results due to insufficient statistics.
networkAnalysis <- subset(networkAnalysis,Watershed!="LakeTahoe")
#Set infinities to null values.
networkAnalysis[networkAnalysis=="-Inf"] <- NA
networkAnalysis[networkAnalysis=="Inf"] <- NA

#Regression between network parameters.
library(Hmisc)
library(corrplot)
library("PerformanceAnalytics")
#Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
chart.Correlation(networkAnalysis[,c("meanLU","cov_nodes","cov_C","cov_Edge","cov_meanStrength","cov_zeta","cov_lr_zeta","cov_M","cov_lr_L")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","con_nodes","con_C","con_Edge","con_meanStrength","con_zeta","con_lr_zeta","con_M","con_lr_L")], histogram=FALSE, method="spearman")

#Check how network topologies change within watersheds.
test <- networkAnalysis
test$diff_meanLU <- ave(test$meanLU, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_cov_nodes <- ave(test$cov_nodes, test$Watershed, FUN=function(x) c(NA,diff(x)))
test$diff_con_nodes <- ave(test$con_nodes, test$Watershed, FUN=function(x) c(NA,diff(x)))
test$diff_cov_C <- ave(test$cov_C, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_con_C <- ave(test$con_C, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_cov_zeta <- ave(test$cov_zeta, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_con_zeta <- ave(test$con_zeta, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_cov_lr_zeta <- ave(test$cov_lr_zeta, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_con_lr_zeta <- ave(test$con_lr_zeta, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_cov_meanStrength <- ave(test$cov_meanStrength, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_con_meanStrength <- ave(test$con_meanStrength, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_cov_M <- ave(test$con_M, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_con_M <- ave(test$cov_M, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_cov_L <- ave(test$cov_L, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_con_L <- ave(test$con_L, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_cov_Cl <- ave(test$cov_Cl, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_con_Cl <- ave(test$con_Cl, test$Watershed, FUN=function(x) c(NA, diff(x)))

chart.Correlation(test[,c("diff_meanLU","diff_cov_nodes","diff_cov_C","diff_cov_meanStrength","diff_cov_zeta","diff_cov_lr_zeta","diff_cov_M","diff_cov_L","diff_cov_Cl")], histogram=FALSE, method="spearman")
chart.Correlation(test[,c("diff_meanLU","diff_con_nodes","diff_con_C","diff_con_meanStrength","diff_con_zeta","diff_con_lr_zeta","diff_con_M","diff_con_L","diff_con_Cl")], histogram=FALSE, method="spearman")
