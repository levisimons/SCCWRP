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
GISBiochemData <- read.table("GISBiochemData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Ensure that all sites have a CSCI value.
GISBiochemData <- subset(GISBiochemData, CSCI != "NA")

#Randomly subsample the site data by individual sample location-date (UniqueID).
SampleSize=45
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
  print(ID)
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
  print (paste(repMax,repNum,year,sep=" "))
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
      print(paste(k,repLabel,sep=" "))
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
filename

#Read in eLSA output.
#Compute network statistics of the likeliest association networks between taxa.
library(igraph)
library(network)
networkFile = "N45S13R8M0.969482509039757Network.txt"
networkdata <- read.delim(networkFile,header=TRUE, sep="\t",as.is=T,check.names=FALSE)
#Filter out association network data based on P scores, for the local similarity
#between two factors, with values less than 0.05.
networkdata <- filter(networkdata, P <= 0.01)
names(networkdata)[names(networkdata)=="LS"]<-"weight"
#Filter network data based on local similarity scores.
#networkdata <- subset(networkdata,networkdata$weight>0)
#Get unique identifiers for algal, invertebrate, and chemical measurement types.
algae <- subset(GISBiochemData,GISBiochemData$MeasurementType=="Algal relative abundance")
algaeID <- unique(algae$FinalID)
insect <-subset(GISBiochemData,GISBiochemData$MeasurementType=="Invertebrate relative abundances" | GISBiochemData$MeasurementType=="Invertebrate relative abundance")
insectID <- unique(insect$FinalID)
chem <- subset(GISBiochemData,GISBiochemData$MeasurementType!="Algal relative abundance" & GISBiochemData$MeasurementType!="Invertebrate relative abundances" & GISBiochemData$MeasurementType!="Invertebrate relative abundance")
chemID <- unique(chem$FinalID)

#Define a 'not in' function.
'%!in%' <- function(x,y)!('%in%'(x,y))
#Remove some subset of chemical and biological factors as nodes from the network.
#networkdata1 <- subset(networkdata,networkdata$X %!in% chemID & networkdata$Y %in% chemID)
#networkdata2 <- subset(networkdata,networkdata$Y %!in% chemID & networkdata$X %in% chemID)
#networkdata <- rbind(networkdata1,networkdata2)
networkdata <- subset(networkdata,networkdata$X %!in% chemID)
networkdata <- subset(networkdata,networkdata$Y %!in% chemID)

#Generate network graph and begin calculating network parameters.
networkgraph=graph.data.frame(networkdata,directed=FALSE)
l <- layout_with_fr(networkgraph)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
par(mfrow=c(1,1), mar=c(0,0,0,0))
#Option 1: plot just the structure of the network.
V(networkgraph)$label <- ""
plot(networkgraph,vertex.size=3)
#Option 2: plot the network and weight links by LS scores and node size by number of links.
#plot(networkgraph,rescale=F,layout=l*1.0,vertex.size=5*degree(networkgraph),edge.width=abs(E(networkgraph)$weight*10),edge.color=ifelse(E(networkgraph)$weight > 0, "blue","red"))
# Calculate number of groups and the modularity of the network.
cluster_edge_betweenness(networkgraph, weights=NULL,directed=FALSE)
# Calculate the average network path length
mean_distance(networkgraph,directed=FALSE)
# Calculate the clustering coefficient
transitivity(networkgraph,type="globalundirected",isolate="zero")
# Generate adjacency matrix of relative taxa abundance correlations
adj= as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
# Get the number of unique network edges
network.edgecount(adj)
# Get the number of nodes
network.size(adj)
# Get the network density.
network.density(adj)

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
gamma = lambda_network_1/lambda_rand_1
gamma

