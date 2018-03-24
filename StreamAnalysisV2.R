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

#Read in California Stream Condition Index data
csciData <- read.csv("csci_scored_sites_tbl.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
csciData <- filter(csciData,CSCI!="NA")
names(csciData)[names(csciData)=="StationCode"]<-"SampleStationID"
names(csciData)[names(csciData)=="SAMPLEDATE"]<-"SampleDate"
#Make the date format uniform.
csciData$SampleDate <- as.Date(csciData$SampleDate,format="%m/%d/%y")
#Add in qualifier columns based on California Streams Condition Index.  The cutoff is 0.79.
csciData$CSCIQualifier <- ifelse(csciData$CSCI >= 0.79, "Healthy","Disturbed")
csciData$CSCIQualNum <- ifelse(csciData$CSCI >= 0.79, 1,0)
#Add in qualifier column based on total nitrogen concentrations.  The cutoff is 0.42mg/L.
#csciData$TNQualifier <- ifelse(csciData$`Total_Nitrogen (mg/L)`<=0.42,"Healthy","Disturbed")
#Add in qualifier column based on total phosphate concentrations.  The cutoff is 0.03mg/L.
#csciData$TPQualifier <- ifelse(csciData$`Total_Phosphorous (mg/L)`<=0.03,"Healthy","Disturbed")
#Parse out the year the CSCI was calculated by site.
csciData$Year <- year(csciData$SampleDate)
#Keep only the first replicate.  Most sites only have one replicate.
csciData <- filter(csciData,REPLICATE==1)

#Merge CSCI data onto site data.
GISBiochemcsciData <- merge(GISBiochemData,csciData,by=c("SampleStationID","Year"))

#Read in network statistics generated under different chemical parameter selection windows.
chemNetworks <- read.csv("SCCWRPNetworkAnalysisChemicalThresholdsSummary.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)

#Assign network parameters based on total nitrogen concentrations.
csciData$TN_Con_l_rL <-ifelse(csciData$`Total_Nitrogen (mg/L)`>0.42,chemNetworks[chemNetworks$Prefix=='HTN','Con_l_rL'],chemNetworks[chemNetworks$Prefix=='LTN','Con_l_rL'])  
csciData$TN_Con_l_rCl <-ifelse(csciData$`Total_Nitrogen (mg/L)`>0.42,chemNetworks[chemNetworks$Prefix=='HTN','Con_l_rCl'],chemNetworks[chemNetworks$Prefix=='LTN','Con_l_rCl'])
csciData$TN_Con_l_rM <-ifelse(csciData$`Total_Nitrogen (mg/L)`>0.42,chemNetworks[chemNetworks$Prefix=='HTN','Con_l_rM'],chemNetworks[chemNetworks$Prefix=='LTN','Con_l_rM'])  
csciData$TN_Cov_l_rL <-ifelse(csciData$`Total_Nitrogen (mg/L)`>0.42,chemNetworks[chemNetworks$Prefix=='HTN','Cov_l_rL'],chemNetworks[chemNetworks$Prefix=='LTN','Cov_l_rL'])  
csciData$TN_Cov_l_rCl <-ifelse(csciData$`Total_Nitrogen (mg/L)`>0.42,chemNetworks[chemNetworks$Prefix=='HTN','Cov_l_rCl'],chemNetworks[chemNetworks$Prefix=='LTN','Cov_l_rCl'])  
csciData$TN_Cov_l_rM <-ifelse(csciData$`Total_Nitrogen (mg/L)`>0.42,chemNetworks[chemNetworks$Prefix=='HTN','Cov_l_rM'],chemNetworks[chemNetworks$Prefix=='LTN','Cov_l_rM'])  

#Note that this dataframe still contains all sites with a CSCI.
#It still needs to be subsetted so that only sites where a particular parameter,
#such as total nitrogen, are measured.

#Merge the CSCI data by site and date with that of chemical, biological, and land usage data by site.
#BioChemcsciData <- merge(GISBiochemData,csciData,by=c("SampleStationID","Year"))

#Logistic regression between network parameters and CSCI
library(aod)
library(glmm)
library(rcompanion)
logReg <- glm(formula = CSCIQualNum ~ TN_Con_l_rL+TN_Con_l_rCl+TN_Con_l_rM+TN_Cov_l_rL+TN_Cov_l_rCl+TN_Cov_l_rM, data = csciData, family=binomial(link="logit"))
#Determine pseudo-r^2 and p for logistic model.
nagelkerke(logReg)
#Summary statistics for logistic model.
#log(p/(1-p)) = k*x + b.
#The probability of the binomial model being in a state of 1, versus 0, is p.
#First decide a value of p, such as 0.95 for 95% probability.
#The estimate of variable, in this case l_rL, is k.
#The estimate of the intercept is b.
#Solve for the value of the target variable which gives a probability of p that the binomial model is in a state of 1 versus 0.
#x = -1*([log(p/(1-p))-b]/k) for a given value of p, b, and k.
summary(logReg)

#Find the average value of the parameters most strongly correlated to the top
#two principal components which describe variations in chemical parameter space
#given changes in land usage intensity.  These averages will be used to split
#the site data according to sites above or below the average parameter values.
#parameter1 = strongly correlated parameter to axis composed of
#principal components 1 and 2.
i=0
parameter1 <- data.frame()
parameter1Name <- "Nitrate + Nitrite as N"
parameterSites <- data.frame()
for(site in unique(GISBiochemData$UniqueID)){
  GISBiochemDataSite <- GISBiochemData[GISBiochemData$UniqueID == site,]
  if(parameter1Name %in% GISBiochemDataSite$FinalID){
    i=i+1
    tmp1 <- GISBiochemDataSite[which(GISBiochemDataSite$FinalID==parameter1Name),]
    if(tmp1$Measurement[1]>0){
      parameter1 <- rbind(parameter1,tmp1$Measurement[1])
      parameterSites <- rbind(parameterSites,tmp1)
    }
    #print(paste(tmp1))
  }
}

#Get all of the unique sample site location-year pairings where a particular
#parameter was measured.
uniqueSiteYears <- unique(parameterSites[,c('SampleStationID','Year')])

#Split site data based on the average value of the two parameters.
#Given the most significant chemical factors related to changes in land usage
#subset site data based on those values.
#H1 = high parameter 1 concentration.  L1 = low parameter 1 concentration.
GISBiochemDataH1 <- data.frame()
GISBiochemDataL1 <- data.frame()
for(site in unique(GISBiochemData$UniqueID)){
  GISBiochemDataSite <- GISBiochemData[GISBiochemData$UniqueID == site,]
  if(parameter1Name %in% GISBiochemDataSite$FinalID){
    tmp1 <- GISBiochemDataSite[which(GISBiochemDataSite$FinalID==parameter1Name),]
    if(tmp1$Measurement[1]>parameter1Ave){
      print(paste("H1: ",site,parameter1Name,tmp1$Measurement[1]))
      GISBiochemDataH1 <- rbind(GISBiochemDataH1,GISBiochemDataSite)
    }
    if(tmp1$Measurement[1]<parameter1Ave){
      print(paste("L1: ",site,parameter1Name,tmp1$Measurement[1]))
      GISBiochemDataL1 <- rbind(GISBiochemDataL1,GISBiochemDataSite)
    }
  }
}

#Select a geographic subset data frame from the total merged data set.
selected <- GISBiochemDataH1
suffix <- "HighNitrateNitrite"

#Initialize a data frame where the rows are all of the unique measurements for a given
#subset of the data.
#Order the data frame by measurement name.
eLSAInput <- as.data.frame(unique(selected$FinalID))
colnames(eLSAInput)<-c("FinalID")
eLSAInput <- as.data.frame(eLSAInput[order(as.character(eLSAInput$FinalID)),])
colnames(eLSAInput)<-c("FinalID")

#Add the relative taxa abundances by column to a new dataframe.
#The rows are the unique taxa in a given subset of data.
for(ID in unique(selected$UniqueID)){
  tmp <- filter(selected, UniqueID == ID)[,c(3,4,6)]
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
  print (paste(repMax,repNum,sep=" "))
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

#If you want to determine the average parameter values across the subsetted data set.
library(matrixStats)
eLSAAverage <- eLSAInput
eLSAAverage <- as.data.frame(sapply(eLSAAverage,as.numeric))
eLSAAverage$FinalID <- eLSAInput$FinalID
eLSAAverage$mean <- rowMeans(eLSAAverage[,2:ncol(eLSAAverage)],na.rm=TRUE)
eLSAAverage$median <- rowMedians(as.matrix(eLSAAverage[,2:ncol(eLSAAverage)]),na.rm=TRUE)
eLSAAverage$STDEV <- rowSds(as.matrix(eLSAAverage[,2:ncol(eLSAAverage)]),na.rm=TRUE)

#If you want to output the parameter averages for a given geographic subset of data.
chemID <- unique(chemData$FinalID)
eLSAAverage <- subset(eLSAAverage,eLSAAverage$FinalID %in% chemID)[,c(-2:-(ncol(eLSAAverage)-3))]
colnames(eLSAAverage) <- c("FinalID",paste("mean",suffix),"median","STDEV")
write.table(eLSAAverage,paste("eLSAAverage",suffix,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)

#Aggreate all of the mean physical parameter averages into a single data frame
#for analysis.  Make sure you've already generated these subsetted files first.
suffixList = c("HD1K","MD1K","LD1K","HD5K","MD5K","LD5K","HDWS","MDWS","LDWS")
means <- data.frame(chemID)
colnames(means) <- c("FinalID")
for(suffix in suffixList){
  print(paste("eLSAAverage",suffix,".txt",sep=""))
  parameter <- read.delim(paste("eLSAAverage",suffix,".txt",sep=""),header=TRUE, sep="\t",as.is=T,check.names=FALSE)
  print(parameter[,-c(3:4)])
  parameter <- data.frame(parameter[,-c(3:4)])
  means <- join(means,parameter,by="FinalID")
}
means <- t(means)
colnames(means) <- means[1,]
means <- means[-1,]
means <- as.data.frame(means)

#Output dataframe for use in eLSA.
#Note that the the data needs to have at least two location replicates per time point
#and that the number of replicates per time point needs to be uniform.
#This may involve subsampling data depending on the variation in the number of replicates per time point.
#The first character in an eLSA formatted file needs to be #.
names(eLSAInput)[names(eLSAInput)=="FinalID"]<-"#FinalID"
write.table(eLSAInput,paste("eLSAInput",suffix,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#At this point insert columns with all NA values for each year in order to even
#out the number of replicates per year.

#Read in eLSA output.
#Compute network statistics of the likeliest association networks between taxa.
library(igraph)
library(network)
suffix <- "LDWS"
networkdata <- read.delim(paste("eLSAOutput",suffix,".txt",sep=""),header=TRUE, sep="\t",as.is=T,check.names=FALSE)
#Filter out association network data based on P scores, for the local similarity
#between two factors, with values less than 0.05.
networkdata <- filter(networkdata, P <= 0.01)
names(networkdata)[names(networkdata)=="LS"]<-"weight"
#Filter network data based on local similarity scores.
networkdata <- subset(networkdata,networkdata$weight<0)
algaeID <- unique(algaeData$FinalID)
insectID <- unique(insectData$FinalID)
chemID <- unique(chemData$FinalID)

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

#Plot a corresponding random graph given a number of edges and nodes.
RandGraph <- erdos.renyi.game(network.size(adj),network.edgecount(adj),type="gnm")
plot(RandGraph, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 3, xlab = "Random Network: G(N,L) model")

library(vcd)
library(MASS)
# Get degree distribution of network.
DDN <- degree(networkgraph)
# Fit a poisson distribution to the link distribution of the network
curveFit <- fitdistr(DDN,"exponential")
# Get the fit parameters for the distribution.
# Scale-free networks are exponential and random networks are Poisson distributions.
coef(curveFit)
# Get the log-likelihood for this fit
logLik(curveFit)
# Histogram of node degree distribution.
hist(degree(networkgraph, mode="all"), breaks=1:vcount(networkgraph)-1, main="Histogram of node degree")
# Find the largest clique within the network.
largest_cliques(as.undirected(networkgraph, mode= "collapse", edge.attr.comb=list(weight="sum", "ignore")))
# List nodes by their number of links.
degree(networkgraph)
