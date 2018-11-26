require(plyr)
require(dplyr)
require(ape)
require(vegan)
require(data.table)

#Create co-occurrence networks by watershed and land use decile in the SCCWRP data set.
#setwd("~/Desktop/SCCWRP")
setwd("~/panfs/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage values.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Remove duplicate measures.
GISBioData <- GISBioData[!duplicated(GISBioData[,c("UniqueID","FinalID","Count")]),]
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
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
#Merge in functional feeding groups into sample data.
GISBioData <- join(GISBioData,FFG[,c("FinalID","FunctionalFeedingGroup")],by=c("FinalID"))
FFGCounts <- na.omit(as.data.frame(unique(GISBioData$FunctionalFeedingGroup)))
colnames(FFGCounts) <- c("FunctionalFeedingGroups")
#How many samples per watershed?
groupNum = 25
selectNum = 30
#Select watersheds with a large enough set of samples for analysis.
watersheds <- as.data.frame(table(SCCWRP$Watershed))
colnames(watersheds) <- c("Watershed","Samples")
GISBioData <- join(GISBioData,watersheds,by=c("Watershed"))
#Get samples only found in more heavily sampled watersheds.
GISBioDataLargeWS <- subset(GISBioData,Samples>=selectNum)

LUquantile <- quantile(GISBioDataLargeWS$LU_2000_5K,probs=seq(0,1,0.2))#To get land use quantiles.

set.seed(1)
#Find watersheds with a sufficient number of samples.  Randomly subsample groups of samples within
#a watershed and land use quantile for co-occurrence network generation using the eLSA program.
#Run 100 permutations to generate all of the raw input data to create co-occurrence networks.
#Also generate a shell script containing all of the eLSA commands to generate co-occurrence networks on the cluster.
commandLines <- data.frame()
for(iteration in 1:100){
  for(watershed in unique(GISBioDataLargeWS$Watershed)){
    for(i in 1:length(LUquantile)){
      if(i>1){
        LULow <- LUquantile[i-1]
        LUHigh <- LUquantile[i]
        MidLU <- 0.5*(LULow+LUHigh)
        localSubset <- subset(GISBioDataLargeWS,Watershed==watershed & LU_2000_5K >= LULow & LU_2000_5K <= LUHigh)
        if(length(unique(localSubset$UniqueID)) >= selectNum){
          selected <- filter(localSubset, UniqueID %in% sample(unique(localSubset$UniqueID),groupNum))
          tmp <- selected[,c("UniqueID","LU_2000_5K")]
          tmp <- tmp[!duplicated(tmp$UniqueID),]
          meanLU = mean(na.omit(tmp$LU_2000_5K))
          localTaxa <- unique(selected$FinalID) #Taxa found in a set of samples defined by target watershed and land use.
          
          #Initialize a data frame where the rows are all of the unique measurements for a given
          #subset of the data.
          #Order the data frame by measurement name.
          selected <- arrange(selected,Year,UniqueID)
          eLSAInput <- as.data.frame(localTaxa)
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
          eLSAInput[is.na(eLSAInput)] <- 0
          
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
          filename = paste("LUSweepCAWatershedPerm100",watershed,"Years",spotNum,"Reps",repNum,"MeanLU",meanLU,"nTaxa",length(localTaxa),".txt",sep="")
          commandLine <- data.frame()
          commandLine[1,1] = paste("lsa_compute ",filename," -r ",repNum," -s ",spotNum," ",gsub(".txt","Network.txt",filename),";",sep="")
          print(commandLine[1,1])
          commandLines <- rbind(commandLines,commandLine)
          write.table(eLSAInput,filename,quote=FALSE,sep="\t",row.names = FALSE)
        }
      }
    }
  } 
}
#Output eLSA commands into a single shell script.
write.table(commandLines,"eLSAPerm100.sh",quote=FALSE,sep="\t",row.names = FALSE,col.names = FALSE)
#If the full file of eLSA commands already exists, just read it in here instead of generating it from above.
commandLines <- read.table("eLSAPerm100.sh", header=FALSE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

processorNum = 500 #How many network generation processes will be run in parallel?

##The header dataframe represents the lines to include in each shell script to run them as part of SBATCH on the cluster.
## The header requests a job that runs for 300 hours using one processor and node per job, with a memory request of 2GB per job.
header <- data.frame()
header[1,1] <- "#!/bin/sh"
header[2,1] <- "#SBATCH -t300:00:00"
header[3,1] <- "#SBATCH -n1"
header[4,1] <- "#SBATCH -c1"
header[5,1] <- "#SBATCH --mem=2000m"
header[6,1] <- "#SBATCH -pcegs"
header[7,1] <- "cd ~/panfs/SCCWRP"

#Split eLSA commands into separate shell scripts to run in parallel on a cluster.
for(i in 1:processorNum){
  #print(paste((nrow(commandLines)/processorNum)*(i-1)+1, (nrow(commandLines)/processorNum)*i))
  lowRow <- (nrow(commandLines)/processorNum)*(i-1)+1
  highRow <- (nrow(commandLines)/processorNum)*i
  commandFile <- as.data.frame(commandLines[lowRow:highRow,])
  colnames(commandFile) <- c("V1")
  commandFile <- rbind(header,commandFile) #Add the header to each shell script to run as part of a batch
  filename <- paste("eLSA",i,".sh",sep="")
  write.table(commandFile,filename,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}

