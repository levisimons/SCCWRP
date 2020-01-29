rm(list=ls())
require("plyr")
require(dplyr)
require(stringr)
require(tidyr)
require(naniar)
require(taxize)

options(ENTREZ_KEY="3f9e29c7c03842f72cf7523e34390d9f2208")
wd <- "~/Desktop/SCCWRP/Metagenomics/"
setwd(wd)

#Read in metagenomic count tables and format them as presence/absence tables.
communityInputRawPlate1 <- read.table("18SV9P1TableWithTaxonomy.txt", header=T, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")
communityInputRawPlate2 <- read.table("18SV9P2TableWithTaxonomy.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=F, encoding = "UTF-8")
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
uniqueOTUs <- suppressWarnings(separate(uniqueOTUs,'ConsensusLineage',c("Rank1", "Rank2","Rank3","Rank4","Rank5","Rank6","Rank7"),sep=";", extra="drop"))
uniqueOTUs$Rank7 <- trimws(uniqueOTUs$Rank7,which="left") #Remove starting blank space from genus names
uniqueOTUs <- uniqueOTUs[!duplicated(uniqueOTUs$OTUID),]
#Filter out ambiguous taxonomies
uniqueOTUs <- uniqueOTUs[uniqueOTUs$Rank1!="Unassigned" & uniqueOTUs$Rank1!="Ambiguous_taxa",]
ambiguousList <- c("Incertae Sedis","metagenome","sp.","environmental","eukaryote","uncultured","soil","Ambiguous_taxa","group","cf.","aff.","gen.","marine","cf","unidentified","Uncultured","invertebrate")
ambiguousList <- as.list(ambiguousList)
uniqueOTUs <- data.frame(lapply(uniqueOTUs, trimws), stringsAsFactors = FALSE)
uniqueOTUs <- replace_with_na_all(data=uniqueOTUs,condition=~.x %in% as.list(ambiguousList))
#Get the furthest resolved taxonomic level.
uniqueOTUs$LeafTaxa <- apply(uniqueOTUs[,!colnames(uniqueOTUs) %in% c("FullTaxonomy")], 1, function(x) tail(na.omit(x), 1))

#Read in SCCWRP BMI taxonomic family-level ids in order to subset read data and focus analyses on just BMIs.
BMIList <- read.table("bug_STE.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8",na.strings=c("","NA"))
BMIList <- unique(BMIList$FamilyNames)
BMIList <- na.omit(BMIList)

BMITaxonomies <- data.frame()
i=0
for(name in BMIList){
  tmp <- classification(name,db="ncbi",rows=1)
  #Sys.sleep(0.5)
  if(nrow(as.data.frame(tmp[1]))>1){
    tmp <- as.data.frame(tmp[1])
    colnames(tmp) <- c("taxa","rank","id")
    tmp <- tmp[tmp$rank!="no rank",]
    if(nrow(tmp)>1){
      tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
      colnames(tmp) <- as.character(unlist(tmp["rank",]))
      tmp <- tmp[!row.names(tmp) %in% c("rank"),]
      rownames(tmp) <- name
      tmp$LeafTaxa <- name
      BMITaxonomies <- dplyr::bind_rows(BMITaxonomies,tmp)
      i=i+1
      print(paste("BMIs",i,length(BMIList)))
    } else{
      tmp <- classification(name,db="gbif",rows=1)
      #Sys.sleep(0.5)
      if(nrow(as.data.frame(tmp[1]))>1){
        tmp <- as.data.frame(tmp[1])
        colnames(tmp) <- c("taxa","rank","id")
        tmp <- tmp[tmp$rank!="no rank",]
        if(nrow(tmp)>1){
          tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
          colnames(tmp) <- as.character(unlist(tmp["rank",]))
          tmp <- tmp[!row.names(tmp) %in% c("rank"),]
          rownames(tmp) <- name
          tmp$LeafTaxa <- name
          BMITaxonomies <- dplyr::bind_rows(BMITaxonomies,tmp)
          i=i+1
          print(paste("BMIs",i,length(BMIList)))
        }
      }
    }
  }
}

#Subset 18Sv9 OTU table to only contain taxonomic BMI data from SCCWRP
BMITaxonomies <- na.omit(unique(unlist(BMITaxonomies)))
BMITaxonomies <- BMITaxonomies[!grepl("\\d",BMITaxonomies)] #Remove numerical taxa terms.
uniqueBMIs <- uniqueOTUs[grep(paste(BMITaxonomies,collapse="|"),uniqueOTUs$FullTaxonomy),]
#
BMIList <- na.omit(unique(uniqueBMIs$LeafTaxa))
BMITaxonomies <- data.frame()
i=0
for(name in BMIList){
  tmp <- classification(name,db="ncbi",rows=1)
  #Sys.sleep(0.5)
  if(nrow(as.data.frame(tmp[1]))>1){
    tmp <- as.data.frame(tmp[1])
    colnames(tmp) <- c("taxa","rank","id")
    tmp <- tmp[tmp$rank!="no rank",]
    if(nrow(tmp)>1){
      tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
      colnames(tmp) <- as.character(unlist(tmp["rank",]))
      tmp <- tmp[!row.names(tmp) %in% c("rank"),]
      rownames(tmp) <- name
      tmp$LeafTaxa <- name
      BMITaxonomies <- dplyr::bind_rows(BMITaxonomies,tmp)
      i=i+1
      print(paste("BMIs",i,length(BMIList)))
    } else{
      tmp <- classification(name,db="gbif",rows=1)
      #Sys.sleep(0.5)
      if(nrow(as.data.frame(tmp[1]))>1){
        tmp <- as.data.frame(tmp[1])
        colnames(tmp) <- c("taxa","rank","id")
        tmp <- tmp[tmp$rank!="no rank",]
        if(nrow(tmp)>1){
          tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
          colnames(tmp) <- as.character(unlist(tmp["rank",]))
          tmp <- tmp[!row.names(tmp) %in% c("rank"),]
          rownames(tmp) <- name
          tmp$LeafTaxa <- name
          BMITaxonomies <- dplyr::bind_rows(BMITaxonomies,tmp)
          i=i+1
          print(paste("BMIs",i,length(BMIList)))
        }
      }
    }
  }
}
BMITaxonomies <- dplyr::left_join(BMITaxonomies,uniqueBMIs,by=c("LeafTaxa"))
BMITaxonomies <- BMITaxonomies[BMITaxonomies$kingdom=="Metazoa",]
BMITaxonomies <- BMITaxonomies[c("OTUID","kingdom","phylum","class","order","family","genus","species")]

write.table(BMITaxonomies,"BMITaxonomies18SV9.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Read in SCCWRP algal taxonomic phyla-level ids in order to subset read data and focus analyses on just algae.
AlgaeList <- read.table("algae_STE.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8",na.strings=c("","NA"))
AlgaeList <- na.omit(unique(AlgaeList$Phylum))
AlgaeList <- na.omit(AlgaeList)
AlgaeList <- AlgaeList[!(AlgaeList %in% c("Eukaryota","Cyanobacteria","Bacteria"))]
AlgalTaxonomies <- data.frame()
i=0
for(name in AlgaeList){
  tmp <- classification(name,db="ncbi",rows=1)
  #Sys.sleep(0.5)
  if(nrow(as.data.frame(tmp[1]))>1){
    tmp <- as.data.frame(tmp[1])
    colnames(tmp) <- c("taxa","rank","id")
    tmp <- tmp[tmp$rank!="no rank",]
    if(nrow(tmp)>1){
      tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
      colnames(tmp) <- as.character(unlist(tmp["rank",]))
      tmp <- tmp[!row.names(tmp) %in% c("rank"),]
      rownames(tmp) <- name
      tmp$LeafTaxa <- name
      AlgalTaxonomies <- dplyr::bind_rows(AlgalTaxonomies,tmp)
      i=i+1
      print(paste("Algae",i,length(AlgaeList)))
    } else{
      tmp <- classification(name,db="gbif",rows=1)
      #Sys.sleep(0.5)
      if(nrow(as.data.frame(tmp[1]))>1){
        tmp <- as.data.frame(tmp[1])
        colnames(tmp) <- c("taxa","rank","id")
        tmp <- tmp[tmp$rank!="no rank",]
        if(nrow(tmp)>1){
          tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
          colnames(tmp) <- as.character(unlist(tmp["rank",]))
          tmp <- tmp[!row.names(tmp) %in% c("rank"),]
          rownames(tmp) <- name
          tmp$LeafTaxa <- name
          AlgalTaxonomies <- dplyr::bind_rows(AlgalTaxonomies,tmp)
          i=i+1
          print(paste("Algae",i,length(AlgaeList)))
        }
      }
    }
  }
}

#Subset 18Sv9 OTU table to only contain taxonomic algal data from SCCWRP
AlgalTaxonomies <- na.omit(unique(unlist(AlgalTaxonomies)))
AlgalTaxonomies <- AlgalTaxonomies[!grepl("\\d",AlgalTaxonomies)]
uniqueAlgae <- uniqueOTUs[grep(paste(AlgalTaxonomies,collapse="|"),uniqueOTUs$FullTaxonomy),]
#
AlgaeList <- na.omit(unique(uniqueAlgae$LeafTaxa))
AlgalTaxonomies <- data.frame()
i=0
for(name in AlgaeList){
  tmp <- classification(name,db="ncbi",rows=1)
  #Sys.sleep(0.5)
  if(nrow(as.data.frame(tmp[1]))>1){
    tmp <- as.data.frame(tmp[1])
    colnames(tmp) <- c("taxa","rank","id")
    tmp <- tmp[tmp$rank!="no rank",]
    if(nrow(tmp)>1){
      tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
      colnames(tmp) <- as.character(unlist(tmp["rank",]))
      tmp <- tmp[!row.names(tmp) %in% c("rank"),]
      rownames(tmp) <- name
      tmp$LeafTaxa <- name
      AlgalTaxonomies <- dplyr::bind_rows(AlgalTaxonomies,tmp)
      i=i+1
      print(paste("Algae",i,length(AlgaeList)))
    } else{
      tmp <- classification(name,db="gbif",rows=1)
      #Sys.sleep(0.5)
      if(nrow(as.data.frame(tmp[1]))>1){
        tmp <- as.data.frame(tmp[1])
        colnames(tmp) <- c("taxa","rank","id")
        tmp <- tmp[tmp$rank!="no rank",]
        if(nrow(tmp)>1){
          tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
          colnames(tmp) <- as.character(unlist(tmp["rank",]))
          tmp <- tmp[!row.names(tmp) %in% c("rank"),]
          rownames(tmp) <- name
          tmp$LeafTaxa <- name
          AlgalTaxonomies <- dplyr::bind_rows(AlgalTaxonomies,tmp)
          i=i+1
          print(paste("Algae",i,length(AlgaeList)))
        }
      }
    }
  }
}
AlgalTaxonomies[] <- lapply(AlgalTaxonomies,as.character)
AlgalTaxonomies <- AlgalTaxonomies[!AlgalTaxonomies$kingdom %in% c("Chromista","Metazoa"),]
AlgalTaxonomies <- dplyr::left_join(AlgalTaxonomies,uniqueAlgae,by=c("LeafTaxa"))
AlgalTaxonomies <- AlgalTaxonomies[c("OTUID","kingdom","phylum","class","order","family","genus","species")]
write.table(AlgalTaxonomies,"AlgaeTaxonomies18SV9.txt",quote=FALSE,sep="\t",row.names = FALSE)
