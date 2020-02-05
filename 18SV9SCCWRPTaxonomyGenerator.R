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

#Get all animal OTUs.
animalOTUs <- uniqueOTUs[uniqueOTUs$Rank4=="Metazoa (Animalia)",]
animalOTUs <- animalOTUs[!is.na(animalOTUs$OTUID),]
# Get full taxonomies for all animal OTUs.
AnimalTaxonomies <- data.frame()
i=0
for(name in unique(animalOTUs$LeafTaxa)){
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
      AnimalTaxonomies <- dplyr::bind_rows(AnimalTaxonomies,tmp)
      i=i+1
      print(paste("Animal",i,length(unique(animalOTUs$LeafTaxa))))
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
          AnimalTaxonomies <- dplyr::bind_rows(AnimalTaxonomies,tmp)
          i=i+1
          print(paste("Animal",i,length(unique(animalOTUs$LeafTaxa))))
        }
      }
    }
  }
}

#Organize animal taxonomies.
rankList <- c("LeafTaxa","species","genus","family","superfamily","suborder","order","infraclass","subclass","class","subphylum","phylum","kingdom","superkingdom")
AnimalTaxonomies <- AnimalTaxonomies[,rankList]
AnimalTaxonomies[] <- lapply(AnimalTaxonomies,as.character)
AnimalTaxonomies <- AnimalTaxonomies[AnimalTaxonomies$kingdom=="Metazoa" & !is.na(AnimalTaxonomies$phylum),]

#Read in SCCWRP BMI taxonomic family-level ids in order to subset read data and focus analyses on just BMIs.
BMIList <- read.table("bug_STE.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8",na.strings=c("","NA"))
BMIList <- data.frame(lapply(BMIList, trimws), stringsAsFactors = FALSE)
BMIList <- unique(BMIList$FamilyNames)
BMIList <- na.omit(BMIList)

#Subset 18Sv9 OTU table to only contain taxonomic BMI data from SCCWRP
BMITaxonomies <- AnimalTaxonomies[apply(AnimalTaxonomies, 1, function(x) any(x %in% BMIList)),]
uniqueBMIs <- uniqueOTUs[grep(paste(unique(BMITaxonomies$LeafTaxa),collapse="|"),uniqueOTUs$FullTaxonomy),]
uniqueBMIs <- dplyr::left_join(uniqueBMIs,BMITaxonomies,by=c("LeafTaxa"))
rankList <- c("OTUID",rankList)
uniqueBMIs <- uniqueBMIs[,rankList]
write.table(uniqueBMIs,"BMITaxonomies18SV9.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Get all plant OTUs.
plantOTUs <- uniqueOTUs[uniqueOTUs$Rank4 %in% c("Chlorophyta","Charophyta"),]
plantOTUs <- plantOTUs[plantOTUs$LeafTaxa!="Chlorophyta" & plantOTUs$LeafTaxa!="Charophyta",]
# Get full taxonomies for all animal OTUs.
PlantTaxonomies <- data.frame()
i=0
for(name in unique(plantOTUs$LeafTaxa)){
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
      PlantTaxonomies <- dplyr::bind_rows(PlantTaxonomies,tmp)
      i=i+1
      print(paste("Plants",i,length(unique(plantOTUs$LeafTaxa))))
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
          PlantTaxonomies <- dplyr::bind_rows(PlantTaxonomies,tmp)
          i=i+1
          print(paste("Plant",i,length(unique(PlantOTUs$LeafTaxa))))
        }
      }
    }
  }
}
#Organize plant taxonomies.
rankList <- c("LeafTaxa","varietas","species","subgenus","genus","subtribe","tribe","subfamily","family","order","superorder","subclass","class","subphylum","phylum","kingdom","superkingdom")
PlantTaxonomies <- PlantTaxonomies[,rankList]
PlantTaxonomies[] <- lapply(PlantTaxonomies,as.character)
PlantTaxonomies <- PlantTaxonomies[PlantTaxonomies$kingdom=="Viridiplantae",]

#Read in SCCWRP algal taxonomic ids in order to subset read data and focus analyses on just algae.
AlgaeList <- read.table("algae_STE.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8",na.strings=c("","NA"))
AlgaeList <- AlgaeList[AlgaeList$Kingdom=="Eukaryota",]
AlgaeList <- AlgaeList[, !duplicated(colnames(AlgaeList))]
AlgaeList <- dplyr::rename(AlgaeList,LeafTaxa = FinalIDassigned)
AlgaeList$FinalID <- NULL
AlgaeList <- AlgaeList[!duplicated(AlgaeList),]
AlgaeList <- na.omit(unique(AlgaeList$LeafTaxa))

#Subset 18Sv9 OTU table to only contain taxonomic BMI data from SCCWRP
AlgaeTaxonomies <- PlantTaxonomies[apply(PlantTaxonomies, 1, function(x) any(x %in% AlgaeList)),]
uniqueAlgae <- uniqueOTUs[grep(paste(unique(AlgaeTaxonomies$LeafTaxa),collapse="|"),uniqueOTUs$FullTaxonomy),]
uniqueAlgae <- dplyr::left_join(uniqueAlgae,AlgaeTaxonomies,by=c("LeafTaxa"))
rankList <- c("OTUID",rankList)
uniqueAlgae <- uniqueAlgae[,rankList]
write.table(uniqueAlgae,"BMITaxonomies18SV9.txt",quote=FALSE,sep="\t",row.names = FALSE)
