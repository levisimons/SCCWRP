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

wd <- "/home/cmb-07/sn1/alsimons/SCCWRP"
wd <- "~/Desktop/SCCWRP/Metagenomics/"
setwd(wd)

#API Key for using ncbi.
options(ENTREZ_KEY="3f9e29c7c03842f72cf7523e34390d9f2208")

#Read in metagenomic count tables and format them as presence/absence tables.
communityInputRawPlate1 <- read.table("rbcLP1TableWithTaxonomy.txt", header=T, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")
communityInputRawPlate2 <- read.table("rbcLP2TableWithTaxonomy.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=F, encoding = "UTF-8")
#Remove spaces from column names for the count tables.
communityInputRawPlate1 <- dplyr::rename(communityInputRawPlate1, OTUID = `OTU ID`)
communityInputRawPlate2 <- dplyr::rename(communityInputRawPlate2, OTUID = `OTU ID`)

#Get a list of unique OTU names from count tables and convert to a data frame.
uniqueOTUs <- rbind(communityInputRawPlate1[,c("OTUID","ConsensusLineage")],communityInputRawPlate2[,c("OTUID","ConsensusLineage")])
uniqueOTUs$OTUID <- as.character(uniqueOTUs$OTUID)

#Removed unassigned taxa.
uniqueOTUs <- uniqueOTUs[uniqueOTUs$ConsensusLineage!="Unassigned",]
#Split OTU names into Domain through Species.
uniqueOTUs$FullTaxonomy <- uniqueOTUs$ConsensusLineage
uniqueOTUs <- suppressWarnings(separate(uniqueOTUs,'ConsensusLineage',c("Rank1", "Rank2","Rank3","Rank4","Rank5","Rank6","Rank7"),sep=";", extra="drop"))
uniqueOTUs <- as.data.frame(lapply(uniqueOTUs,function(x) trimws(x, which="left") ))
uniqueOTUs <- data.frame(lapply(uniqueOTUs, as.character), stringsAsFactors=FALSE)
rankList <- grep("Rank",colnames(uniqueOTUs),value=T)

#Filter to only retain diatoms.
uniqueOTUs <- uniqueOTUs[uniqueOTUs$Rank3=="Bacillariophyceae",]
# Get full taxonomies for all diatom OTUs.
DiatomTaxonomies <- data.frame()
i=0
for(name in unique(uniqueOTUs$Rank7)){
  tmp <- classification(name,db="ncbi",rows=1)
  Sys.sleep(0.5)
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
      DiatomTaxonomies <- dplyr::bind_rows(DiatomTaxonomies,tmp)
      i=i+1
      print(paste("Diatom",i,length(unique(uniqueOTUs$Rank7))))
    } else{
      tmp <- classification(name,db="gbif",rows=1)
      Sys.sleep(0.5)
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
          DiatomTaxonomies <- dplyr::bind_rows(DiatomTaxonomies,tmp)
          i=i+1
          print(paste("Diatom",i,length(unique(uniqueOTUs$Rank7))))
        }
      }
    }
  }
}

#Filter out non-diatom taxonomies
DiatomTaxonomies <- DiatomTaxonomies[DiatomTaxonomies$class=="Bacillariophyceae",]
DiatomTaxonomies[] <- lapply(DiatomTaxonomies,as.character)
#Merge OTUIDs with taxonomies
DiatomTaxonomies$LeafTaxa <- trimws(DiatomTaxonomies$LeafTaxa,which="left") #Remove starting blank space from genus names
uniqueDiatoms <- dplyr::left_join(uniqueOTUs,DiatomTaxonomies,by=c("Rank7"="LeafTaxa"))
write.table(uniqueDiatoms,"DiatomTaxonomiesrbcL.txt",quote=FALSE,sep="\t",row.names = FALSE)
