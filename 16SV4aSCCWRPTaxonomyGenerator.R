rm(list=ls())
require("plyr")
require(dplyr)
require(stringr)
require(tidyr)
require(naniar)
require(taxize)
require(textclean)

##Only run this once to generate full taxonomy files for prokaryotic samples.
#Read in metagenomic count tables and format them as presence/absence tables.
#Use 16SV4a reads.

options(ENTREZ_KEY="3f9e29c7c03842f72cf7523e34390d9f2208")
wd <- "~/Desktop/SCCWRP/Metagenomics/"
setwd(wd)

#Read in metagenomic count tables and format them as presence/absence tables.
communityInputRawPlate1 <- read.table("16SV4aP1TableWithTaxonomy.txt", header=T, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")
communityInputRawPlate2 <- read.table("16SV4aP2TableWithTaxonomy.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=F, encoding = "UTF-8")
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
uniqueOTUs <- data.frame(lapply(uniqueOTUs, trimws), stringsAsFactors = FALSE) #Filter out whitespace for text management.
#Filter out ambiguous taxonomies
uniqueOTUs <- uniqueOTUs[uniqueOTUs$Rank1!="Unassigned",]
rankList <- grep("Rank",colnames(uniqueOTUs),value=T)
#Get initial list of taxonomic names.
#Define list of ambiguous taxonomic terms.
ambiguousList <- c("Incertae Sedis","metagenome","environmental","organism",
                   "eukaryote","uncultured","soil","group","marine","cf","unidentified","Uncultured",
                   "invertebrate"," bacterium","Unassigned","clone","compost","symbiont","manure","Chloroplast",
                   "Unknown","agricultural","algae"," alpha","cluster","anaerobic",
                   "Antarctic"," beta","bioreactor","candidate","Chloroplast",
                   "clade"," delta","denitrifying","Ambiguous",
                   "diatom","division","endolithic","endophytic","enrichment","epidermidis",
                   " epsilon","Family XI","Family XII","Family XIII","Family XVIII","FBP",
                   "fecal","flagellate","forest","fragile"," gamma","Green","groundwater","gut","hgcI clade",
                   "JGI","lake","lineage","Lineage IIa","Lineage IIb","Lineage IIc",
                   "Lineage IV","low","microorganism","minor","permafrost","phototrophic","prokaryote",
                   "RFB","SCGC","sediment","sludge","subdivision","Termite","UW","Candidatus")
#Determine full list of ambiguous taxonomic names which contain at least one of the ambiguous terms or a non-alphabetic character.
uniqueOTUs <- data.frame(lapply(uniqueOTUs, trimws), stringsAsFactors = FALSE)
uniqueOTUs[,colnames(uniqueOTUs) %in% rankList] <- as.data.frame(lapply(uniqueOTUs[,colnames(uniqueOTUs) %in% rankList], function(x) replace(x, grep(paste0(ambiguousList,collapse="|"), x), NA)))

#Remove remaining ambiguous taxa terms with NA
tmp <- uniqueOTUs[,colnames(uniqueOTUs) %in% rankList]
tmp[tmp=="bacterium"] <- NA
uniqueOTUs[,colnames(uniqueOTUs) %in% rankList] <- tmp[,colnames(tmp) %in% rankList]

#Get the furthest resolved taxonomic level.
uniqueOTUs$LeafTaxa <- apply(uniqueOTUs[,!colnames(uniqueOTUs) %in% c("FullTaxonomy")], 1, function(x) tail(na.omit(x), 1))

#Go through unique taxa names, resolved to the furthest taxonomic level, and try to obtain the full taxonomies.
BacteriaList <- unique(uniqueOTUs$LeafTaxa)
i=0
BacteriaTaxonomies <- data.frame()
for(name in BacteriaList){
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
      BacteriaTaxonomies <- dplyr::bind_rows(BacteriaTaxonomies,tmp)
      i=i+1
      print(paste("Bacteria",i,length(BacteriaList)))
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
          BacteriaTaxonomies <- dplyr::bind_rows(BacteriaTaxonomies,tmp)
          i=i+1
          print(paste("Bacteria",i,length(BacteriaList)))
        }
      }
    }
  }
}

#Subset the 16Sv4a table to only contain resolved taxonomic data.
BacteriaTaxonomies <- BacteriaTaxonomies[BacteriaTaxonomies$superkingdom=="Bacteria" | BacteriaTaxonomies$superkingdom=="Archaea",]
#rankList <- colnames(BacteriaTaxonomies)
BacteriaTaxonomies[] <- lapply(BacteriaTaxonomies,as.character)
#tmp <- trimws(na.omit(unique(unlist(BacteriaTaxonomies))),which="both")
#for(term in grep(" bacterium",tmp,value=T)){
#  BacteriaTaxonomies[BacteriaTaxonomies==term] <- NA
#}
#Filter out eukaryotes, organize taxonomic columns, and write output.
#BacteriaTaxonomies <- BacteriaTaxonomies[!BacteriaTaxonomies$superkingdom== "Eukaryota",]
BacteriaTaxonomies <- dplyr::left_join(uniqueOTUs,BacteriaTaxonomies,by="LeafTaxa")
BacteriaTaxonomies <- BacteriaTaxonomies[!is.na(BacteriaTaxonomies$superkingdom),]
#rankList <- c(c("OTUID","LeafTaxa"),rankList)
#BacteriaTaxonomies <- BacteriaTaxonomies[,rankList]
#BacteriaTaxonomies$LeafTaxa.1 <- NULL
write.table(BacteriaTaxonomies,"BacteriaTaxonomies16SV4a.txt",quote=FALSE,sep="\t",row.names = FALSE)
#
BacteriaTaxonomies <- read.table("BacteriaTaxonomies16SV4a.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")
