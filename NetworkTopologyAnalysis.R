#Regression between network parameters.
require(Hmisc)
require(corrplot)
require("PerformanceAnalytics")
require(relaimpo)
require(ggmap)
require(maps)
require(mapview)
require(mapdata)
require(munsell)
require(leaflet)
require(devtools)
require(webshot)
require(viridis)
require(stringr)

setwd("~/Desktop/SCCWRP")
# Generate co-occurrence networks from relative abundance of benthic macroinvertebrates
# in stream communities grouped on a watershed by upstream land use basis using the following script:
# https://github.com/levisimons/SCCWRP/blob/master/StreamNetworkGeneratorVcluster.R
# Read in topological analysis of that co-occurrence network data generated with this script:
# https://github.com/levisimons/SCCWRP/blob/master/NetworkAnalysisVcluster.R
networkAnalysis <- read.table("LUSweepCAWatershedPerm100Networks.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Calculate network connectance.
networkAnalysis$C <- (2*networkAnalysis$networkEdgecount)/(networkAnalysis$N*(networkAnalysis$N-1))
#Calculate the degree heterogeneity of the equivalent random network: http://barabasi.com/f/624.pdf
networkAnalysis$zeta_rand <- sqrt((networkAnalysis$N*networkAnalysis$K - networkAnalysis$K^2)/(networkAnalysis$N*networkAnalysis$K^4))
#Get the watershed name for each network
networkAnalysis$Watershed <- str_match(networkAnalysis$filename, "LUSweepCAWatershedPerm100(.*?)Years(.*?)Reps(.*?)MeanLU(.*?)nTaxa(.*?)Network.txt")[,2]

#Search for watersheds by frequency.
watershedTable <- as.data.frame(table(networkAnalysis$Watershed))
watershedSubsample <- subset(watershedTable,watershedTable$Freq>=100)

#Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
chart.Correlation(networkAnalysis[,c("meanLU","zeta_Cov","cov_M","cov_ModGroups","cov_C","meanStrength_Cov","cov_N","cov_k","networkEdgecountCov","lambda_network_m")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta_Con","con_M","con_ModGroups","con_C","meanStrength_Con","con_N","con_k","networkEdgecountCon")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta","K","M","ModGroups","C","Pm","lambda_network_m","N","meanStrength","networkEdgecount")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta","M","C","Pm","lambda_network_m","meanStrength")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta_Cov","cov_M","l_cov_rM","cov_C","meanStrength_Cov","lambda_network_m")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta_Con","con_M","l_con_rM","con_C","meanStrength_Con","lambda_network_m")], histogram=FALSE, method="spearman")

#Testing models relating network parameters
TopologyModel <- lm(lambda_network_m ~ Pm+zeta+M+C+meanStrength, data=networkAnalysis)
summary(aov(TopologyModel))
printCoefmat(coef(summary(step(TopologyModel)))) #Determine parameters to drop from model via AIC scores.
calc.relimp(TopologyModel,type="lmg",rela=FALSE)

#What factors are underlying variations in network topology relationships?
networkModel <- lm(lambda_network_m ~ N+M+C+zeta+meanStrength,data=networkAnalysis)
summary(aov(networkModel))
printCoefmat(coef(summary(step(networkModel))))
calc.relimp(networkModel,type="lmg",rela=FALSE)
networkPlot <- ggplot(networkAnalysis,aes(x=ModGroups,y=M,color=lambda_network_m))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=zeta))
networkPlot+xlab("ModGroups")+ylab("M")+scale_color_gradientn("lambda_network_m",colours = plasma(10))

#What factors are underlying variations in taxonomic richness?
SCCWRPRaw <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
tmp <- as.data.frame(table(SCCWRPRaw$Watershed))
colnames(tmp) <- c("Watershed","Nsamples")
tmp$Watershed <- as.character(tmp$Watershed)
SCCWRPRaw <- join(SCCWRPRaw,tmp)
SCCWRP <- filter(SCCWRPRaw,Nsamples>=30)
TaxaModel <- glm(nTaxa ~ Watershed+Year+altitude+LU+Watershed*Year*altitude*LU, data=SCCWRP)
summary(aov(TaxaModel))
summary(TaxaModel)
printCoefmat(coef(summary(step(TaxaModel)))) #Determine parameters to drop from model via AIC scores.
calc.relimp(TaxaModel,type="lmg",rela=TRUE)

#How does taxonomic richness vary with year between low and high altitude sites?
SCCWRPLowAltitude <- subset(SCCWRP,SCCWRP$altitude < median(SCCWRP$altitude))
SCCWRPHighAltitude <- subset(SCCWRP,SCCWRP$altitude >= median(SCCWRP$altitude))
dev.off()
plot(SCCWRPLowAltitude$Year,SCCWRPLowAltitude$CSCI,type='p',xlim=c(min(SCCWRP$Year),max(SCCWRP$Year)), ylim=c(min(SCCWRP$CSCI,na.rm=TRUE),max(SCCWRP$CSCI,na.rm=TRUE)),xlab="Year",ylab="Number of unique taxa per sample",col="blue")
abline(lm(CSCI~Year,data=SCCWRPLowAltitude),col="blue")
par(new=T)
plot(SCCWRPHighAltitude$Year,SCCWRPHighAltitude$CSCI,type='p',xlim=c(min(SCCWRP$Year),max(SCCWRP$Year)), ylim=c(min(SCCWRP$CSCI,na.rm=TRUE),max(SCCWRP$CSCI,na.rm=TRUE)),xlab="Year",ylab="Number of unique taxa per sample",col="red")
abline(lm(CSCI~Year,data=SCCWRPHighAltitude),col="red")
legend("topright", c(paste("Above",median(SCCWRP$altitude),"m"),paste("Below",median(SCCWRP$altitude),"m")),col=c("red","blue"),lwd=5)
cor.test(SCCWRPLowAltitude$Year,SCCWRPLowAltitude$CSCI)
cor.test(SCCWRPHighAltitude$Year,SCCWRPHighAltitude$CSCI)

#Run through correlations between upstream land use and topology
varList <- c("meanLU","M","meanStrength","zeta","C","N","K","networkEdgecount")
networkCorrelations <- data.frame()
#Find watersheds that fit within certain geographic bounds
watershedsByGeography <- subset(SCCWRP, Longitude >= mean(SCCWRP$Longitude))
#Subset networks by geography
networkAnalysisSubset <- subset(networkAnalysis, Watershed %in% unique(watershedsByGeography$Watershed))
for(var in varList){
  topologyCorrelation <- cor.test(networkAnalysisSubset$meanLU,networkAnalysisSubset[,var],method="pearson")
  resilienceCorrelation <- cor.test(networkAnalysisSubset$lambda_network_m,networkAnalysisSubset[,var],method="pearson")
  print(paste("meanLU and",var,topologyCorrelation$estimate,topologyCorrelation$p.value))
  print(paste("lambda_network_m and",var,resilienceCorrelation$estimate,resilienceCorrelation$p.value))
  dat <- data.frame()
  dat[1,1] <- var
  dat[1,2] <- topologyCorrelation$estimate
  dat[1,3] <- topologyCorrelation$p.value
  #Resilience is -1* the real component of the most positive eigenvalue
  dat[1,4] <- -1*resilienceCorrelation$estimate
  dat[1,5] <- resilienceCorrelation$p.value
  networkCorrelations <- rbind(networkCorrelations,dat)
}
colnames(networkCorrelations) <- c("Variable","CorrelationWithMeanLU","pCorrelationWithMeanLU","CorrelationWithResilience","pCorrelationWithResilience")
#write.table(networkCorrelations,"NetworkTopologyCorrelationsLUSweepCAWatershedPerm100Networks.txt",quote=FALSE,sep="\t",row.names = FALSE)


#To map SCCWRP data.
MapCoordinates <- watershedsByGeography
#MapCoordinates <- subset(MapCoordinates,MapCoordinates$Watershed %in% watershedSubsample$Var1)
MapCoordinates <- MapCoordinates[!is.na(MapCoordinates$Latitude) & !is.na(MapCoordinates$Longitude),]
CalMap = leaflet(MapCoordinates) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=plasma(10),domain=MapCoordinates$LU)
CalMap %>% addCircleMarkers(color = ~ColorScale(LU), fill = TRUE,radius=1,fillOpacity = 0.1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  leaflet::addLegend(position="topright", pal=ColorScale,values=~LU,title="% Upstream\nLand Use")
