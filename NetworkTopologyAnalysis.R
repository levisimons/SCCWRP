#Regression between network parameters.
library(Hmisc)
library(corrplot)
library("PerformanceAnalytics")
library(relaimpo)

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

#Run through correlations between upstream land use and topology
varList <- c("meanLU","Pm","M","meanStrength","zeta","C","N","K","networkEdgecount","cov_k","con_k","cov_N","con_N","networkEdgecountCov","networkEdgecountCon","cov_C","con_C","cov_M","con_M","meanStrength_Cov","meanStrength_Con","zeta_Cov","zeta_Con","l_cov_rL","l_con_rL")
networkCorrelations <- data.frame()
for(var in varList){
  topologyCorrelation <- cor.test(networkAnalysis$meanLU,networkAnalysis[,var],method="pearson")
  resilienceCorrelation <- cor.test(networkAnalysis$lambda_network_m,networkAnalysis[,var],method="pearson")
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
write.table(networkCorrelations,"NetworkTopologyCorrelationsLUSweepCAWatershedPerm100Networks.txt",quote=FALSE,sep="\t",row.names = FALSE)

#What factors are underlying variations in network topology relationships?
networkModel <- lm(lambda_network_m ~ M+C+zeta+meanStrength,data=networkAnalysis)
summary(aov(networkModel))
printCoefmat(coef(summary(step(networkModel))))
calc.relimp(networkModel,type="lmg",rela=FALSE)

#What factors are underlying variations in taxonomic richness?
SCCWRPRaw <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
tmp <- as.data.frame(table(SCCWRPRaw$Watershed))
colnames(tmp) <- c("Watershed","Nsamples")
tmp$Watershed <- as.character(tmp$Watershed)
SCCWRP <- join(SCCWRPRaw,tmp)
SCCWRP <- filter(SCCWRP,Nsamples>=30)
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
