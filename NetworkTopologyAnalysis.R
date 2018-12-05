#Regression between network parameters.
library(Hmisc)
library(corrplot)
library("PerformanceAnalytics")
library(relaimpo)

setwd("~/Desktop/SCCWRP")
networkAnalysis <- read.table("LUSweepCAWatershedPerm100Networks.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
chart.Correlation(networkAnalysis[,c("meanLU","zeta_Cov","l_cov_rL","l_cov_rM","cov_C","meanStrength_Cov","Pm","lambda_network_m","Network_size")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta_Con","l_con_rL","l_con_rM","con_C","meanStrength_Con","Pm","lambda_network_m","Network_size")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta","l_rM","M","lambda_network_m","Network_size")], histogram=FALSE, method="spearman")

#Testing models relating network parameters
TopologyModel <- lm(lambda_network_m ~ con_M+cov_M+cov_C+con_C+l_cov_rL+l_con_rL+meanStrength_Cov+meanStrength_Con+Pm+Network_size+zeta_Cov+zeta_Con, data=networkAnalysis)
summary(TopologyModel)
printCoefmat(coef(summary(step(TopologyModel)))) #Determine parameters to drop from model via AIC scores.
calc.relimp(TopologyModel,type="lmg",rela=TRUE)

#Run through correlations between upstream land use and topology
varList <- c("meanLU","Pm","Network_size","cov_C","con_C","cov_M","con_M","meanStrength_Cov","meanStrength_Con","zeta_Cov","zeta_Con","l_cov_rL","l_con_rL")
for(var in varList){
  topologyCorrelation <- cor.test(networkAnalysis$meanLU,networkAnalysis[,var])
  resilienceCorrelation <- cor.test(networkAnalysis$lambda_network_m,networkAnalysis[,var])
  print(paste("meanLU and",var,topologyCorrelation$estimate,topologyCorrelation$p.value))
  print(paste("lambda_network_m and",var,resilienceCorrelation$estimate,resilienceCorrelation$p.value))
}
