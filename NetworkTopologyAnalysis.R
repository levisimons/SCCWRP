#Regression between network parameters.
library(Hmisc)
library(corrplot)
library("PerformanceAnalytics")
library(relaimpo)

setwd("~/Desktop/SCCWRP")
networkAnalysis <- read.table("LUSweepCAWatershedPerm100Networks.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
chart.Correlation(networkAnalysis[,c("meanLU","zeta_Cov","l_cov_rL","cov_M","cov_C","meanStrength_Cov","cov_N","cov_k","networkEdgecountCov","lambda_network_m")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta_Con","l_con_rL","con_M","con_C","meanStrength_Con","con_N","con_k","networkEdgecountCon","lambda_network_m")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta","K","M","Pm","lambda_network_m","N","meanStrength","networkEdgecount")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta_Cov","cov_M","cov_C","meanStrength_Cov","lambda_network_m")], histogram=FALSE, method="spearman")
chart.Correlation(networkAnalysis[,c("meanLU","zeta_Con","con_M","con_C","meanStrength_Con","lambda_network_m")], histogram=FALSE, method="spearman")

#Testing models relating network parameters
TopologyModel <- glm(lambda_network_m ~ con_M+cov_M+cov_C+con_C+l_cov_rL+l_con_rL+meanStrength_Cov+meanStrength_Con+Pm+zeta_Cov+zeta_Con, data=networkAnalysis)
summary(TopologyModel)
printCoefmat(coef(summary(step(TopologyModel)))) #Determine parameters to drop from model via AIC scores.
calc.relimp(TopologyModel,type="lmg",rela=TRUE)

#Run through correlations between upstream land use and topology
varList <- c("meanLU","Pm","N","K","networkEdgecount","cov_k","con_k","cov_N","con_N","networkEdgecountCov","networkEdgecountCon","cov_C","con_C","cov_M","con_M","meanStrength_Cov","meanStrength_Con","zeta_Cov","zeta_Con","l_cov_rL","l_con_rL")
networkCorrelations <- data.frame()
for(var in varList){
  topologyCorrelation <- cor.test(networkAnalysis$meanLU,networkAnalysis[,var],method="spearman")
  resilienceCorrelation <- cor.test(networkAnalysis$lambda_network_m,networkAnalysis[,var],method="spearman")
  print(paste("meanLU and",var,topologyCorrelation$estimate,topologyCorrelation$p.value))
  print(paste("lambda_network_m and",var,resilienceCorrelation$estimate,resilienceCorrelation$p.value))
  dat <- data.frame()
  dat[1,1] <- var
  dat[1,2] <- topologyCorrelation$estimate
  dat[1,3] <- topologyCorrelation$p.value
  dat[1,4] <- resilienceCorrelation$estimate
  dat[1,5] <- resilienceCorrelation$p.value
  networkCorrelations <- rbind(networkCorrelations,dat)
}
colnames(networkCorrelations) <- c("Variable","CorrelationWithMeanLU","pCorrelationWithMeanLU","CorrelationWithResilience","pCorrelationWithResilience")
write.table(networkCorrelations,"NetworkTopologyCorrelationsLUSweepCAWatershedPerm100Networks.txt",quote=FALSE,sep="\t",row.names = FALSE)
