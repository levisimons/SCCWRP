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

TopologyModel <- lm(lambda_network_m ~ meanLU+M+cov_C+con_C+l_cov_rL+l_con_rL+meanStrength_Cov+meanStrength_Con+Pm+Network_size+zeta_Cov+zeta_Con, data=networkAnalysis)
summary(TopologyModel)
calc.relimp(TopologyModel,type="lmg",rela=TRUE)
