#generate simulation scenario tables
library(xtable)

sc <- c(1:18)
nGenes <- rep(c(10000,1000),each=9)
nSamples <- rep(rep(c(8,4,16),each=3),2)
pDiff <- rep(c(0.1,0.3,0.01),6)

sc_table <- cbind(sc, nGenes, nSamples, pDiff)
xtable(sc_table)
