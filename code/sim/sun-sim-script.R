#setwd("~/Desktop/ebayes_de/sim")
source("../sim/help_func.R")
real_parms <- readRDS("../sim/trimmed_parms.rds")

##################################################################
#sc18: nGenes=1000, nSample=16, pDiff=0.01, n_sim=5
##################################################################

sc18_alter_pval <- NBsim_scenario(params=real_parms,
                                 nGenes=1000, 
                                 nSample=16, 
                                 pDiff=0.01, nSim=5)

saveRDS(sc18_alter_pval, "results/sc18_alter_pval.rds")

sc18_alter_pval <- readRDS("results/sc18_alter_pval.rds")

#add ebayes
##############################################
#sc18
nGenes=1000; nSample=16; pDiff=0.01
diffPerc = pDiff*100; nRep=nSample/2
##############################################
#sim5
i=5;
data_file <- paste(paste("data/sim_genes",nGenes,"g",nRep, "pDiff",diffPerc, i, sep="_"),"rds",sep=".")

d = as.data.frame(readRDS(data_file))
d_count <- d[,-ncol(d)]
n = ncol(d_count)/2
names(d_count) = c(paste("B73_",1:n, sep=''), paste('Mo17_',1:n, sep=''))
variety = factor(gsub("_[0-9]{1,2}", "", names(d_count)), levels=c("B73","Mo17"))

model = stan_model("../common/normal.stan")
hyperparameters = get_hyperparameters(d=d_count, variety)
saveRDS(hyperparameters, paste0("data/sim_genes_", nGenes,"_g_",nRep, 
                                "_pDiff_",diffPerc, "_", i,"_hyperparameters.rds"))
# Run individual gene analyses
if (parallel <- require(doMC)) {
  registerDoMC()
} else if (parallel <- require(doParallel)) {
  cl = makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
}

parallel <- require(doMC)

start_time <- Sys.time()
analysis = adply(d_count,
                 1,
                 function(x) single_gene_analysis(x),
                 .id = 'gene',
                 .parallel = parallel,
                 .paropts = list(.export=c('single_gene_analysis','model','hyperparameters'), 
                                 .packages='rstan'))
end_time <- Sys.time()
end_time - start_time

rownames(analysis) = rownames(d_count)

saveRDS(analysis, file=paste0("results/","sim_genes_",nGenes,"_g_",nRep, "_pDiff_",diffPerc, "_",i,"_analysis",
                              ".rds"))

ebayes_analysis <- analysis %>% 
  mutate(p_value = if_else(prob_LDE<prob_HDE, prob_LDE, prob_HDE))
true_de <- d$true_de

ebayes_pval <- as.data.frame(cbind(ebayes_analysis$p_value, true_de))
colnames(ebayes_pval) <- c("ebayes_pval", "de")
saveRDS(ebayes_pval, file=paste0("results/","sim_genes_",nGenes,"_g_",nRep, "_pDiff_",diffPerc, "_",i,"_ebayes_pval",
                                 ".rds"))

#plot ROC curve for sc18_sim1

alter_pval <- sc18_alter_pval[[i]]
all_pval <- cbind(alter_pval, ebayes_pval$ebayes_pval)
colnames(all_pval)[ncol(all_pval)] <- c("ebayes_pval")

saveRDS(all_pval, file=paste0("results/","sim_genes_",nGenes,"_g_",nRep, "_pDiff_",diffPerc, "_",i,"_pval",
                                 ".rds"))
sc18_sim5_auc <- plot_roc_all(all_result=all_pval, name="sc18 sim5")
saveRDS(sc18_sim5_auc, file=paste0("results/","sim_genes_",nGenes,"_g_",nRep, "_pDiff_",diffPerc, "_",i,"_auc",
                                  ".rds"))


sc18_auc <- rbind(sc18_sim1_auc, 
                 sc18_sim2_auc,
                 sc18_sim3_auc,
                 sc18_sim4_auc,
                 sc18_sim5_auc)
rownames(sc18_auc) <- paste0("sim",1:5)

saveRDS(sc18_auc, file = paste0("results/","sim_genes_",nGenes,"_g_",nRep, "_pDiff_",diffPerc,"_auc",
                               ".rds"))


