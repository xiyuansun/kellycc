#supplemental script

#setwd("~/Desktop/kellycc/code")
##########################################################
#Step 0: Load the required packages
##########################################################
library(MASS)
library(edgeR)
library(DESeq)
library(DESeq2)
library(sSeq)
library(EBSeq)
library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(Paschold2012)
library(ROCR)
library(rstan)
library(reshape2)

source("./sim/help_func.R")
###########################################################
#Step 1: Read in the Paschold2012 RNA-Seq Count Data
###########################################################

# 2 varieties: B73 and Mo17
# 4 replicates per variety
real_count <- readRDS("./real/data/paschold.rds") 

real_split_names <- unlist(strsplit(colnames(real_count), "[_]"))
real_sample_names <- real_split_names[seq_along(real_split_names)%%2 !=0]
# 2 level factor conditions
real_cond <- factor(real_sample_names)
colnames(real_count) <- real_sample_names
# get gene names and group ready for parameter estimation
real_group <- real_cond
real_geneid <- rownames(real_count)


###########################################################
#Step 2: Trimmed real data
###########################################################

# select genes with average count >=1 
# and at most two zero read counts 
# for each variety across the four replicates

# use trim_genes function in help_func.R

trimmed_real_data <- trim_genes(counts=real_count, group=real_group, 
                                geneid=real_geneid, z.allow.tr = 2,
                                mean.thr = exp(1))

trimmed_data <- trimmed_real_data[[1]]

saveRDS(trimmed_data, "./real/data/trimmed_real_data.rds")

###########################################################
#Step 3: Estimate parameters based on trimmed real data
###########################################################

trimmed_data <- readRDS("./real/data/trimmed_real_data.rds")
trimmed_cond <- factor(colnames(trimmed_data))
# DGEList object of the trimmed real data
trimmed_object <- trimmed_data %>% DGEList()
# design matrix of trimmed real data
trimmed_design <- model.matrix(~trimmed_cond)
rownames(trimmed_design) <- colnames(trimmed_object)
# dispersion estimation
trimmed_dispsCR <- dispCoxReidInterpolateTagwise(trimmed_object$counts, 
                                                 design=trimmed_design,
                                                 offset=getOffset(trimmed_object), 
                                                 dispersion=.1,
                                                 trend=FALSE, AveLogCPM=NULL, 
                                                 min.row.sum=5,
                                                 prior.df=0, span=0.3, grid.npts=15,
                                                 grid.range=c(-8,8))

#lib size of each sample
trimmed_sample_data = data.frame(trimmed_cond)
trimmed_sample_data$libsize = log(colSums(trimmed_object$counts))
trimmed_libsize = trimmed_sample_data$libsize
# trimmed real data fold-change
nofit = 1000000
trimmed_fc = matrix(nrow=dim(trimmed_object$counts)[1], ncol=2)

for(i in 1:dim(trimmed_object$counts)[1]) {
  f <- negative.binomial(link="log",theta=1/trimmed_dispsCR[i])
  tryCatch({glm(trimmed_object$counts[i,] ~ real_cond + 0, offset=trimmed_libsize, family=f) -> fit},
           warning=function(w) {assign('nofit', c(nofit, i), parent.env(environment()))})
  trimmed_fc[i,] <- fit$coefficients
}

# updated the trimmed real data DGEList object with nofit
trimmed_object <- trimmed_data[-nofit,]%>% DGEList()%>%calcNormFactors()

trimmed_parms <- list(y=trimmed_object, fc=trimmed_fc[-nofit,], 
                      dispsCR = trimmed_dispsCR[-nofit], 
                      sample_data=trimmed_sample_data, 
                      nofit=nofit)

saveRDS(trimmed_parms, "./sim/trimmed_parms.rds")

#####################################################################
#Step 4: Simulate
###############################################################

# load the estimated parameters based on trimmed real data
real_parms <- readRDS("./sim/trimmed_parms.rds")

# Example: simulation scenario 18
#sc18: nGenes=1000, nSample=16, pDiff=0.01, n_sim=5

sc18_alter_pval <- NBsim_scenario(params=real_parms,
                                  nGenes=1000, 
                                  nSample=16, 
                                  pDiff=0.01, nSim=5)

saveRDS(sc18_alter_pval, "./sim/results/sc18_alter_pval.rds")

sc18_alter_pval <- readRDS("./sim/results/sc18_alter_pval.rds")

#add eBayes
##############################################
#sc18
nGenes=1000; nSample=16; pDiff=0.01
diffPerc = pDiff*100; nRep=nSample/2
##############################################
#sim5
i=5;
data_file <- paste(paste("./sim/data/sim_genes",nGenes,"g",nRep, "pDiff",diffPerc, i, sep="_"),"rds",sep=".")

d = as.data.frame(readRDS(data_file))
d_count <- d[,-ncol(d)]
n = ncol(d_count)/2
names(d_count) = c(paste("B73_",1:n, sep=''), paste('Mo17_',1:n, sep=''))
variety = factor(gsub("_[0-9]{1,2}", "", names(d_count)), levels=c("B73","Mo17"))

model = stan_model("./common/normal.stan")
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

saveRDS(analysis, file=paste0("./sim/results/","sim_genes_",nGenes,"_g_",nRep, "_pDiff_",diffPerc, "_",i,"_analysis",
                              ".rds"))

ebayes_analysis <- analysis %>% 
  mutate(p_value = if_else(prob_LDE<prob_HDE, prob_LDE, prob_HDE))
true_de <- d$true_de

ebayes_pval <- as.data.frame(cbind(ebayes_analysis$p_value, true_de))
colnames(ebayes_pval) <- c("ebayes_pval", "de")
saveRDS(ebayes_pval, file=paste0("./sim/results/","sim_genes_",nGenes,"_g_",nRep, "_pDiff_",diffPerc, "_",i,"_ebayes_pval",
                                 ".rds"))

#plot ROC curve for sc18_sim1

alter_pval <- sc18_alter_pval[[i]]
all_pval <- cbind(alter_pval, ebayes_pval$ebayes_pval)
colnames(all_pval)[ncol(all_pval)] <- c("ebayes_pval")

saveRDS(all_pval, file=paste0("./sim/results/","sim_genes_",nGenes,"_g_",nRep, "_pDiff_",diffPerc, "_",i,"_pval",
                              ".rds"))
sc18_sim5_auc <- plot_roc_all(all_result=all_pval, name="sc18 sim5")
saveRDS(sc18_sim5_auc, file=paste0("./sim/results/","sim_genes_",nGenes,"_g_",nRep, "_pDiff_",diffPerc, "_",i,"_auc",
                                   ".rds"))


sc18_auc <- rbind(sc18_sim1_auc, 
                  sc18_sim2_auc,
                  sc18_sim3_auc,
                  sc18_sim4_auc,
                  sc18_sim5_auc)
rownames(sc18_auc) <- paste0("sim",1:5)

saveRDS(sc18_auc, file = paste0("./sim/results/","sim_genes_",nGenes,"_g_",nRep, "_pDiff_",diffPerc,"_auc",
                                ".rds"))













