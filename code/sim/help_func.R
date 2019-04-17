#help functions

##########################################################################
#load packages
##########################################################################
library(Paschold2012)

library(MASS)
library(rstan)
library(reshape2)

library(edgeR)
library(DESeq)
library(DESeq2)
library(sSeq)
library(EBSeq)
library(baySeq)
library(ShrinkBayes)


library(tidyr)
library(plyr)
library(dplyr)
library(tibble)

library(ROCR)
library(ggplot2)


###############################################################################
#trim real datasets
###############################################################################

#counts is a num_gene x num_sample matrix
#group is a num_sample length factor mapping columns to groups
#geneid is a num_gene length vector rows to genes, i.e. 1:G or descriptive names
#z.allow.tr is maximum zeros allowable for any treatment(gene)
#mean.thr is a minimum average count for a gene

trim_genes <- function(counts, group, geneid, z.allow.tr = 2, mean.thr = exp(1)){
  
  G = length(geneid)
  key_for_id = cbind(geneid,1:G)
  counts = as.matrix(counts)
  
  #which genes have a treatment with excessive zeros?
  flag_zeros <- NULL
  for(i in group){
    flags  = apply(counts, MARGIN=1,FUN= function(x){
      #indicator for group i (for all genes) 
      sum(x[which(group==i)]==0) > z.allow.tr       
    })
    flag_zeros = c(flag_zeros, which(flags))
  }                  
  
  #Which genes have too low average expression?
  flags2 = apply(counts, MARGIN=1, FUN = function(x) mean(x) < mean.thr)
  flag_low = which(flags2)
  allflags = sort(unique(c(flag_zeros,flag_low)))
  
  trimmed_data <- list(counts[-allflags,],geneid[-allflags])
  return(trimmed_data)
}


edgeR.fit = function(counts, group){
  dge <- DGEList(counts = counts, group = group)
  design <- model.matrix(~ -1 + group)
  dge <- calcNormFactors(dge)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  glmFit(dge, design)
}

edgeR.ord = function(fit){
  apply(fit$coef, 1, function(x){
    if(x["groupB73"] > x["groupMo17"])
      return("LDE")
    else if(x["groupB73"] < x["groupMo17"])
      return("HDE")
    else
      return("nonDE")
  })
}



##########################################################################################
#simulation helper functions
##########################################################################################

#input fitted parameters from preliminary data 
#output NB count simulated data

KS_NBsim <-function(params,nGenes, nSample, pDiff, seed){
    set.seed(seed)
    nRep=nSample/2 #number of replicates in each condition
    fc <-params$fc; dispsCR <- params$dispsCR
    libsizes=params$y@.Data[[2]]$lib.size
    
    #randomly select nGenes with replacement
    id_r <- sample(nrow(fc), nGenes, replace = FALSE) #index of selected 10000 genes
    id_r <- sort(id_r)
    fc_r <- fc[id_r,]
    disps_r <- dispsCR[id_r]
    libsizes_r = runif(nSample, min=min(libsizes), max=max(libsizes))
    
    #de genes index
    indDE <- sample(id_r, floor(pDiff*nGenes), replace=FALSE)
    indDE <- sort(indDE)
    indNonDE<-id_r[!id_r %in% indDE]
    
    true_de <- rep(FALSE, nrow(fc))
    true_de[indDE] <- TRUE
    true_de[indNonDE]<-FALSE
    true_de <- true_de[id_r]
    
    
    #for nonDE genes
    #set two group means the same
    #mu_g1 = mu_g2
    #mu_g1j = libsize_(1j)*exp(fc[g,1])
    #mu_g2j = libsize_(2j)*exp(fc[g,2])
    mu = matrix(nrow=length(disps_r), ncol=nSample)
    
    for (i in 1:length(disps_r)){
      for (j in 1:nSample){
        mu[i,j] <- rnbinom(1, mu = libsizes_r[j]*exp(ifelse(j <= nSample/2, fc_r[i,1],fc_r[i,2])), 
                           size = 1/disps_r[i])
      }
    }
    
    group_mean = matrix(nrow=length(disps_r), ncol=2)
    grp <- seq(1, ncol(mu), by=nRep)
    group_mean<- sapply(grp, function(x) rowSums(mu[ , x:(x+nRep-1)])/nRep)
    
    
    row_expression <- as.matrix((group_mean[,1]+group_mean[,2])/2,ncol=1)
    
    group_mean[!true_de, ] <- c(row_expression[!true_de], row_expression[!true_de])
    
    m = matrix(nrow = length(disps_r), ncol = nSample)
    
    for(i in 1:length(disps_r)) {
      for(j in 1:nSample) {
        m[i,j]=rnbinom(1, mu = ifelse(j <= nSample/2, group_mean[i,1],group_mean[i,2]), 
                       size = 1/disps_r[i])
      }
    }
    
    label <- paste0(c(rep("A_",nRep),rep("B_",nRep)),rep(1:nRep,2))
    
    colnames(m) = label
    rownames(m) = paste0("g",1:length(disps_r))
    
    result <- as.data.frame(cbind(m,true_de))
    return(result)
  }


NBsim_scenario <- function(params,nGenes, nSample, pDiff, nSim){
  nRep = nSample/2
  diffPerc = pDiff*100
  
  
  #auc_final_result <- matrix(nrow=nSim,ncol=5)
  sim_final_result <- list()
  
  for(seed in 1:nSim){
    
    simdata_name <- paste("sim_genes",nGenes,"g",nRep, "pDiff",diffPerc, seed, sep="_")
    
    simdata <- KS_NBsim(params=params,nGenes=nGenes, 
                        nSample=nSample, pDiff=pDiff, seed=seed)
    
    save_filename <- paste(simdata_name, "rds", sep=".")
    
    path <- "../sim/data/"
    
    saveRDS(simdata, paste(path,save_filename,sep=""))
    
    
    sim_count = as.matrix(simdata[,-ncol(simdata)])
    sim_true_de <- as.matrix(simdata[,ncol(simdata)],ncol=1)
    colnames(sim_true_de) <- "true_de"
    sim_split_names <- unlist(strsplit(colnames(sim_count), "[_]"))
    sim_cond =sim_split_names[seq_along(sim_split_names)%%2 !=0]
    
    sim_er <- de_eval(rawdata=sim_count, condition=sim_cond, program="edgeR")
    sim_er_df <- as.data.frame(sim_er[[1]])
    colnames(sim_er_df) <- paste("er", colnames(sim_er_df),sep="_")
    
    sim_ds2 <- de_eval(rawdata=sim_count, condition=sim_cond, program="DESeq2")
    sim_ds2_df <- as.data.frame(sim_ds2[[1]])
    colnames(sim_ds2_df) <- paste("ds2", colnames(sim_ds2_df),sep="_")
    
    sim_ds <- de_eval(rawdata=sim_count, condition=sim_cond, program="DESeq")
    sim_ds_df <- as.data.frame(sim_ds[[1]])
    colnames(sim_ds_df) <- paste("ds", colnames(sim_ds_df),sep="_")
    
    
    sim_ss <- de_eval(rawdata=sim_count, condition=sim_cond, program="sSeq")
    sim_ss_df <- as.data.frame(sim_ss[[1]])
    colnames(sim_ss_df) <- paste("ss", colnames(sim_ss_df),sep="_")
    
    sim_eb <- de_eval(rawdata=sim_count, condition=sim_cond, program="EBSeq")
    sim_eb_df <- as.data.frame(sim_eb[[1]])
    colnames(sim_eb_df) <- paste("eb", colnames(sim_eb_df),sep="_")
    
    
    sim_results <- as.data.frame(cbind(sim_er_df, sim_ds2_df, 
                                           sim_ds_df, sim_ss_df,
                                           sim_eb_df))
    sim_result <- sim_results[, grep('_padj', names(sim_results))]
    sim_result_tb <- cbind(sim_result, sim_true_de)
    sim_result_tb[is.na(sim_result_tb)] <- 1
    
    sim_final_result[[seed]] <- sim_result_tb
    
    sim_result_name <- paste(simdata_name, "pval_result", sep="_")
    result_filename <- paste(sim_result_name, "rds", sep=".")
    saveRDS(sim_result_tb, paste(path,result_filename,sep=""))
    
    
  }
  
  return(sim_final_result)

}


#######################################################################################
#empirical Bayesian de analysis helper functions
#######################################################################################

get_hyperparameters = function(d, variety) {
  
  # phi, alpha, delta parameterization
  ebayes_design = cbind(1,
                 c(1,-1)[as.numeric(variety)])
  
  # GLM fit using edgeR
  fit = d %>% 
    DGEList() %>%
    calcNormFactors %>%
    estimateCommonDisp %>%
    estimateGLMTagwiseDisp(ebayes_design) %>%
    glmFit(ebayes_design)
  
  # Calculate gene-specific estimates for phi, alpha, and psi
  hat = data.frame(gene = 1:length(fit$dispersion),
                   phi   = fit$coefficients[,1] + mean(fit$offset[1,]),
                   alpha = fit$coefficients[,2],
                   psi   = log(fit$dispersion))
  
  ss = hat %>%
    melt(id.vars='gene', variable.name = 'parameter') %>%
    group_by(parameter) %>%
    summarize(mean=mean(value), sd=sd(value))
  
  
  list(eta_phi   = ss$mean[ss$parameter=='phi'],
       eta_alpha = ss$mean[ss$parameter=='alpha'],
       eta_psi   = ss$mean[ss$parameter=='psi'],
       sigma_phi   = ss$sd[ss$parameter=='phi'],
       sigma_alpha = ss$sd[ss$parameter=='alpha'],  
       sigma_psi   = ss$sd[ss$parameter=='psi'],
       c = fit$offset[1,] - mean(fit$offset[1,]))
}



single_gene_analysis = function(counts) {
  diverge = TRUE
  attempt = 1
  
  while (diverge) {
    r = sampling(model, 
                 data = c(list(S       = length(counts), 
                               count   = as.numeric(counts),
                               variety = as.numeric(factor(gsub("_[0-9]{1,2}", "", colnames(counts)), levels=c("B73","Mo17")))),
                          hyperparameters),
                 pars = c("phi","alpha","psi","LDE","HDE"),
                 iter = 2000*2^(attempt-1),
                 thin = 2^(attempt-1))
    
    # Check PSRF for (lack of) convergence
    s = summary(r)$summary
    diverge = any(s[complete.cases(s),"n_eff"] < 1000)
    attempt = attempt + 1
    
  }
  
  alpha_hat = s[rownames(s) == "alpha","mean"] 
  
  
  data.frame(
    phi      = s[rownames(s) == "phi",  "mean"],
    alpha    = s[rownames(s) == "alpha","mean"],
    psi      = s[rownames(s) == "psi",  "mean"],
    prob_LDE = s[rownames(s) == "LDE",  "mean"],
    prob_HDE = s[rownames(s) == "HDE",  "mean"])
}

##########################################################
#eBayes DE analysis
##########################################################
ebayes_de_eval <- function(counts, model, hyperparm){
  
  d_count = counts
  analysis = adply(d_count,
                   1,
                   function(x) single_gene_analysis(x),
                   .id = 'gene',
                   .parallel = parallel,
                   .paropts = list(.export=c('single_gene_analysis','model','hyperparameters'), 
                                   .packages='rstan'))
  
  rownames(analysis) = rownames(d_count)
  
  ebayes_analysis <- analysis %>% 
    mutate(p_value = if_else(prob_LDE<prob_HDE, prob_LDE, prob_HDE))
  
  ebayes_pval <- as.data.frame(ebayes_analysis$p_value, true_de)
  colnames(ebayes_pval) <- c("ebayes_pval")
  
  return(ebayes_pval)
  
}


###########################################################################################
#alternative de analysis methods helper functions
###########################################################################################

de_eval <- function(rawdata, condition, program) {
  pval_list = list()
  
  if(program=="DESeq2") {
    #DESeq2#
    
      dds <- DESeqDataSetFromMatrix(countData = rawdata, colData = data.frame(condition), design = ~condition)
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds)
      res <- results(dds)
      pval = res$pval
      padj = res$padj
      res = cbind(pval, padj)
      pval_list[["ds2"]] <- as.matrix(res)
      rm(dds, res, pval, padj)
      gc()
    
  } else if(program=="DESeq") {
    #DESeq#
    
      DESeq_cds = newCountDataSet(rawdata, condition)
      DESeq_cds = estimateSizeFactors(DESeq_cds)
      DESeq_cds = estimateDispersions(DESeq_cds)
      pval = nbinomTest(DESeq_cds, "A", "B", pvals_only=TRUE)
      padj = p.adjust( pval, method="BH")
      res = cbind(pval, padj)
      pval_list[["ds"]] <- as.matrix(res)
      rm(DESeq_cds, res, pval, padj)
      gc()
    
  } else if(program=="edgeR") {
    #edgeR#
    
      edgeR_cds = DGEList(rawdata, group = condition )
      edgeR_cds = calcNormFactors( edgeR_cds )
      edgeR_cds = estimateCommonDisp( edgeR_cds )
      edgeR_cds = estimateTagwiseDisp( edgeR_cds )
      res = exactTest(edgeR_cds, pair = c("A","B"))$table
      pval = res$PValue
      padj = p.adjust( pval, method="BH")
      res = cbind(pval, padj)
      pval_list[["er"]] <- as.matrix(res)
      rm(edgeR_cds,res, pval, padj)
      gc()
    
  } else if(program=="sSeq") {
    #sSeq#
    
      as.character(condition) -> sSeq_condition
      res <- nbTestSH(rawdata, sSeq_condition, condA = unique(sSeq_condition)[1],condB = unique(sSeq_condition)[2])
      pval = res$pval
      padj = p.adjust(pval, method="BH")
      res = cbind(pval, padj)
      pval_list[["ss"]] <- as.matrix(res)
      rm(res, sSeq_condition, pval, padj)
      gc()
    
  } else if(program=="EBSeq") {
    #EBSeq
    
      Sizes = MedianNorm(rawdata)
      EBOut = EBTest(Data = rawdata, Conditions = condition,sizeFactors = Sizes, maxround = 5)
      data.frame(pval=1-GetPP(EBOut)) -> temp0
      temp1 = rawdata
      merge(temp1, temp0, all.x=TRUE, by.x=0, by.y=0)-> temp2
      pval = temp2[,"pval"]
      names(pval) = temp2[,"Row.names"]
      pval = pval[rownames(rawdata)]
      padj = pval
      res = cbind(pval, padj)
      pval_list[["eb"]] <- as.matrix(res)
      rm(temp0, temp1, temp2, EBOut, Sizes, res, pval, padj)
      gc()
    
  } else {stop("please select a program: DESeq2, DESeq, edgeR, EBSeq or sSeq")}
  pval_list
}


plot_roc_all <- function(all_result, name){
  sim_result_tb <- all_result
  
  sim_ds2_pred <- prediction(sim_result_tb$ds2_padj,sim_result_tb$true_de, label.ordering = c(1,0))
  sim_ds2_perf <- performance(sim_ds2_pred, "tpr","fpr")
  sim_ds2_auc <- performance(sim_ds2_pred, "auc")@y.values[[1]]
  
  sim_ds_pred <- prediction(sim_result_tb$ds_padj,sim_result_tb$true_de, label.ordering = c(1,0))
  sim_ds_perf <- performance(sim_ds_pred, "tpr","fpr")
  sim_ds_auc <- performance(sim_ds_pred, "auc")@y.values[[1]]
  
  sim_er_pred <- prediction(sim_result_tb$er_padj,sim_result_tb$true_de, label.ordering = c(1,0))
  sim_er_perf <- performance(sim_er_pred, "tpr","fpr")
  sim_er_auc <- performance(sim_er_pred, "auc")@y.values[[1]]
  
  sim_ss_pred <- prediction(sim_result_tb$ss_padj,sim_result_tb$true_de, label.ordering = c(1,0))
  sim_ss_perf <- performance(sim_ss_pred, "tpr","fpr")
  sim_ss_auc <- performance(sim_ss_pred, "auc")@y.values[[1]]
  
  
  sim_eb_pred <- prediction(sim_result_tb$eb_padj,sim_result_tb$true_de, label.ordering = c(1,0))
  sim_eb_perf <- performance(sim_eb_pred, "tpr","fpr")
  sim_eb_auc <- performance(sim_eb_pred, "auc")@y.values[[1]]
  
  
  sim_ebayes_pred <- prediction(sim_result_tb$ebayes_pval,sim_result_tb$true_de, label.ordering = c(1,0))
  sim_ebayes_perf <- performance(sim_ebayes_pred, "tpr","fpr")
  sim_ebayes_auc <- performance(sim_ebayes_pred, "auc")@y.values[[1]]
  
  library(ROCR)
  
  plot(sim_ds2_perf, col="blue", main=paste(name,"ROC curves",sep=" "))
  plot(sim_ds_perf, col="orange", add=TRUE)
  plot(sim_er_perf, col="red", add=TRUE)
  plot(sim_ss_perf, col="green",add=TRUE)
  plot(sim_eb_perf, col="brown",add=TRUE)
  plot(sim_ebayes_perf, col="black", add=TRUE)
  
  text <- c(paste("DESeq2 AUC",round(sim_ds2_auc, digits = 2), sep = "="), 
            paste("DESeq AUC",round(sim_ds_auc, digits = 2), sep = "="),
            paste("edgeR AUC",round(sim_er_auc, digits = 2), sep = "="),
            paste("sSeq AUC",round(sim_ss_auc, digits = 2), sep = "="),
            paste("EBSeq AUC",round(sim_eb_auc, digits = 2), sep = "="),
            paste("eBayes AUC", round(sim_ebayes_auc, digits = 2), sep="="))
  col = c("blue","orange","red","green","brown","black")
  lty <- c(1,1,1,1,1,1)
  legend("bottomright",text,text.col=col ,lwd = 1,lty=lty,col= col)
  
  auc_result <- matrix(nrow=1, ncol=6)
  auc_result[1,1]<- sim_ds2_auc
  auc_result[1,2] <- sim_ds_auc
  auc_result[1,3] <- sim_er_auc
  auc_result[1,4] <- sim_ss_auc
  auc_result[1,5] <- sim_eb_auc
  auc_result[1,6] <- sim_ebayes_auc
  colnames(auc_result) <- c("ds2_auc", "ds_auc", "er_auc",
                            "ss_auc", "eb_auc", "ebayes_auc")
  
  return(auc_result)
}









