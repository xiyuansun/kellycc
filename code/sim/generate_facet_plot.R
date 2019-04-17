#prepare result facet plot
setwd("~/Desktop/kellycc/code/sim")
library(dplyr)
library(plyr)
library(reshape2)
library(stringr)
library(ggplot2)



nsample_vec <- c(4,8,16)
pdiff_vec <- c(1,10,30)
ngenes_vec <- c(1000,10000)

auc_result <- c()


for(i in 1:3){
  for(j in 1:3){
    for(k in 1:2){
      nSample = nsample_vec[i]
      pDiff = pdiff_vec[j]
      nGenes = ngenes_vec[k]
      
      auc_res <- readRDS(paste0(getwd(),
                                "/results/sim_genes_",nGenes,
                                "_g_",nSample/2,
                                "_pDiff_",pDiff,"_auc.rds"))
      
      auc_long <- auc_res %>%
        melt()%>%
        mutate(gene_num = nGenes,
               pdiff_val = pDiff,
               sample_num = nSample,
               sim = str_sub(Var1, start= -1),
               method = sub("_auc","",Var2),
               auc_value = value
               ) %>%
        select(-c(Var1, Var2, value))
      
      
      auc_result <- rbind(auc_result, auc_long)
      
      
      
    }
  }
}



saveRDS(auc_result, "~/Desktop/kellycc/code/sim/results/auc_results.rds")

#generate plots
auc_result <- readRDS("~/Desktop/kellycc/code/sim/results/auc_results.rds")
head(auc_result)
auc_table <- auc_result %>%
  mutate(pdiff_fac = paste0("pDiff = ",pdiff_val,"%"),
         sample_fac = paste0("nSample = ", sample_num),
         genes = as.factor(gene_num))

auc_table$ord_sample_fac <- factor(auc_table$sample_fac)

auc_table$ord_sample_fac<- ordered(auc_table$ord_sample_fac, 
                                    levels = c("nSample = 4", "nSample = 8", "nSample = 16"))


sp <- ggplot(auc_table, 
             aes(x=reorder(method,auc_value,mean), 
                 y=auc_value, group=genes)) + 
  geom_point(aes(shape=genes, color=genes)) + 
  scale_shape_manual(values=c(1,8))+
  scale_color_manual(values=c("red", "blue"))+
  xlab("Method")+ylab("AUC Value")+theme_bw()

sp

# Divide by levels of "nSample", in the vertical direction
sp + facet_grid(ord_sample_fac ~ pdiff_fac)










