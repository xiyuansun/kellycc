#test script
#set the code folder as the working directory
setwd("./code/")

# getwd()
# [1] "/Users/xiyuansun/Desktop/kellycc/code"

#################################################################################

# load required packages
library(Paschold2012)

library(MASS)
library(rstan)
library(reshape2)

library(edgeR)
library(DESeq)
library(DESeq2)
library(sSeq)
library(EBSeq)

library(tidyr)
library(plyr)
library(dplyr)
library(tibble)

library(ROCR)
library(ggplot2)

#################################################################################
# load the real count data
real_count <- readRDS("./real/data/paschold.rds") 
# retrieve the sample names
real_split_names <- unlist(strsplit(colnames(real_count), "[_]"))
real_sample_names <- real_split_names[seq_along(real_split_names)%%2 !=0]
# convert samples to 2-level factors
real_cond <- factor(real_sample_names)
colnames(real_count) <- real_sample_names

#                    B73  B73  B73  B73 Mo17 Mo17 Mo17 Mo17
# AC147602.5_FG004    1    0    0    2    0    1    0    0
# AC148152.3_FG001    3    4    6    0    8   17   18   20
# AC148152.3_FG005 2323 1533 1932 1945 2070 1582 2196 1882
# AC148152.3_FG006    2    1    0    2    4    4    5    3
# AC148152.3_FG008    3    3    4    1   31   40   45   49
# AC148167.6_FG001  672  598  728  713  743  655  821  824

# group factors
real_group <- real_cond
# [1] B73  B73  B73  B73  Mo17 Mo17 Mo17 Mo17
# Levels: B73 Mo17

# gene id names
real_geneid <- rownames(real_count)
#################################################################################
# test the help function trim_genes
#' Remove the genes with zero counts in all conditions
#' as well as genes whose maximum counts are less then 5

counts=real_count
group = real_group
geneid = real_geneid
z.allow.tr = 2
mean.thr = exp(1)
  
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
trimmed_data <- list(counts[-allflags,],geneid[-allflags])[[1]]

#saveRDS(trimmed_data, "./real/data/trimmed_real_data.rds")

##############################################################################
# parameter estimation based on trimmed count dataset
trimmed_cond <- factor(colnames(trimmed_data))
# [1] B73  B73  B73  B73  Mo17 Mo17 Mo17 Mo17
# Levels: B73 Mo17

#init the DGEList based on trimmed count data
trimmed_object <- trimmed_data %>% 
  DGEList() %>%
  calcNormFactors()

#define the design matrix
trimmed_design <- model.matrix(~trimmed_cond)
rownames(trimmed_design) <- colnames(trimmed_object)
#beta_{g1} = lambda_{g1}
#beta_{g2} = lambda_{g2} - lambda_{g1}

#            (Intercept) trimmed_condMo17
# B73            1                0
# B73            1                0
# B73            1                0
# B73            1                0
# Mo17           1                1
# Mo17           1                1
# Mo17           1                1
# Mo17           1                1

trimmed_cond <- factor(colnames(trimmed_data))
# DGEList object of the trimmed real data
trimmed_object <- trimmed_data %>% DGEList()
trimmed_design <- model.matrix(~trimmed_cond)
rownames(trimmed_design) <- colnames(trimmed_object)

trimmed_dispsCR <- dispCoxReidInterpolateTagwise(trimmed_object$counts, 
                                                 design=trimmed_design, 
                                                 offset=getOffset(trimmed_object), 
                                                 dispersion=.1, trend=FALSE, AveLogCPM=NULL, 
                                                 min.row.sum=5, prior.df=0, span=0.3, 
                                                 grid.npts=15, grid.range=c(-8,8)) 

#normalized lib size of each sample
trimmed_sample_data = data.frame(trimmed_cond)
trimmed_sample_data$libsize = log(colSums(trimmed_object$counts))
trimmed_libsize = trimmed_sample_data$libsize
# 16.26282 16.12499 16.01508 16.08242 16.11941 16.01232 16.14484 16.08049
nofit = 1000000
trimmed_fc = matrix(nrow=dim(trimmed_object$counts)[1], ncol=2)

for(i in 1:dim(trimmed_object$counts)[1]) {
  f <- negative.binomial(link="log",theta=1/trimmed_dispsCR[i])
  tryCatch({glm(trimmed_object$counts[i,] ~ trimmed_cond + 0, offset=trimmed_libsize, family=f) -> fit},
           warning=function(w) {assign('nofit', c(nofit, i), parent.env(environment()))})
  trimmed_fc[i,] <- fit$coefficients
}

trimmed_object <- DGEList(counts = trimmed_data[-nofit,])

trimmed_parms <- list(y=trimmed_object, fc=trimmed_fc[-nofit,], 
                      dispsCR = trimmed_dispsCR[-nofit], 
                      sample_data=trimmed_sample_data, 
                      nofit=nofit)

#saveRDS(trimmed_parms, "./real/data/trimmed_parms.rds")




