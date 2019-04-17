#prepare real count dataset

real_dataset <- Paschold2012::Paschold2012

de_counts <- real_dataset %>% 
  filter(genotype %in% c("B73","Mo17")) %>% 
  mutate(id=paste(genotype, replicate, sep="_"))%>%
  as.data.frame()%>%
  dplyr::select(GeneID, id, total) %>%
  tidyr::spread(id, total) %>%
  column_to_rownames("GeneID")

paschold <- as.matrix(de_counts) #serve as raw data

saveRDS(paschold, "../real/data/paschold.rds")
