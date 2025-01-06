rm(list=ls())  
library(ggplot2)

setwd("../data/")
adr_dp = read.csv("adr_dp.csv", row.names = 1)
meddra_hierarchy = readRDS("meddra_hierarchy.rds")
meddra_hierarchy$llt_name = tolower(meddra_hierarchy$llt_name)

colnames(adr_dp) = gsub("\\.", " ", colnames(adr_dp))
I = intersect(rownames(adr_dp), colnames(adr_dp))
adr_dp = adr_dp[I, I]

I = intersect(rownames(adr_dp), meddra_hierarchy$llt_name)
meddra_hierarchy = meddra_hierarchy[meddra_hierarchy$llt_name %in% I,]
all_groups = unique(meddra_hierarchy$soc_name)

paired_SOC = c()
unpaired_SOC = c()

for(i in all_groups){
  llt_i = unique(meddra_hierarchy[meddra_hierarchy$soc_name %in% i, 3])
  if(length(llt_i)>5){
    SOC_i = adr_dp[rownames(adr_dp) %in% llt_i, colnames(adr_dp) %in% llt_i]
    paired_SOC = c(paired_SOC, as.matrix(SOC_i))

    SOC_rest = c(as.matrix(adr_dp[llt_i, !(colnames(adr_dp) %in% llt_i)]), 
                 as.matrix(adr_dp[!(rownames(adr_dp) %in% llt_i), llt_i]))
    
    unpaired_SOC = c(unpaired_SOC, as.matrix(SOC_rest))
  
  }
}

wilcox.test(as.matrix(paired_SOC), as.matrix(unpaired_SOC), alternative = "less")$p.value
