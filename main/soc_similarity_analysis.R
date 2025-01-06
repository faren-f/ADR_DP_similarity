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

################################################################################
unpaired_SOC = c()

for(i in all_groups){
  llt_i = unique(meddra_hierarchy[meddra_hierarchy$soc_name %in% i, 3])
  if(length(llt_i)>5){
    SOC_rest = c(as.matrix(adr_dp[llt_i, !(colnames(adr_dp) %in% llt_i)]), 
                 as.matrix(adr_dp[!(rownames(adr_dp) %in% llt_i), llt_i]))
    
    unpaired_SOC = c(unpaired_SOC, as.matrix(SOC_rest))
  }
}
#################################################################################################
boxplot_stat_BG = boxplot.stats(unpaired_SOC)
Q1_BG = boxplot_stat_BG$stats[2]
median_BG = boxplot_stat_BG$stats[3]
Q3_BG = boxplot_stat_BG$stats[4]
stats_BG = c("background", Q1_BG, median_BG, Q3_BG)

wilcox_test_organ = c()
median_groups = c()
Q3_groups = c()
Q1_groups = c()

r = 1
for(i in all_groups){
  llt_i = unique(meddra_hierarchy[meddra_hierarchy$soc_name %in% i, 3])

  if(length(llt_i)>5){
    SOC_i = adr_dp[rownames(adr_dp) %in% llt_i, colnames(adr_dp) %in% llt_i]
    SOC_i = as.matrix(SOC_i)
    
    boxplot_stat = boxplot.stats(SOC_i)
    Q1_groups = rbind(Q1_groups, c(i, boxplot_stat$stats[2]))
    median_groups = rbind(median_groups, c(i, boxplot_stat$stats[3]))
    Q3_groups = rbind(Q3_groups, c(i, boxplot_stat$stats[4]))
    
    numbers_assigned_to_group_labels = c(numbers_assigned_to_group_labels, rep(r+1, length(SOC_i)))
    r = r+2
    
    wilcox_test_organ = rbind(wilcox_test_organ, c(i, wilcox.test(SOC_i, unpaired_SOC, alternative = "less")$p.value))
  }
}

wilcox_test_organ = data.frame(wilcox_test_organ)
colnames(wilcox_test_organ) = c("SOC", "pval")
wilcox_test_organ$pval = as.numeric(wilcox_test_organ$pval)
wilcox_test_organ$qval = p.adjust(wilcox_test_organ$pval, method = "BH")
colnames(wilcox_test_organ) = c("group_name", "pval", "qval")

#####################################
Q1_groups = data.frame(Q1_groups)
colnames(Q1_groups) = c("group_name", "Q1")
Q1_groups$Q1 = as.numeric(Q1_groups$Q1)

median_groups = data.frame(median_groups)
colnames(median_groups) = c("group_name", "median")
median_groups$median = as.numeric(median_groups$median)

Q3_groups = data.frame(Q3_groups)
colnames(Q3_groups) = c("group_name", "Q3")
Q3_groups$Q3 = as.numeric(Q3_groups$Q3)

stats_rest = cbind(Q1_groups, median_groups$median, Q3_groups$Q3, wilcox_test_organ$pval, wilcox_test_organ$qval)
colnames(stats_rest) = c("group", "Q1", "median", "Q3", "pval", "adj_pval")

###sort by median
sorted_data_stat = stats_rest[order(stats_rest$median, decreasing = TRUE), ]

stats_BG = c(stats_BG, 1, 1)

sorted_data_stat = rbind(sorted_data_stat, stats_BG)
sorted_data_stat$Q1 = as.numeric(sorted_data_stat$Q1)
sorted_data_stat$median = as.numeric(sorted_data_stat$median)
sorted_data_stat$Q3 = as.numeric(sorted_data_stat$Q3)
sorted_data_stat$pval = as.numeric(sorted_data_stat$pval)
sorted_data_stat$adj_pval = as.numeric(sorted_data_stat$adj_pval)

rownames(sorted_data_stat) = sorted_data_stat$group
significance = c()
for(i in rownames(sorted_data_stat)){
  
  if(sorted_data_stat[i,5] >0.05){
    significance = c(significance, "")
    
  }else if((sorted_data_stat[i,5] <0.05) & (sorted_data_stat[i,5] >0.01)){
    significance = c(significance, "*")
    
  }else if((sorted_data_stat[i,5] <0.01) & (sorted_data_stat[i,5] >0.001)){
    significance = c(significance, "**")
    
  }else if(sorted_data_stat[i,5] <0.001){
    significance = c(significance, "***")
    
  }
}                                                                       

#################plot
data = data.frame(
  category = sorted_data_stat$group,
  median = sorted_data_stat$median,
  Q1 = sorted_data_stat$Q1,
  Q3 = sorted_data_stat$Q3,
  significance = significance,
  color = c(rep("#27E2C8",21), "#8FAADC")
  
)

data$category = factor(data$category, levels = data$category)

p = ggplot(data, aes(x = median, y = category)) +
  geom_point(aes(color = I(color)), size = 2) + 
  geom_errorbarh(aes(xmin = Q1, xmax = Q3), height = 0.2, linewidth = 0.3) +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(panel.grid.major.x = element_blank())+
  theme(panel.grid.minor.x = element_blank())+
  geom_text(aes(x = Q3 + 0.1, label = significance), hjust = -0.5)+
  geom_vline(xintercept = as.numeric(stats_BG[3]), linetype = "dashed", color = "red", linewidth = 0.2)  # Line at OR=1
  


