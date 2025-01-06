rm(list=ls()) 
library(ggplot2)

setwd("../data/")
meddra_hierarchy = readRDS("meddra_hierarchy.rds")
disease_symptom = readRDS("disease_symptom.rds")
adr_dp = read.csv("adr_dp.csv", row.names = 1)

colnames(adr_dp) = gsub("\\.", " ", colnames(adr_dp))
I = intersect(rownames(adr_dp), colnames(adr_dp))
adr_dp = adr_dp[I, I]

###########
all_groups = unique(meddra_hierarchy$soc_name)
adr_dp_mat = as.matrix(adr_dp)
unknown_pairs = c(adr_dp_mat[upper.tri(adr_dp_mat)], adr_dp_mat[lower.tri(adr_dp_mat)])

###convert adr_name to meddra_id
meddra_name = unique(disease_symptom[,c(3,4)])
rownames(meddra_name) = meddra_name$symptom_name

I = intersect(rownames(adr_dp), rownames(meddra_name))
meddra_name = meddra_name[I,]
adr_dp = adr_dp[I,I]
rownames(adr_dp) = meddra_name$meddra_id
colnames(adr_dp) = meddra_name$meddra_id

##############################################
boxplot_stat_BG = boxplot.stats(unknown_pairs)
Q1_BG = boxplot_stat_BG$stats[2]
median_BG = boxplot_stat_BG$stats[3]
Q3_BG = boxplot_stat_BG$stats[4]
stats_BG = c("background", Q1_BG, median_BG, Q3_BG)

wilcox_test = c()
Median = c()
Q3 = c()
Q1 = c()
for(i in all_groups){
  llt_i = unique(meddra_hierarchy[meddra_hierarchy$soc_name %in% i, 3])
  
  if(length(llt_i)>5){
    adr_dp_i = adr_dp[rownames(adr_dp) %in% llt_i, colnames(adr_dp) %in% llt_i]
    adr_dp_mat = as.matrix(adr_dp_i)
    known_pairs_i = diag(adr_dp_mat)
    boxplot_stat = boxplot.stats(known_pairs_i)
    Q1 = rbind(Q1, c(i, boxplot_stat$stats[2]))
    Median = rbind(Median, c(i, boxplot_stat$stats[3]))
    Q3 = rbind(Q3, c(i, boxplot_stat$stats[4]))
    
    wilcox_test = rbind(wilcox_test, c(i, wilcox.test(known_pairs_i, unknown_pairs, alternative = "less")$p.value))
  }else{
    next
  }
}

wilcox_test = data.frame(wilcox_test)
wilcox_test$X2 = as.numeric(wilcox_test$X2)
wilcox_test$qval = p.adjust(wilcox_test$X2, method = "BH")
colnames(wilcox_test) = c("group_name", "pval", "qval")

###### median, upper, lower
Q1 = data.frame(Q1)
colnames(Q1) = c("group_name", "Q1")
Q1$Q1 = as.numeric(Q1$Q1)

Median = data.frame(Median)
colnames(Median) = c("group_name", "median")
Median$median = as.numeric(Median$median)

Q3 = data.frame(Q3)
colnames(Q3) = c("group_name", "Q3")
Q3$Q3 = as.numeric(Q3$Q3)


stats_rest = cbind(Q1, Median$median, Q3$Q3, wilcox_test$pval, wilcox_test$qval)
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
  color = c(rep("#F28BEE",21), "#8FAADC")
  
)

data$category = factor(data$category, levels = data$category)

ggplot(data, aes(x = median, y = category)) +
  geom_point(aes(color = I(color)), size = 2) + 
  geom_errorbarh(aes(xmin = Q1, xmax = Q3), height = 0.2, linewidth = 0.3) +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(panel.grid.major.x = element_blank())+
  theme(panel.grid.minor.x = element_blank())+
  geom_text(aes(x = Q3 + 0.1, label = significance), hjust = -0.5) + 
  geom_vline(xintercept = as.numeric(stats_BG[3]), linetype = "dashed", color = "red", linewidth = 0.2)  # Line at OR=1




