rm(list=ls()) 
library(ggplot2)

setwd("../data/")
data_SOC_drugs = readRDS("average_drugs_in_each_SOC.rds")
data_SOC_diseases = readRDS("average_diseases_in_each_SOC.rds")
wilcox_test_organ = readRDS("wilkox_in_each_SOC.rds")

SOC_data = merge(wilcox_test_organ, data_SOC_drugs, by.x = "group_name", by.y = "SOC")
SOC_data = merge(SOC_data, data_SOC_diseases, by.x = "group_name", by.y = "SOC")
SOC_data$ave_DrugsDiseases = rowMeans(as.matrix(cbind(SOC_data$mean_drugs, SOC_data$mean_diseases)))

wilcox.test(SOC_data$N_adrs[SOC_data$pval<0.05], SOC_data$N_adrs[SOC_data$pval>0.05], alternative = "two.sided")
wilcox.test(SOC_data$mean_drugs[SOC_data$pval<0.05], SOC_data$mean_drugs[SOC_data$pval>0.05], alternative = "two.sided")
wilcox.test(SOC_data$mean_diseases[SOC_data$pval<0.05], SOC_data$mean_diseases[SOC_data$pval>0.05], alternative = "two.sided")
wilcox.test(SOC_data$ave_DrugsDiseases[SOC_data$pval<0.05], SOC_data$ave_DrugsDiseases[SOC_data$pval>0.05], alternative = "greater")

######## Plots
SOC_data$significance = ifelse(SOC_data$pval < 0.05, "pval < 0.05", "pval > 0.05")

p = ggplot(SOC_data, aes_string(x = "significance", y = "N_adrs")) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "darkblue", size = 3) + # Adjust the size here
  theme(
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)# Remove minor grid lines
  )

print(p)

########
p = ggplot(SOC_data, aes_string(x = "significance", y = "mean_drugs")) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "darkblue", size = 3) + # Adjust the size here
  theme(
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)# Remove minor grid lines
  )

print(p)

########
p = ggplot(SOC_data, aes_string(x = "significance", y = "mean_diseases")) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "darkblue", size = 3) + # Adjust the size here
  theme(
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)# Remove minor grid lines
  )

print(p)

########
p = ggplot(SOC_data, aes_string(x = "significance", y = "ave_DrugsDiseases")) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "darkblue", size = 3) + # Adjust the size here
  theme(
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)# Remove minor grid lines
  )
print(p)






