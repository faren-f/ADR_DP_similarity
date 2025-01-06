rm(list=ls()) 
library(ggplot2)

setwd("../data/")

N_tardis = readRDS("Overlap_with_tardis.rds")
N_smit = readRDS("Overlap_with_smit.rds")
N_kuhn = readRDS("Overlap_with_kuhn.rds")
N_curated = readRDS("Overlap_with_curated.rds")
N_literature = readRDS("Overlap_with_pubmed.rds")

df1_node2vec_smit = readRDS("node2vec_adr_gsea_results_with_Smit.rds")
df1_RGCN_smit = readRDS("RGCN_1000run_adr_gsea_results_with_Smit.rds")

node2vec_smit = sum(df1_node2vec_smit$qval<0.05)/N_smit
RGCN_smit = sum(df1_RGCN_smit$qval<0.05)/N_smit

#################################
df1_node2vec_tardis = readRDS("node2vec_adr_gsea_results_with_tardis.rds")
df1_RGCN_tardis = readRDS("RGCN_1000run_adr_gsea_results_with_tardis.rds")

node2vec_tardis = sum(df1_node2vec_tardis$qval<0.05)/N_tardis
RGCN_tardis = sum(df1_RGCN_tardis$qval<0.05)/N_tardis

#################################
df2_node2vec_curated = readRDS("node2vec_dp_gsea_results_with_curated.rds")
df2_RGCN_curated = readRDS("RGCN_1000run_dp_gsea_results_with_curated.rds")

node2vec_curated = sum(df2_node2vec_curated$qval<0.05)/N_curated
RGCN_curated = sum(df2_RGCN_curated$qval<0.05)/N_curated

#################################
df1_node2vec_kuhn = readRDS("node2vec_adr_gsea_results_with_literature_Kuhn.rds")
df1_RGCN_kuhn = readRDS("RGCN_1000run_adr_gsea_results_with_literature_Kuhn.rds")

node2vec_kuhn = sum(df1_node2vec_kuhn$qval<0.05)/N_kuhn
RGCN_kuhn = sum(df1_RGCN_kuhn$qval<0.05)/N_kuhn

#################################
df2_node2vec_literature = readRDS("node2vec_dp_gsea_results_with_literature.rds")
df2_RGCN_literature = readRDS("RGCN_1000run_dp_gsea_results_with_pubmed.rds")

node2vec_literature = sum(df2_node2vec_literature$qval<0.05)/N_literature
RGCN_literature = sum(df2_RGCN_literature$qval<0.05)/N_literature

############################################
values_node2vec = c(node2vec_smit, NA, node2vec_tardis, NA, node2vec_curated, NA, node2vec_literature, NA, node2vec_kuhn)
values_RGCN = c(RGCN_smit, NA, RGCN_tardis, NA, RGCN_curated, NA, RGCN_literature, NA, RGCN_kuhn)

databases = c("Smit", "", "Tardis", "", "Curated", "", "Literature", "", "Kuhn")

df1 = data.frame(methods = "node2vec", databases = databases, values = values_node2vec)
df3 = data.frame(methods = "RGCN", databases = databases, values = values_RGCN)

data = rbind(df1, df3)
data$methods = factor(data$methods, levels = c("node2vec", "RGCN"))
data$databases = factor(data$databases, levels = databases)

color1 = rgb(0.8, 0.2, 0.6, alpha = 0.6) 
color2 = rgb(0.3, 0.8, 0.5, alpha = 0.6)  
colors = c(color1, color2)
  
ggplot(data, aes(x = databases, y = values, fill = methods)) +
  geom_col(position = position_dodge(width = 1), width = 0.8, color = "black") +
  scale_fill_manual(values = colors) +
  theme_classic() +
  scale_y_continuous(n.breaks = 10) +
  labs(x = "", y = "", fill = "Methods") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.ticks.x = element_blank())









