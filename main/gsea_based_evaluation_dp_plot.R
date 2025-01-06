rm(list=ls()) 
library(fgsea)
library(ggplot2)


setwd("../data/")
benchmark = readRDS("symptom_gene_curated.rds")
adr_gene = read.csv("adr_gene.csv", row.names = 1)
dp_gene = read.csv("dp_gene.csv", row.names = 1)
drug_se = unique(readRDS("drug_se.rds")[,2:3])

rownames(drug_se) = drug_se$se_name
I = intersect(rownames(adr_gene), rownames(dp_gene))
adr_gene = adr_gene[I, ]
dp_gene = dp_gene[I, ]
identical(rownames(adr_gene), rownames(dp_gene))

I = intersect(drug_se$se_name, rownames(adr_gene))
drug_se = drug_se[I,]
adr_gene = adr_gene[I,]
dp_gene = dp_gene[I,]

rownames(adr_gene) = drug_se$meddra_id
rownames(dp_gene) = drug_se$meddra_id

###################
adrs = intersect(benchmark$meddra_id, rownames(adr_gene))
gene_names = colnames(adr_gene)

adr_gene_reduced = adr_gene[adrs,]
dp_gene_reduced = dp_gene[adrs,]
benchmark = unique(benchmark[benchmark$meddra_id %in% adrs,])

###################################################################################
i = "meddra.10010071"   
score_dp_i = t(dp_gene_reduced[i,])[,1]

gene_scores_dp = data.frame(
  gene = gene_names,
  score_set = score_dp_i
)

gene_scores_dp$rank_set = -rank(gene_scores_dp$score_set)
scores_dp = gene_scores_dp$rank_set
names(scores_dp) = gene_scores_dp$gene

######
genes_benchmark = unique(benchmark[benchmark$meddra_id %in% i, 2])
benchmark_genes = list()
benchmark_genes[["benchmark"]] = genes_benchmark

##########################################
# enrichment score plot

enrichment_plot <- plotEnrichment(benchmark_genes[["benchmark"]], scores_dp)

# Extract data from the ggplot object
plot_data <- ggplot_build(enrichment_plot)$data

# Extract the main enrichment line data (smooth curve)
enrichment_line_data <- plot_data[[1]]

# Extract the black tick marks data (ranked list indicators)
tick_marks_data <- plot_data[[2]]

# Inspect column names of tick_marks_data
head(tick_marks_data)  # Verify column names

# Recreate the plot
ggplot() +
  geom_segment(
    data = tick_marks_data,
    aes(x = x, xend = x, y = y, yend = yend),  # Adjust based on column names
    color = "black"
  ) +
  geom_line(
    data = enrichment_line_data,
    aes(x = x, y = y),  # Green enrichment line
    size = 1.5,
    color = "green"
  ) +
  labs(
    title = enrichment_plot$labels$title,
    x = enrichment_plot$labels$x,
    y = enrichment_plot$labels$y
  ) +
  theme_minimal()
