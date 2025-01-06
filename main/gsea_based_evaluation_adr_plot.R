rm(list=ls()) 
library(fgsea)
library(ggplot2)

setwd("../data/")
benchmark = readRDS("adr_gene_tardis.rds")
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

####################################
i = "meddra.10044074"   

score_adr_i = t(adr_gene_reduced[i,])[,1]
gene_scores_adr = data.frame(
  gene = gene_names,
  score_set = score_adr_i
)

gene_scores_adr$rank_set = -rank(gene_scores_adr$score_set)
scores_adr = gene_scores_adr$rank_set
names(scores_adr) = gene_scores_adr$gene


genes_benchmark = unique(benchmark[benchmark$meddra_id %in% i, 2])
benchmark_genes = list()
benchmark_genes[["benchmark"]] = genes_benchmark

##################
# Plot
enrichment_plot = plotEnrichment(benchmark_genes[["benchmark"]], scores_adr)
plot_data = ggplot_build(enrichment_plot)$data
enrichment_line_data <- plot_data[[1]]
tick_marks_data <- plot_data[[2]]

ggplot() +
  geom_line(
    data = enrichment_line_data,
    aes(x = x, y = y),  
    size = 1.5,
    color = "green"
  ) +
  geom_segment(
    data = tick_marks_data,
    aes(x = x, xend = xend, y = y, yend = yend),  
    color = "black",  
    linewidth = 0.2
  ) +
  labs(
    title = enrichment_plot$labels$title,
    x = enrichment_plot$labels$x,
    y = enrichment_plot$labels$y
  ) +
  theme_minimal()


plotEnrichment(benchmark_genes[["benchmark"]], scores_adr) +
  ggtitle("Zoomed-In Enrichment Plot") +
  xlab("Rank in Ordered Dataset") +
  ylab("Enrichment Score") +
  coord_cartesian(xlim = c(0, 500))+
  geom_segment(
    data = tick_marks_data,
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "black",
    linewidth = 0.2
  )

