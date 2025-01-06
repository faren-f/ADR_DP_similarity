rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

library(biomaRt)
library(ReactomePA)

mart = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")


setwd("../data/")
identified_genes = readRDS("identified_genes.rds")


for(i in names(identified_genes)){
  
  genes = identified_genes[[i]][["intersected_genes"]]
  gene_ids = gsub("entrez.", "", genes)
  
  
  reactom_pw = enrichPathway(gene = gene_ids, organism = "human", pvalueCutoff = 0.05,
                             pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10,
                             maxGSSize = 500, readable = FALSE)
  
  if(is.null(reactom_pw) || nrow(reactom_pw@result) == 0 || all(reactom_pw@result$p.adjust > 0.05)){
    next
  }
  
  options(enrichplot.colours = c("#F28BEE", "#27E2C8"))
  print(dotplot(reactom_pw, showCategory = 5) +
          theme(axis.ticks.x = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_line(linewidth=.2, color="lightgray"),
                panel.grid.minor.y = element_line(linewidth=.2, color="lightgray"))+
          ggtitle(i))
  
}









