rm(list=ls()) 
library(fgsea)
library(ggplot2)

setwd("../data/")

# enter benchmark data
benchmark = readRDS("adr_gene_tardis.rds")
drug_se = readRDS("drug_se.rds")[,2:3]
rownames(drug_se) = drug_se$se_name

adr_gene = read.csv("adr_gene.csv", row.names = 1)
dp_gene = read.csv("dp_gene.csv", row.names = 1)

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
df1 = c()
df2 = c()

for(i in adrs){
  print(r); r = r+1
  score_adr_i = t(adr_gene_reduced[i,])[,1]
  score_dp_i = t(dp_gene_reduced[i,])[,1]
  
  gene_scores_adr = data.frame(
    gene = gene_names,
    score_set = score_adr_i
  )
  gene_scores_dp = data.frame(
    gene = gene_names,
    score_set = score_dp_i
  )
  
  gene_scores_adr$rank_set = -rank(gene_scores_adr$score_set)
  gene_scores_dp$rank_set = -rank(gene_scores_dp$score_set)
  
  
  #######
  scores_adr = gene_scores_adr$rank_set
  scores_dp = gene_scores_dp$rank_set
  
  names(scores_adr) = gene_scores_adr$gene
  names(scores_dp) = gene_scores_dp$gene
  
  ######
  genes_benchmark = unique(benchmark[benchmark$meddra_id %in% i, 2])
  benchmark_genes = list()
  benchmark_genes[["benchmark"]] = genes_benchmark
  
  
  fgsea_results1 = fgsea(pathways = benchmark_genes,
                         stats = sort(scores_adr, decreasing = TRUE),
                         minSize = 1,
                         maxSize = 1000,
                         nperm = 10000,
                         scoreType = "pos"
  )
  
  if (nrow(fgsea_results1)){
    df1 = rbind(df1, data.frame(fgsea_results1))
  }else{
    df1 = rbind(df1, data.frame(fgsea_results1)[1,])
  }
  
  
  fgsea_results2 = fgsea(pathways = benchmark_genes,
                         stats = sort(scores_dp, decreasing = TRUE),
                         minSize = 1,
                         maxSize = 500,
                         nperm = 10000,
                         scoreType = "pos"
  )
  if (nrow(fgsea_results2)){
    df2 = rbind(df2, data.frame(fgsea_results2))
  }else{
    df2 = rbind(df2, data.frame(fgsea_results2)[1,])
  }
}

rownames(df1) = adrs
rownames(df2) = adrs

rows_with_na_adr <- apply(df1, 1, function(x) any(is.na(x)))
df1 = df1[!rows_with_na_adr,]

rows_with_na_dp <- apply(df2, 1, function(x) any(is.na(x)))
df2 = df2[!rows_with_na_dp,]


df1$qval = p.adjust(df1$pval, method = "BH")
df2$qval = p.adjust(df2$pval, method = "BH")

saveRDS(df1, "node2vec_adr_gsea_results_with_tardis.rds.rds")
saveRDS(df2, "node2vec_adr_gsea_results_with_tardis.rds.rds")







