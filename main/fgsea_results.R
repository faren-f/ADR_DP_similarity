rm(list = ls())
library(fgsea)

setwd("../data/")
adr_gene = read.table("adr_gene.csv", sep = ",")
dp_gene = read.table("dp_gene.csv", sep = ",")
pathways = gmtPathways("c2.cp.reactome.v2024.1.Hs.entrez.gmt.txt")
meddra_hierarchy = readRDS("meddra_hierarchy.rds")
best_phenos = readRDS("best_phenotypes.rds")
meddra_hierarchy$llt_name = tolower(meddra_hierarchy$llt_name)
gene_names = colnames(adr_gene)
all_sig_phonos = best_phenos[,1]

I = meddra_hierarchy[meddra_hierarchy$llt_code %in% all_sig_phonos, c(1,4,6)]

sig_soc = best_phenos[,3]

I = I[I$soc_name %in% sig_soc,]
adrs = unique(I$llt_name)
####################################
result = list()

for(i in adrs){
  score_adr_i = t(adr_gene[i,])[,1]
  score_dp_i = t(dp_gene[i,])[,1]
  
  gene_scores = data.frame(
    gene = gene_names,
    score_set_1 = score_adr_i,
    score_set_2 = score_dp_i
  )
  
  gene_scores$rank_set_1 = rank(gene_scores$score_set_1)
  gene_scores$rank_set_2 = rank(gene_scores$score_set_2)
  
  
  #######
  scores_adr = gene_scores$rank_set_1
  scores_dp = gene_scores$rank_set_2
  
  names(scores_adr) = gene_scores$gene
  names(scores_dp) = gene_scores$gene
  
  ######
  
  df1 = c()
  df2 = c()
  
  thr = seq(50, 500, 50)

  for (t in thr){
    print(t)
    
    top_genes_adr = gene_scores[gene_scores$rank_set_1 <= t, 1]
    top_genes_dp = gene_scores[gene_scores$rank_set_2 <= t, 1]
    
    top_genes_adr_list = list()
    top_genes_adr_list[["top"]] = top_genes_adr
    
    top_genes_dp_list = list()
    top_genes_dp_list[["top"]] = top_genes_dp
    
    ######
    
    fgsea_results1 = fgsea(pathways = top_genes_dp_list,
                           stats = sort(-scores_adr, decreasing = TRUE),
                           minSize = 1,
                           maxSize = 2000,
                           nperm = 10000,
                           scoreType = "pos"
    )
    
    if (nrow(fgsea_results1)){
      df1 = rbind(df1, data.frame(fgsea_results1))
    }else{
      df1 = rbind(df1, data.frame(fgsea_results1)[1,])
    }
    
    
    fgsea_results2 = fgsea(pathways = top_genes_adr_list,
                           stats = sort(-scores_dp, decreasing = TRUE),
                           minSize = 1,
                           maxSize = 2000,
                           nperm = 10000,
                           scoreType = "pos"
    )
    if (nrow(fgsea_results2)){
      df2 = rbind(df2, data.frame(fgsea_results2))
    }else{
      df2 = rbind(df2, data.frame(fgsea_results2)[1,])
    }
  }
  

  result[[i]] = list()
  result[[i]][["dp2adr"]] = df1
  result[[i]][["adr2dp"]] = df2
}

saveRDS(result, "gsea_results.rds")


