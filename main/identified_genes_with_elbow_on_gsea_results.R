rm(list = ls())

source("../function/elbow_point_detection.R")

setwd("../data/")
result = readRDS("gsea_results.rds")
thr = seq(50, 500, 50)

identified_genes = list()
gene_phenotypes = c()
r = 1
for(i in names(result)){
  print(r); r = r+1
  df1 = result[[i]][["dp2adr"]]
  df1 = df1[!is.na(df1$pval),]
  df2 = result[[i]][["adr2dp"]]
  df2 = df2[!is.na(df2$pval),]
  
  df1$padj = p.adjust(df1$padj)
  df2$padj = p.adjust(df2$padj)
  
  elbow_index1 = kneedle_method(df1$size, -df1$ES, df1$pval, threshold = 3, n_try = 5)
  elbow_index2 = kneedle_method(df2$size, -df2$ES, df2$pval, threshold = 3)
  
  ########elbow detection
  first_big_jump = function(y){
    dy = y[1:(length(y)-1)] - y[2:length(y)]
    ind = which.max(dy)
    return(ind)
  }
  
  #################
  
  le1 = df1$leadingEdge[elbow_index1][[1]]
  le2 = df2$leadingEdge[elbow_index2][[1]]
  I2 = intersect(le1, le2)
  I3 = union(le1, le2)
  
  L1 = length(le1)
  L2 = length(le2)
  L_intersect = length(I2)
  L_union = length(I3)
  
  if(L_intersect>0){
    if(L_intersect==1){
      print(i)
    }
    
    N_genes = data.frame("phenotype" = i, 
                         "len_dp2adr" = L1, 
                         "len_adr2dp" = L2, 
                         "len_intersect" = L_intersect,
                         "len_union" = L_union)
    
    
    identified_genes[[i]] = list()
    identified_genes[[i]][["adr_side"]] = le1
    identified_genes[[i]][["dp_side"]] = le2
    
    identified_genes[[i]][["intersected_genes"]] = I2
    identified_genes[[i]][["union_genes"]] = I3
    identified_genes[[i]][["N_genes"]] = N_genes
    
    gene_phenotypes = rbind(gene_phenotypes, cbind(rep(i, L_intersect), I2))
  }
}

saveRDS(identified_genes, "identified_genes.rds")


