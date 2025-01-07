rm(list = ls())
library(igraph)

setwd("../data/")
conversion_table = read.csv("conversion_proteins.csv")
colnames(conversion_table)[1] = "gene_name"
conversion_table$entrez = paste0("entrez.", conversion_table$entrez)
###########
# Read PPI
PPI = read.csv("PPI.csv")
PPI = unique(PPI)
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)
PPI = as_edgelist(PPI_Graph); dim(PPI)
PPI = rbind(PPI, PPI[,2:1]); dim(PPI)
PPI = unique(PPI); dim(PPI)
colnames(PPI) = c("gene1", "gene2")
PPI_genes = V(PPI_Graph)$name; length(PPI_genes)
dim(PPI)

###########
#Read network 
drug_se = readRDS("drug_se.rds")
drug_target = readRDS("drug_target.rds")
disease_symptom = readRDS("disease_symptom.rds")
disease_gene_protein = readRDS("disease_gene.rds")

ADR_DP = unique(disease_symptom[,2:3]); dim(ADR_DP)
ADRs = unique(ADR_DP$meddra_id)
DPs = unique(ADR_DP$hpo_id)
Drugs = unique(drug_se$drugbank_id)
Diseases = unique(disease_gene_protein$mondo_id)
##################################################################################################
# Node Attributes
ADR_Node = cbind(ADRs, c(0:(length(ADRs)-1)))
colnames(ADR_Node) = c("ADR","ADR_node")

Drug_Node = cbind(Drugs, c(0:(length(Drugs)-1)))
colnames(Drug_Node) = c("drug","drug_node")

DP_Node = cbind(DPs, c(0:(length(ADRs)-1)))
colnames(DP_Node) = c("DP","DP_node")

Disease_Node = cbind(Diseases, c(0:(length(Diseases)-1)))
colnames(Disease_Node) = c("disease","disease_node")

gene_Node = cbind(PPI_genes, c(0:(length(PPI_genes)-1)))
colnames(gene_Node) = c("gene","gene_node")

#################################################
#@ Edge index

# ADR-Drug & DP-Disease
edge_index_ADR_drug = c()
edge_index_DP_disease = c()

r = 1
for(i in ADRs){
  print(r); r=r+1
  
  # ADR
  D = unique(drug_se[drug_se$meddra_id == i, 1])
  edge_index_ADR_drug = rbind(edge_index_ADR_drug, cbind(rep(i, length(D)), D))
  
  # DP
  DP = unique(disease_symptom[disease_symptom$meddra_id %in% i, 2])
  Di = unique(disease_symptom[disease_symptom$hpo_id == DP,1])
  edge_index_DP_disease = rbind(edge_index_DP_disease, cbind(rep(DP, length(Di)), Di))
}

colnames(edge_index_ADR_drug) = c("ADR", "drug")
colnames(edge_index_DP_disease) = c("DP", "disease")


# Drug-Gene edge-index
edge_index_drug_gene = c()
r = 1
for(i in Drugs){
  print(r); r=r+1
  # targets of each drug
  targets_Dr = unique(drug_target[drug_target$drugbank_id == i, 2])
  # intersected targets with the genes in the PPI
  targets_Dr = intersect(PPI_genes, targets_Dr)
  
  edge_index_drug_gene = rbind(edge_index_drug_gene, 
                               cbind(rep(i, length(targets_Dr)), targets_Dr))
}

colnames(edge_index_drug_gene) = c("drug", "gene")

# Disease-Gene edge-index
edge_index_disease_gene = c()
r = 0
for(i in Diseases){
  print(r); r=r+1
  # genes associated with each disaese
  gene_Di = unique(disease_gene_protein[disease_gene_protein$mondo_id == i,2])
  gene_Di = intersect(PPI_genes, gene_Di)
  
  edge_index_disease_gene = rbind(edge_index_disease_gene, 
                                  cbind(rep(i, length(gene_Di)), gene_Di))
}
colnames(edge_index_disease_gene) = c("disease", "gene")

################################################################################
# Assign node number to edges
edge_index_ADR_drug = merge(edge_index_ADR_drug, ADR_Node, by = "ADR")
edge_index_ADR_drug = merge(edge_index_ADR_drug, Drug_Node, by = "drug")
edge_index_ADR_drug = edge_index_ADR_drug[, c(2,1,3,4)]

edge_index_DP_disease = merge(edge_index_DP_disease, DP_Node, by = "DP")
edge_index_DP_disease = merge(edge_index_DP_disease, Disease_Node, by = "disease")
edge_index_DP_disease = edge_index_DP_disease[, c(2,1,3,4)]

edge_index_drug_gene = merge(edge_index_drug_gene, gene_Node, by = "gene")
edge_index_drug_gene = merge(edge_index_drug_gene, Drug_Node, by = "drug")
edge_index_drug_gene = edge_index_drug_gene[, c(1,2,4,3)]
edge_index_drug_gene[,2] = paste0("entrez.", edge_index_drug_gene[,2])

edge_index_disease_gene = merge(edge_index_disease_gene, gene_Node, by = "gene")
edge_index_disease_gene = merge(edge_index_disease_gene, Disease_Node, by = "disease")
edge_index_disease_gene = edge_index_disease_gene[, c(1,2,4,3)]
edge_index_disease_gene[,2] = paste0("entrez.", edge_index_disease_gene[,2])

#################################################
#edge_index_PPI

# Merge the two data frames on the 'gene' column
PPI_merged = merge(merge(PPI, gene_Node, by.x = "gene1", by.y = "gene"),
                   gene_Node, by.x = "gene2", by.y = "gene")

edge_index_PPI = PPI_merged[,c(2,1,3,4)]

colnames(edge_index_PPI) = c("gene1", "gene2", "gene_node1", "gene_node2"); dim(edge_index_PPI)
edge_index_PPI[,1] = paste0("entrez.", edge_index_PPI[,1])
edge_index_PPI[,2] = paste0("entrez.", edge_index_PPI[,2])

gene_Node[,1] = paste0("entrez.", gene_Node[,1])


#####Save data
write.table(edge_index_ADR_drug, "edge_index_drug_adr.csv", sep = ",")
write.table(edge_index_DP_disease, "edge_index_disease_dp.csv", sep = ",")
write.table(edge_index_drug_gene, "edge_index_gene_drug.csv", sep = ",")
write.table(edge_index_disease_gene, "edge_index_gene_disease.csv", sep = ",")
write.table(edge_index_PPI, "edge_index_ppi.csv", sep = ",")

write.table(ADR_Node, "adr_node.csv", sep = ",")
write.table(DP_Node, "dp_node.csv", sep = ",")
write.table(Drug_Node, "drug_node.csv", sep = ",")
write.table(Disease_Node, "disease_node.csv", sep = ",")
write.table(gene_Node, "gene_node.csv", sep = ",")
write.table(ADR_DP, "adr_dp.csv", sep = ",")


