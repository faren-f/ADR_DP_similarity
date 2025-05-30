rm(list = ls())
library(igraph)

#Read network 
setwd("../data/")
drug_se = readRDS("drug_adr_PT_KG.rds")
drug_target = readRDS("drug_target_PT_KG.rds")
disease_symptom = readRDS("disease_dp_PT_KG.rds")
disease_gene_protein = readRDS("disease_gene_PT_KG.rds")

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

ADR_DP = unique(disease_symptom[,2:3]); dim(ADR_DP)
ADRs = unique(ADR_DP$meddra_id)
DPs = unique(ADR_DP$hpo_id)
Drugs = unique(drug_se$drugbank_id)
Diseases = unique(disease_gene_protein$mondo_id)

#################################################
# Node Attributes
L1 = length(ADRs); print(L1)
ADR_Node = cbind(ADRs, c(0:(L1-1))); dim(ADR_Node)
colnames(ADR_Node) = c("ADR","ADR_node")

L2 = L1 + length(DPs); print(L2)
DP_Node = cbind(DPs, c(L1:(L2-1))); dim(DP_Node)
colnames(DP_Node) = c("DP","DP_node")

L3 = L2 + length(Drugs); print(L3)

Drug_Node = cbind(Drugs, c(L2:(L3-1))); dim(Drug_Node)
colnames(Drug_Node) = c("drug","drug_node")

L4 = L3 + length(Diseases); print(L4)

Disease_Node = cbind(Diseases, c(L3:(L4-1)))
colnames(Disease_Node) = c("disease","disease_node")

L5 = L4 + length(PPI_genes); print(L5)

gene_Node = cbind(PPI_genes, c(L4:(L5-1))); dim(gene_Node)
colnames(gene_Node) = c("gene","gene_Node")

#################################################
# Edge index
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
  targets_Dr = unique(drug_target[drug_target$drugbank_id == i, 2])
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


# Assign node number to edges
edge_index_ADR_drug = merge(edge_index_ADR_drug, ADR_Node, by = "ADR")
edge_index_ADR_drug = merge(edge_index_ADR_drug, Drug_Node, by = "drug")
edge_index_ADR_drug = edge_index_ADR_drug[, c(4,3)]
edge_index_ADR_drug = as.matrix(edge_index_ADR_drug)
edge_index_ADR_drug = rbind(edge_index_ADR_drug, edge_index_ADR_drug[,2:1])


edge_index_DP_disease = merge(edge_index_DP_disease, DP_Node, by = "DP")
edge_index_DP_disease = merge(edge_index_DP_disease, Disease_Node, by = "disease")
edge_index_DP_disease = edge_index_DP_disease[, c(4,3)]
edge_index_DP_disease = as.matrix(edge_index_DP_disease)
edge_index_DP_disease = rbind(edge_index_DP_disease, edge_index_DP_disease[,2:1])


edge_index_drug_gene = merge(edge_index_drug_gene, gene_Node, by = "gene")
edge_index_drug_gene = merge(edge_index_drug_gene, Drug_Node, by = "drug")
edge_index_drug_gene = edge_index_drug_gene[, c(4,3)]
edge_index_drug_gene = as.matrix(edge_index_drug_gene)
edge_index_drug_gene = rbind(edge_index_drug_gene, edge_index_drug_gene[,2:1])


edge_index_disease_gene = merge(edge_index_disease_gene, gene_Node, by = "gene")
edge_index_disease_gene = merge(edge_index_disease_gene, Disease_Node, by = "disease")
edge_index_disease_gene = edge_index_disease_gene[, c(4,3)]
edge_index_disease_gene = as.matrix(edge_index_disease_gene)
edge_index_disease_gene = rbind(edge_index_disease_gene, edge_index_disease_gene[,2:1])

#################################################
#edge_index_PPI
# Merge the two data frames on the 'gene' column
PPI_merged = merge(merge(PPI, gene_Node, by.x = "gene1", by.y = "gene"),
                   gene_Node, by.x = "gene2", by.y = "gene")

edge_index_PPI = PPI_merged[,c(3,4)]

colnames(edge_index_PPI) = c("gene_node1", "gene_node2"); dim(edge_index_PPI)

######combine all the edge indexes
colnames(edge_index_ADR_drug) = c("r1", "r2")
colnames(edge_index_DP_disease) = c("r1", "r2")
colnames(edge_index_drug_gene) = c("r1", "r2")
colnames(edge_index_disease_gene) = c("r1", "r2")
colnames(edge_index_PPI) = c("r1", "r2")
edge_index_homo = rbind(edge_index_ADR_drug, edge_index_DP_disease, 
                        edge_index_drug_gene, edge_index_disease_gene,
                        edge_index_PPI)


gene_Node[,1] = paste0("entrez.", gene_Node[,1])
edge_index_homo[,1] = as.numeric(edge_index_homo[,1])
edge_index_homo[,2] = as.numeric(edge_index_homo[,2])

####add names to nodes
ADR_names = unique(disease_symptom[disease_symptom$meddra_id %in% ADR_Node[,1],3:4])
ADR_Node = merge(ADR_Node, ADR_names, by.x = 'ADR', by.y = 'meddra_id')
ADR_Node = ADR_Node[order(as.numeric(ADR_Node[,2]), decreasing = FALSE),]

DP_names = unique(disease_symptom[disease_symptom$hpo_id %in% DP_Node[,1],c(2,4)])
DP_Node = merge(DP_Node, DP_names, by.x = 'DP', by.y = 'hpo_id')
DP_Node = DP_Node[order(as.numeric(DP_Node[,2]), decreasing = FALSE),]

Drug_names = unique(drug_se[drug_se$drugbank_id %in% Drug_Node[,1],c(1,4)])
Drug_Node = merge(Drug_Node, Drug_names, by.x = 'drug', by.y = 'drugbank_id')
Drug_Node = Drug_Node[order(as.numeric(Drug_Node[,2]), decreasing = FALSE),]

Disease_names = unique(disease_symptom[disease_symptom$mondo_id %in% Disease_Node[,1],c(1,5)])
Disease_Node = merge(Disease_Node, Disease_names, by.x = 'disease', by.y = 'mondo_id')
Disease_Node = Disease_Node[order(as.numeric(Disease_Node[,2]), decreasing = FALSE),]

gene_Node = cbind(gene_Node, gene_Node[,1])


colnames(ADR_Node) = c('c1', 'c2', 'c3')
colnames(DP_Node) = c('c1', 'c2', 'c3')
colnames(Drug_Node) = c('c1', 'c2', 'c3')
colnames(Disease_Node) = c('c1', 'c2', 'c3')
colnames(gene_Node) = c('c1', 'c2', 'c3')

Conversion_table = rbind(ADR_Node, DP_Node, Drug_Node, Disease_Node, gene_Node)
colnames(Conversion_table) = c("Names", "Nodes", "Symbole")

#####Save data
write.table(edge_index_homo, "edge_index_homogenious_PT.csv", sep = ",")
write.table(Conversion_table, "conversion_table_homogenious_PT.csv", sep = ",")


