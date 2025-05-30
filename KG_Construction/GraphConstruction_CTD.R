rm(list=ls())

setwd("../data/")
drug_se = readRDS("drug_adr.rds")
drug_target = readRDS("drug_target.rds")
disease_symptom = readRDS("disease_dp.rds")
disease_gene_protein = readRDS("disease_gene_CTD.rds")

dis_merge = merge(disease_gene_protein, disease_symptom, by = "mondo_id")
dr_dis_merge = merge(dis_merge, drug_se, by.x = "meddra_id")
dr_dis_merge_all = merge(dr_dis_merge, drug_target, by.x = "drugbank_id")

######### 
drug_se = unique(dr_dis_merge_all[,c(1,2,10,11)])
drug_target = unique(dr_dis_merge_all[,c(1,12,13)])
disease_symptom = unique(dr_dis_merge_all[,c(3,7,2,8,9)])
colnames(disease_symptom)[5] = "disease_name"
disease_gene_protein = unique(dr_dis_merge_all[,c(3,4,5,6)])
colnames(disease_gene_protein)[3] = "disease_name"

saveRDS(drug_se, "drug_se_CTD.rds")
saveRDS(disease_symptom, "disease_symptom_CTD.rds")
saveRDS(disease_gene_protein, "disease_gene_protein_CTD.rds")
saveRDS(drug_target, "drug_target_CTD.rds")






