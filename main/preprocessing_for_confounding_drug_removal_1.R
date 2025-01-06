rm(list=ls())

setwd("../data/")
drug_adr = readRDS("drug_se.rds")
Cortellis_ta = readRDS("Cortellis_ta.rds")

I = intersect(drug_adr$drug_name, rownames(Cortellis_ta))
drug_adr = drug_adr[drug_adr$drug_name %in% I,]
Cortellis_ta = Cortellis_ta[I,]

saveRDS(drug_adr, "drug_adr_Cortellis.rds")

meddra_hierarchy = unique(readRDS("meddra_hierarchy.rds")[,3:4])
drug_adr = merge(drug_adr, meddra_hierarchy, by.x = "meddra_id", by.y = "llt_code") 

saveRDS(drug_adr, "drug_adr_withSOC.rds")
