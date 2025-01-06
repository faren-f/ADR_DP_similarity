rm(list=ls())

setwd("../data/")
Cortellis_ta = readRDS("Cortellis_ta.rds")
drug_adr = readRDS("drug_adr_withSOC.rds")

# remove cardiac-related drugs
drugs_organ_i_KG = unique(drug_adr[drug_adr$soc_name %in% "Cardiac disorders",4])
drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$heart == 1])

# drugs_organ_i_KG = unique(drug_adr[drug_adr$soc_name %in% "Nervous system disorders",4])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$nervous == 1])

# drugs_organ_i_KG = unique(drug_adr[drug_adr$soc_name %in% "Metabolism and nutrition disorders",4])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$metabolism == 1])

# drugs_organ_i_KG = unique(drug_adr[drug_adr$soc_name %in% "Psychiatric disorders",4])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$mental == 1])

# drugs_organ_i_KG = unique(drug_adr[drug_adr$soc_name %in% "Renal and urinary disorders",4])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$urologic == 1])

# drugs_organ_i_KG = unique(drug_adr[drug_adr$soc_name %in% "Reproductive system and breast disorders",4])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$repro == 1])

# drugs_organ_i_KG = unique(drug_adr[drug_adr$soc_name %in% "Endocrine disorders",4])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$endocrine == 1])

# drugs_organ_i_KG = unique(drug_adr[drug_adr$soc_name %in% "Musculoskeletal and connective tissue disorders",4])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$musculoskeletal == 1])

# drugs_organ_i_KG = unique(drug_adr[drug_adr$soc_name %in% "Hepatobiliary disorders",4])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$hepatobiliary == 1])


# drugs that we should remove from our KG
I = intersect(drugs_organ_i_KG, drugs_With_organ_i_Indication)
drug_adr_organ = drug_adr[!drug_adr$drug_name %in% I,]
drug_adr_organ = drug_adr_organ[,c(2, 1, 3, 4)]
saveRDS(drug_adr_organ, "drug_adr_cardiacIndicatedDrugs_Removed.rds")


