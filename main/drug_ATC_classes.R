rm(list=ls()) 
library(stringr)
setwd("../data/")

drug_adr = readRDS("drug_adr.rds")

Cortellis_ta = readRDS("Cortellis_ta.rds")
Cortellis_ta2 = Cortellis_ta
Cortellis_ta2$name = rownames(Cortellis_ta)

drug_classes_pubchem = read.csv("PubChem_ATC_Code.csv")
drug_classes_pubchem$cmpdname = tolower(drug_classes_pubchem$cmpdname)
drug_classes_pubchem = unique(drug_classes_pubchem[,c(2,3,40)])

drug_class_all = read.csv("drugbank_annotation.csv")
drug_class_all = unique(drug_class_all)
##########################
organs = c("Cardiac disorders", "Metabolism and nutrition disorders", "Psychiatric disorders", 
          "Renal and urinary disorders", "Reproductive system and breast disorders", "Hepatobiliary disorders")

# change this line for different organs
organ = organs[6]
drugs_organ_i_KG = unique(drug_adr[drug_adr$soc_name %in% organ,4])

drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$heart == 1])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$nervous == 1])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$metabolism == 1])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$mental == 1])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$urologic == 1])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$repro == 1])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$endocrine == 1])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$musculoskeletal == 1])
# drugs_With_organ_i_Indication = unique(rownames(Cortellis_ta)[Cortellis_ta$hepatobiliary == 1])


I = intersect(drugs_organ_i_KG, drugs_With_organ_i_Indication)
drug_class_adr_indication_nervous = unique(drug_class_all[drug_class_all$drug_name %in% I,]) 
drug_adr_removed = drug_adr[!drug_adr$drug_name %in% I,]


drug_adr_organ = drug_adr_removed[drug_adr_removed$soc_name %in% organ,]
drug_adr_rest = drug_adr_removed[!(drug_adr_removed$soc_name %in% organ),]

drug_adr_organ = drug_adr_organ[,c(2, 1, 3, 4, 5)]
drug_adr_rest = drug_adr_rest[,c(2, 1, 3, 4, 5)]

drug_adr_organ$drugbank_id = gsub("drugbank.", "", drug_adr_organ$drugbank_id)
drug_adr_rest$drugbank_id = gsub("drugbank.", "", drug_adr_rest$drugbank_id)


############# drug classes organ
all_drugs = unique(drug_adr_organ$drugbank_id)


drug_class_organ = c()
r = 1
for(i in all_drugs){
  print(r); r = r+1
  flag = unlist(lapply(drug_classes_pubchem$cmpdsynonym, function(x){grepl(i, x)[[1]]}))
  if(sum(flag)>=1){
    drug_class_organ = rbind(drug_class_organ, cbind(rep(i, sum(flag)), drug_classes_pubchem$annotation[flag]))
  }else{
    drug_name = unique(drug_adr_organ[drug_adr_organ$drugbank_id %in% i,4])
    flag2 = drug_classes_pubchem$cmpdname %in% drug_name
    flag3 = unlist(lapply(drug_classes_pubchem$cmpdsynonym, function(x){grepl(paste0("\\b", drug_name, "\\b"), x, ignore.case = TRUE)[[1]]}))
    
    if(sum(flag2)>=1){
      ann = drug_classes_pubchem[drug_classes_pubchem$cmpdname %in% drug_name, 3]
      drug_class_organ = rbind(drug_class_organ, cbind(rep(i, sum(flag2)), ann))
    }
    else if(sum(flag3)>=1){
      drug_class_organ = rbind(drug_class_organ, cbind(rep(i, sum(flag3)), drug_classes_pubchem$annotation[flag3]))
    }else{
      drug_class_organ = rbind(drug_class_organ, c(i, i))
      
    }
  }
}

drug_class_organ = data.frame(drug_class_organ)
drug_class_organ = unique(drug_class_organ)

colnames(drug_class_organ) = c("drugbank_id", "annotation")

drug_class_organ = merge(drug_adr_organ[,c(1,4)], drug_class_organ, by = "drugbank_id")
drug_class_organ = unique(drug_class_organ)

######extract ATC pattern
atc_pattern <- "\\b[A-Z]\\d{2}([A-Z]*)?( > \\b[A-Z]\\d{2}([A-Z]*)?)*\\b"

extract_atc <- function(text) {
  
  entries <- unlist(strsplit(text, "\\|")) 
  atc_entries <- grep(atc_pattern, entries, value = TRUE) 
  
  return(atc_entries)
}

drug_class_organ$atc_hierarchy_list <- lapply(drug_class_organ$annotation, function(x){extract_atc(x[[1]])})

atc_hierarchy_organ <- data.frame(
  original_row = rep(seq_along(drug_class_organ$annotation), lengths(drug_class_organ$atc_hierarchy_list)),
  drugbank_id = rep(drug_class_organ$drugbank_id, lengths(drug_class_organ$atc_hierarchy_list)),
  drug_name = rep(drug_class_organ$drug_name, lengths(drug_class_organ$atc_hierarchy_list)),
  atc_hierarchy = unlist(drug_class_organ$atc_hierarchy_list),
  annotiation = rep(drug_class_organ$annotation, lengths(drug_class_organ$atc_hierarchy_list)),
  stringsAsFactors = FALSE
)

##############seperate ATC levels

split_hierarchy <- function(hierarchy) {
  levels <- strsplit(hierarchy, " > ")[[1]]
  
  first_level <- ifelse(length(levels) >= 1, levels[1], NA)
  second_level <- ifelse(length(levels) >= 2, levels[2], NA)
  third_level <- ifelse(length(levels) >= 3, levels[3], NA)
  fourth_level <- ifelse(length(levels) >= 4, levels[4], NA)
  
  return(data.frame(
    First_Level = first_level,
    Second_Level = second_level,
    Third_Level = third_level,
    Fourth_Level = fourth_level,
    stringsAsFactors = FALSE
  ))
}

atc_levels <- do.call(rbind, lapply(atc_hierarchy_organ$atc_hierarchy, split_hierarchy))
atc_hierarchy_organ = cbind(atc_hierarchy_organ, atc_levels)



################################### drug classes rest
all_drugs = unique(drug_adr_rest$drugbank_id)

drug_class_rest = c()
r = 1
for(i in all_drugs){
  print(r); r = r+1
  flag = unlist(lapply(drug_classes_pubchem$cmpdsynonym, function(x){grepl(i, x)[[1]]}))
  if(sum(flag)>=1){
    drug_class_rest = rbind(drug_class_rest, cbind(rep(i, sum(flag)), drug_classes_pubchem$annotation[flag]))
  }else{
    drug_name = unique(drug_adr_rest[drug_adr_rest$drugbank_id %in% i,4])
    flag2 = drug_classes_pubchem$cmpdname %in% drug_name
    flag3 = unlist(lapply(drug_classes_pubchem$cmpdsynonym, function(x){grepl(paste0("\\b", drug_name, "\\b"), x, ignore.case = TRUE)[[1]]}))
    
    if(sum(flag2)>=1){
      ann = drug_classes_pubchem[drug_classes_pubchem$cmpdname %in% drug_name, 3]
      drug_class_rest = rbind(drug_class_rest, cbind(rep(i, sum(flag2)), ann))
    }
    else if(sum(flag3)>=1){
      drug_class_rest = rbind(drug_class_rest, cbind(rep(i, sum(flag3)), drug_classes_pubchem$annotation[flag3]))
    }else{
      drug_class_rest = rbind(drug_class_rest, c(i, i))
      
    }
  }
}

drug_class_rest = data.frame(drug_class_rest)
drug_class_rest = unique(drug_class_rest)

colnames(drug_class_rest) = c("drugbank_id", "annotation")

drug_class_rest = merge(drug_adr_rest[,c(1,4)], drug_class_rest, by = "drugbank_id")
drug_class_rest = unique(drug_class_rest)

######extract ATC pattern
atc_pattern <- "\\b[A-Z]\\d{2}([A-Z]*)?( > \\b[A-Z]\\d{2}([A-Z]*)?)*\\b"

extract_atc <- function(text) {
  
  entries <- unlist(strsplit(text, "\\|")) # Split the text by '|' separator
  atc_entries <- grep(atc_pattern, entries, value = TRUE) # Find entries that match ATC codes and their descriptions
  
  return(atc_entries)
}

drug_class_rest$atc_hierarchy_list <- lapply(drug_class_rest$annotation, function(x){extract_atc(x[[1]])})

atc_hierarchy_rest <- data.frame(
  original_row = rep(seq_along(drug_class_rest$annotation), lengths(drug_class_rest$atc_hierarchy_list)),
  drugbank_id = rep(drug_class_rest$drugbank_id, lengths(drug_class_rest$atc_hierarchy_list)),
  drug_name = rep(drug_class_rest$drug_name, lengths(drug_class_rest$atc_hierarchy_list)),
  atc_hierarchy = unlist(drug_class_rest$atc_hierarchy_list),
  annotiation = rep(drug_class_rest$annotation, lengths(drug_class_rest$atc_hierarchy_list)),
  stringsAsFactors = FALSE
)


##############seperate ATC levels

split_hierarchy <- function(hierarchy) {
  # Split by '>'
  levels <- strsplit(hierarchy, " > ")[[1]]
  
  # Extract the first, second, third, and fourth levels if they exist
  first_level <- ifelse(length(levels) >= 1, levels[1], NA)
  second_level <- ifelse(length(levels) >= 2, levels[2], NA)
  third_level <- ifelse(length(levels) >= 3, levels[3], NA)
  fourth_level <- ifelse(length(levels) >= 4, levels[4], NA)
  
  # Return as a data frame row
  return(data.frame(
    First_Level = first_level,
    Second_Level = second_level,
    Third_Level = third_level,
    Fourth_Level = fourth_level,
    stringsAsFactors = FALSE
  ))
}

# Apply the function to all ATC hierarchies
atc_levels <- do.call(rbind, lapply(atc_hierarchy_rest$atc_hierarchy, split_hierarchy))
atc_hierarchy_rest = cbind(atc_hierarchy_rest, atc_levels)

###########Comparison organ and rest

I = intersect(atc_hierarchy_organ$Second_Level, atc_hierarchy_rest$Second_Level)
I1 = unique(atc_hierarchy_organ$Second_Level)

#keep only atc_hierarchy_rest that are in both atc_hierarchy_rest and atc_hierarchy_organ 
atc_hierarchy_rest_intesected = atc_hierarchy_rest[atc_hierarchy_rest$Second_Level %in% I,]
atc_hierarchy_organ_intesected = atc_hierarchy_organ[atc_hierarchy_organ$Second_Level %in% I,]


N_organ = data.frame(table(atc_hierarchy_organ_intesected$Second_Level))
N_rest = data.frame(table(atc_hierarchy_rest_intesected$Second_Level))

N_all_organ = length(unique(drug_adr_organ$drugbank_id))
N_all_rest = length(unique(drug_adr_rest$drugbank_id))

fisher_test_results <- lapply(N_organ$Var1, function(x) {
  count_organ <- N_organ[N_organ$Var1 == x, "Freq"]
  count_rest <- N_rest[N_rest$Var1 == x, "Freq"]
  
  # Prepare the contingency table
  contingency_table <- matrix(
    c(
      count_organ, N_all_organ - count_organ,
      count_rest, N_all_rest - count_rest
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  test_result <- fisher.test(contingency_table)
  
  # Return the p-value and odds ratio
  return(data.frame(
    Level = x,
    P_Value = test_result$p.value,
  ))
})

# Combine results into a single data frame
fisher_test_results_df <- do.call(rbind, fisher_test_results)
fisher_test_results_df$qval = p.adjust(fisher_test_results_df$P_Value, method = "BH")

write.table(fisher_test_results_df, "fisher_test_results_cardiac.csv", sep = ",")





