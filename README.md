# Phenotypic Similarity of Adverse Drug Reactions and Disease Phenotypes: A Bridge to Mechanistic Discovery

## Project Description
Adverse drug reactions (ADRs) remain a significant challenge in drug development and patient safety, often due to our limited understanding of the underlying mechanisms. In this study, we investigated whether ADRs and disease phenotypes (DPs) that share similar clinical manifestations also share molecular mechanisms. 

## Main Findings
We constructed a comprehensive knowledge graph (KG) and applied graph representation learning techniques to quantify the mechanistic similarities between ADRs and DPs. Our analysis revealed substantial mechanistic overlaps within specific system organ classes (SOCs), including cardiac, psychiatric, metabolic, and renal disorders. This approach advocates for the prioritization of drugs that are less likely to result in ADRs, ultimately contributing to the development of safer and more targeted interventions.

## Install packages:

In *R*
```
igraph(2.0.3)
ggplot2(3.5.1)
fgsea(1.28.0)
clusterProfiler(4.10.1)
org.Hs.eg.db(3.18.0)
biomaRt(2.58.2)
ReactomePA(1.46.0)

```
In *Python*

To recreate the Conda environment:

```bash
conda env create -f environment.yml


In the following, the step-by-step instructions to run the pipeline and obtain the results are described:

## Knowledge Graph (KG) Construction

- **Homogeneous KG:**  
  Run the script `KG_Construction/GraphConstruction_homo.R` to construct a homogeneous KG for the **Node2Vec** and **GraphConv** models.

- **Heterogeneous KG:**  
  Run the script `KG_Construction/GraphConstruction_hetro.R` to construct a heterogeneous KG for the **RGCN** model.

- **Random KG:**  
  Run the script `KG_Construction/GraphConstruction_homo_random.R` to construct a random KG.

- **Reduced KG:**  
  Run the script `KG_Construction/GraphConstruction_reducedKG.R` to construct the reduced_KG.

- **KG with Removed Confounding Drugs:**  
  Run the script `KG_Construction/GraphConstruction_confounderRemoval.R` to construct a KG after removing confounding drugs.


*Hint:* Due to licensing restrictions, we cannot redistribute the drug–protein interaction data derived from DrugBank. Researchers interested in accessing this data can apply for an academic license directly through [DrugBank](https://www.drugbank.ca/).


## Main analysis

### Graph Representation Leraning

- **Node2vec:**  
  Run the script `main/` to generate embedding using **Node2vec** model.

- **GraphConv:**  
  Run the script `main/` to generate embedding using **GraphConv** model.

- **RGCN:**  
  Run the script `main/` to generate embedding using **RGCN** model.



### Embedding Evaluation

  Run the script `main/gsea_based_evaluation_adr_plot.R` and `gsea_based_evaluation_dp_plot.R` to reproduce Figure 1C.
  Run the script `main/embedding_based_biological_evaluation.R` and `main/embedding_based_biological_evaluation_plot.R` respectively to reproduce Figure 1D.
  
### Similarity Analysis

- **Phenotype similarity analysis.** 
  Run the script `main/phenotype_similarity_soc.R` to reproduce p-value of Phenotype similarity analysis for SOC similarities
  Run the script `main/phenotype_similarity_llt.R` to reproduce p-value of Phenotype similarity analysis for LLT similarities 

- **SOC similarity-based analysis.** 
  Run the script `main/soc_similarity_analysis.R` to reproduce p-values of SOC similarity-based analysis and Figure 2E

- **LLT similarity-based analysis.** 
  Run the script `main/llt_similarity_analysis.R` to reproduce p-values of LLT similarity-based analysis and Figure 2F

- **Pathway enrichment analysis**
  Run the script `main/fgsea_results.R`, `identified_genes_with_elbow_on_gsea_results.R`, and `pathway_enrichment_on_identified_genes.R` respectively to reproduce Figure 3

- **Number of Nodes**
  Run the script `main/number_of_nodes_associated_with_SOCs.R` to reproduce Figure 2G-I

### Robustness analysis

- **Robustness against confounding effects**
  Run the script `main/preprocessing_for_confounding_drug_removal_1.R`, `preprocessing_for_confounding_drug_removal_2` respectively to reproduce analysis related to Robustness against confounding drugs

- **Analysis of drug classes**
  Run the script `main/drug_ATC_classes.R` to reproduce analysis related to the enrichment of drug classes in each SOC









