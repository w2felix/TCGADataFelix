## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning = FALSE--------------------------------------
library(TCGADataFelix)

  # for the expressionset functions:
  library(Biobase)
  # for the modification of the output
  library(ggplot2)

## ------------------------------------------------------------------------
filepath <- "../../../../Projects/TCGA-Elisa/"

# The filename of the clinical data
clinical_data <- "Clinical_KIRC/KIRC.clin.merged.txt"  

# The filename of the expression data
expression_data <- "RNA_KIRC/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"

# The Folder of the mutation data (each patient has one file)
mutation_data <- "KIRC.Mutation_Packager_Oncotated_Calls/"

