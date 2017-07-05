library(TCGADataFelix)

### EDIT HERE:
# filepath for TCGA patient data downloaded from GDAC firebrowse (http://gdac.broadinstitute.org/):
# Have one folder for clinical data, one for expression data and one for the mutational data

filepath <- "../../../Projects/TCGA-Elisa/"
# if filepath = working directory use:
# filepath <- ""

# Filenames and Folder where the clinical and expression data is stored
clinical_data <- "Clinical_KIRC/KIRC.clin.merged.txt"
expression_data <- "RNA_KIRC/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"

# for the mutations, the Oncotated Calls are downloaded from GDAC firebrowse (http://gdac.broadinstitute.org/)
# now you need to state the FOLDER where the mutations are in:

mutation_data <- "KIRC.Mutation_Packager_Oncotated_Calls/"

### END EDIT HERE

      ## load the expression Data:
      Eset <- build_TCGA_Eset(clinical_file = paste(filepath,clinical_data,sep=""),
                                 expression_file = paste(filepath,expression_data,sep=""))
      ## add the mutation data
      mutset <- fetch_TCGA_mutations(path = paste(filepath,mutation_data,sep=""), Eset = Eset)

        # added columns to the clinical data, stored in pData(Eset)
        # mutation: a vector containing the names of the mutation
        # mutation_classification: a vector containing the mutation classifications
        # variant_type: a vector containing the variant types

### Edit here:

gene <- "FOXD1"
mutation <- "PBRM1"

### End of edit

      ## make new columns for the mutation
      pData(mutset)[,paste(mutation,"status")] <- ifelse(grepl(mutation,pData(mutset)$mutation) & pData(mutset)$mutation != "Silent","mutated","wt")

      ## example to make other categories, e.g. age:
      pData(mutset)[,"Age category"] <- ifelse(pData(mutset)$age > 65,"old","young")

### Adapt here to your needs:

Survival_adaptable(x = c(gene), Eset = mutset,
                   #additional = "patient.neoplasm_histologic_grade",
                   additional = paste(mutation,"status"),
                   value = 0,
                   p.val=F,
                   xlabel="Days",
                   legend_position = "top",
                   average = "median",
                   optimal=F,
                   plot_cutpoint=F
                   )

Survival_adaptable(x = c("patient.neoplasm_histologic_grade"), Eset = mutset,
                   value = 0,
                   p.val=TRUE,
                   xlabel="Days",
                   legend_position = "right",
                   average = "median",
                   optimal=T,
                   plot_cutpoint=F,
                   additional = "patient.gender",
                   exclude=c("gx","g1","male"))

Survival_adaptable(x = c("patient.gender"), Eset = mutset,
                   value = 0,
                   p.val=TRUE,
                   xlabel="Days",
                   legend_position = "right",
                   average = "median",
                   optimal=T,
                   plot_cutpoint=F,
                   additional = paste(mutation,"status"),
                   exclude=c("gx","g1")
                   )

Survival_adaptable(x = c("Age category"), Eset = mutset,
                   value = 0,
                   p.val=TRUE,
                   xlabel="Days",
                   legend_position = "right",
                   average = "median",
                   optimal=T,
                   plot_cutpoint=F,
                   additional = paste(mutation,"status")
                   #exclude=c("gx","g1")
                   )

smoothCoxph_man(x = c("BAMBI", "CD44"), Eset=mutset,
                #exclude_values=c("pathologic_stage", "Stage II", "Stage I", "Stage IV"),
                logrisk=T, average="median")
