# test for geneset:

## load the expression Data:
testset <- build_TCGA_Eset(clinical_file = "../../../Projects/TCGA-Elisa/Clinical_KIRC/KIRC.clin.merged.txt",
                           expression_file = "../../../Projects/TCGA-Elisa/RNA_KIRC/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt")

## add the mutation data
mutset <- fetch_TCGA_mutations(path = "../../../Projects/TCGA-Elisa/KIRC.Mutation_Packager_Oncotated_Calls/", Eset = testset)


# mutation: a vector containing the names of the mutation
# mutation_classification: a vector containing the mutation classifications
# variant_type: a vector containing the variant types

## make new columns for the mutation
pData(mutset)[,"IDH1 status"] <- ifelse(grepl("IDH1",pData(mutset)$mutation) & pData(mutset)$mutation != "Silent","mutated","wt")



Survival_adaptable(x = c("BCAT1"), Eset = mutset,
                   #additional = "patient.neoplasm_histologic_grade",
                   additional = "IDH1 status",
                   value = 0,
                   p.val=F,
                   xlabel="Days",
                   legend_position = "top",
                   average = "median",
                   optimal=F,
                   plot_cutpoint=F
                   )

Survival_adaptable(x = c("patient.neoplasm_histologic_grade"), Eset = testset,
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
                   additional = "mutations",
                   exclude=c("gx","g1"))

smoothCoxph_man(x = c("BAMBI", "CD44"), Eset=testset,
                #exclude_values=c("pathologic_stage", "Stage II", "Stage I", "Stage IV"),
                logrisk=T, average="median")
