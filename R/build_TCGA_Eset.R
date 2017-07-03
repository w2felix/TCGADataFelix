#' Create an Expression Set from TCGA Data, downloaded from https://gdac.broadinstitute.org/
#'
#' @param clinical_file PATH of the clinical data from TCGA, e.g.: "../../../Projects/TCGA-Elisa/Clinical/LAML.clin.merged.txt"
#' @param expression_file PATH of the expression data from TCGA, e.g.: "../../../Projects/TCGA-Elisa/RNA/LAML.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
#'
#' @importFrom dplyr %>%
#' @return ExpressionSet
#' @export
#'
#' @examples
#' \dontrun{
#' build_TCGA_Eset(clinical_file = "../../../Projects/TCGA-Elisa/Clinical/LAML.clin.merged.txt",
#' expression_file = "../../../Projects/TCGA-Elisa/RNA/LAML.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt")
#' }
#'
build_TCGA_Eset <- function (clinical_file,
                             expression_file) {
#clinical_file <- "../../../Projects/TCGA-Elisa/Clinical_KIRC/KIRC.clin.merged.txt"
#expression_file <- "../../../Projects/TCGA-Elisa/RNA_KIRC/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"

### Import clinical Data
clinical <- as.data.frame(t(readr::read_delim(clinical_file,
                                       "\t", escape_double = FALSE, trim_ws = TRUE)))

clinical %>% dplyr::mutate_if(is.factor, as.character) -> clinical1

colnames(clinical1) <- clinical1[1,]
clinical1 <- clinical1[-1,]
rownames(clinical1) <- toupper(clinical1$patient.bcr_patient_barcode)

# View(clinical1)

### Import Expression Data
expressiondata <- as.data.frame(readr::read_delim(expression_file,
                                           "\t", escape_double = FALSE, trim_ws = TRUE))
expressiondata <- expressiondata[-1,]
rownames(expressiondata) <- expressiondata$`Hybridization REF`
expressiondata <- expressiondata[,-1]

# first I remove genes whose expression is == 0 in more than 50% of the samples:
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}
remove_data <- rem(expressiondata)
expression_removed <- expressiondata[-remove_data,]

### end remove empty genes

## 14+15 position of tcga barcode sayswhat kind of tumor sample it is
# e.g.: 01 for tumor, or 11 for normal
# -> position 14: 0 = tumor, 1 = normal -> take normal

# table(substr(colnames(expression_removed),14,14))

# 0      1
# 534   72

# -> Now take only the tumor samples!

expression_tumor <- expression_removed[,substr(colnames(expression_removed),14,14)=="0"]

### normalize data:

# apply voom function from limma package to normalize the data

vm <- function(x){
  x <- t(apply(x,1,as.numeric))
  ex <- limma::voom(x,plot=F)
  return(ex$E)
}

expression_removed_vm  <- vm(expression_tumor)

# is data normally distributed? approximately...
# hist(expression_removed_vm)

# Z-Score the data -> scaling!
expression_removed_vm_zscore <- t(scale(t(expression_removed_vm), center = TRUE, scale = TRUE))

# is data normally distributed? approximately...
# hist(expression_removed_vm_zscore)

## rename columns:

colnames(expression_removed_vm_zscore) <- gsub('\\.','-',substr(colnames(expression_tumor),1,12))

# do clinical data fit the expression data?

# colnames(expression_removed_vm_zscore) %in% rownames(clinical1)

# set the rownames keeping only gene name (otherwise the old rownames are still in the expression_removed_vm)
# take care not to have duplicates:

new_genenames <- sapply(rownames(expression_removed_vm_zscore), function(x) unlist(strsplit(x,'\\|'))[[1]])
new_genenames_filtered <- ifelse(new_genenames=="?",rownames(expression_removed_vm_zscore),new_genenames)
new_genenames_filtered_2 <- ifelse(duplicated(new_genenames_filtered),rownames(expression_removed_vm_zscore),new_genenames_filtered)

rownames(expression_removed_vm_zscore) <- new_genenames_filtered_2
expressionData <- expression_removed_vm_zscore


#### match the patient ID in clinical data with the colnames of expression_removed_vm_zscore

# create vector for death censoring
# table(clinical1$patient.vital_status)

# alive  dead
# 67   133

clinical1$X_OS_IND <- ifelse(clinical1$patient.vital_status == 'alive', 0,1)


ind_keep <- grep('days_to_death',colnames(clinical1))
death <- as.matrix(clinical1[,ind_keep])

death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(clinical1))
View(fl)
fl <- as.matrix(clinical1[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum(is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

# and put everything together
all_clin <- data.frame(death_collapsed,fl_collapsed)
colnames(all_clin) <- c('death_days', 'followUp_days')
rownames(all_clin) <- rownames(clinical1)

all_clin %>% dplyr::mutate_if(is.factor, as.character) -> all_clin_1
all_clin_1 %>% dplyr::mutate_if(is.character, as.numeric) -> all_clin_2

all_clin_2$new_death <- c()
for (i in 1:length(all_clin_2$death_days)){
  all_clin_2$new_death[i] <- ifelse(is.na(all_clin_2$death_days[i]),
                                    all_clin_2$followUp_days[i],all_clin_2$death_days[i])
}
rownames(all_clin_2) <- rownames(clinical1)

clinical1$X_OS <- all_clin_2$new_death

# add age:
clinical1$age <- floor(-as.numeric(clinical1$patient.days_to_birth)/365.2422)


View(clinical1[,c("X_OS","X_OS_IND","patient.follow_ups.follow_up-2.days_to_last_followup","patient.follow_ups.follow_up.days_to_death", "patient.follow_ups.follow_up.vital_status" )])


##### make expression set:

#### take identical patient samples ####
# get patient IDs that are in both clinical and expression data

names_clinical <- rownames(clinical1)
names_genomic <- colnames(expressionData)

names_match_genomic <- match(names_clinical,names_genomic)
na_names.match_genomic <- stats::na.omit(names_match_genomic)
na_c_names.match_genomic <- as.numeric(na_names.match_genomic)

expData_ordered <- expressionData[, na_c_names.match_genomic]
names_genomic_ordered <- colnames(expData_ordered)

names_match_clinical <- match(names_genomic_ordered,names_clinical)
na_names.match_clinical <- stats::na.omit(names_match_clinical)
na_c_names.match_clinical <- as.numeric(na_names.match_clinical)

clinical_ordered <- clinical1[na_c_names.match_clinical, ]

# Add a column of a gene with self defined cut off values

#### make Expression Set ####
# Combine information in an expression set object



Eset <- Biobase::ExpressionSet(expData_ordered, phenoData = Biobase::AnnotatedDataFrame(clinical_ordered))

Eset


}
