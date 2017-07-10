#' Create an Expression Set from TCGA Data, downloaded from https://gdac.broadinstitute.org/
#'
#' @param clinical_file PATH of the clinical data from TCGA, e.g.: "../../../Projects/TCGA-Elisa/Clinical/LAML.clin.merged.txt"
#' @param expression_file PATH of the expression data from TCGA, e.g.: "../../../Projects/TCGA-Elisa/RNA/LAML.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
#' @param source Define the source of the data. "cBio" and "cBio_z_score" for cBioportal Data, "firehose" for Firehose data.
#' @importFrom dplyr %>%
#' @return ExpressionSet
#' @export
#'
#' @examples
#' \dontrun{
#' build_TCGA_Eset(clinical_file <- "../../../Projects/TCGA-Elisa/tcga_KIRC_cBio/data_bcr_clinical_data_patient.txt",
#' expression_file <- "../../../Projects/TCGA-Elisa/tcga_KIRC_cBio/Kdata_RNA_Seq_v2_mRNA_median_Zscores.txt",
#' source = "cBio_z_score",
#' normalization = FALSE)
#' }
#'
build_TCGA_Eset <- function (clinical_file,
                             expression_file,
                             source) {
#clinical_file <- "../../../Projects/TCGA-Elisa/tcga_KIRC_cBio/data_bcr_clinical_data_patient.txt"
#expression_file <- "../../../Projects/TCGA-Elisa/tcga_KIRC_cBio/data_RNA_Seq_v2_mRNA_median_Zscores.txt"

  #expression_file <- "../../../Projects/TCGA-Elisa/tcga_KIRC_cBio/data_RNA_Seq_v2_expression_median.txt"

if(source=="firehose") {

  ### Import clinical Data
  clinical <- as.data.frame(t(readr::read_delim(clinical_file,
                                                "\t", escape_double = FALSE, trim_ws = TRUE)))

  clinical %>% dplyr::mutate_if(is.factor, as.character) -> clinical1

  colnames(clinical1) <- clinical1[1,]
  clinical1 <- clinical1[-1,]
  rownames(clinical1) <- toupper(clinical1$patient.bcr_patient_barcode)

  # remove emtpy columns

  new_pData <- sapply(clinical1, function (k) all(is.na(k)))
  clinical2 <- clinical1[!new_pData]

  # View(clinical2)

  ### Import Expression Data
  # expression_file = "../../../Projects/TCGA-Elisa/RNA_KIRC/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
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

  vm <- function(x){
    x <- t(apply(x,1,as.numeric))
    ex <- limma::voom(x,plot=F)
    return(ex$E)
  }
  expression_removed_vm  <- vm(expression_tumor)

} else if(source=="cBio") {

  # Importing the clinical Data
  clinical <- as.data.frame(readr::read_delim(clinical_file,
                                                "\t", escape_double = FALSE, trim_ws = TRUE,skip = 4))

  #clinical %>% dplyr::mutate_if(is.factor, as.character) -> clinical1
  rownames(clinical) <- clinical$PATIENT_ID



  ### Import Expression Data
  expressiondata <- as.data.frame(readr::read_delim(expression_file,
                                                    "\t", escape_double = FALSE, trim_ws = TRUE))
 #View(expressiondata)
  rownames(expressiondata) <- expressiondata$Hugo_Symbol
  expressiondata <- expressiondata[,-c(1:2)]
  colnames(expressiondata) <- substr(colnames(expressiondata),1,12)

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
  vm <- function(x){
    x <- t(apply(x,1,as.numeric))
    ex <- limma::voom(x,plot=F)
    return(ex$E)
  }
  expression_removed_vm  <- vm(expression_removed)

} else if(source=="cBio_z_score") {

  # Importing the clinical Data
  clinical <- as.data.frame(readr::read_delim(clinical_file,
                                              "\t", escape_double = FALSE, trim_ws = TRUE,skip = 4))

  #clinical %>% dplyr::mutate_if(is.factor, as.character) -> clinical1
  rownames(clinical) <- clinical$PATIENT_ID



  ### Import Expression Data
  expressiondata <- as.data.frame(readr::read_delim(expression_file,
                                                    "\t", escape_double = FALSE, trim_ws = TRUE))
  #View(expressiondata)
  rownames(expressiondata) <- expressiondata$Hugo_Symbol
  expressiondata <- expressiondata[,-c(1:2)]
  colnames(expressiondata) <- substr(colnames(expressiondata),1,12)

  # first I remove genes whose expression is == NA in more than 50% of the samples, rest of NA is removed by the subsequent scripts:
  rem <- function(x){
    x <- as.matrix(x)
    x <- t(apply(x,1,as.numeric))
    r <- apply(x,1,function(x) length(which(is.na(x))))
    remove <- r > dim(x)[2]*0.5
    return(remove)
  }

  remove_data <- rem(expressiondata)
  expression_removed <- expressiondata[!remove_data,]

  ### end remove empty genes
  # as the data in the other options maybe voom normalized:

  expression_removed_vm <- as.matrix(expression_removed)
  # View(expression_removed_vm)
} else {
  stop("You have to define the source of the TCGA/GDC Data by source = \"cBio\" for \'data_RNA_Seq_v2_expression_median.txt\', source = \"cBio_z_score\" for \'data_RNA_Seq_v2_mRNA_median_Zscores.txt\' or for Firehose data: source = \"firehose\"")
}

### scaling the data:
  # is data normally distributed? approximately...
  # hist(as.matrix(expression_removed_vm),breaks = 10000,xlim = c(-5,5))

  # Z-Score the data -> scaling!
if(source=="cBio_z_score"){
  expression_removed_vm_zscore <- expression_removed_vm
} else {
  expression_removed_vm_zscore <- t(scale(t(expression_removed_vm), center = TRUE, scale = TRUE))
}

# hist(expression_removed_vm_zscore,xlim = c(-5,5))


## rename columns:

colnames(expression_removed_vm_zscore) <- gsub('\\.','-',substr(colnames(expression_removed),1,12))

# set the rownames keeping only gene name (otherwise the old rownames are still in the expression_removed_vm)
# take care not to have duplicates:

if(source=="firehose"){
  new_genenames <- sapply(rownames(expression_removed_vm_zscore), function(x) unlist(strsplit(x,'\\|'))[[1]])
  new_genenames_filtered <- ifelse(new_genenames=="?",rownames(expression_removed_vm_zscore),new_genenames)
  new_genenames_filtered_2 <- ifelse(duplicated(new_genenames_filtered),rownames(expression_removed_vm_zscore),new_genenames_filtered)
  rownames(expression_removed_vm_zscore) <- new_genenames_filtered_2
  }

expressionData <- expression_removed_vm_zscore


if(source=="firehose"){

  # create vector for death censoring
  # table(clinical2$patient.vital_status)

  # alive  dead
  # 67   133

  clinical2$X_OS_IND <- ifelse(clinical2$patient.vital_status == 'alive', 0,1)

  ind_keep <- grep('days_to_death',colnames(clinical2))
  death <- as.matrix(clinical2[,ind_keep])

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
  ind_keep <- grep('days_to_last_followup',colnames(clinical2))
  # View(fl)
  fl <- as.matrix(clinical2[,ind_keep])
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
  rownames(all_clin) <- rownames(clinical2)

  all_clin %>% dplyr::mutate_if(is.factor, as.character) -> all_clin_1
  all_clin_1 %>% dplyr::mutate_if(is.character, as.numeric) -> all_clin_2

  all_clin_2$new_death <- c()
  for (i in 1:length(all_clin_2$death_days)){
    all_clin_2$new_death[i] <- ifelse(is.na(all_clin_2$death_days[i]),
                                      all_clin_2$followUp_days[i],all_clin_2$death_days[i])
  }
  rownames(all_clin_2) <- rownames(clinical2)

  clinical2$X_OS <- all_clin_2$new_death

  # add age:
  clinical2$age <- floor(-as.numeric(clinical2$patient.days_to_birth)/365.2422)
  clinical <- clinical2
  # View(clinical2[,c("X_OS","X_OS_IND","patient.follow_ups.follow_up-2.days_to_last_followup","patient.follow_ups.follow_up.days_to_death", "patient.follow_ups.follow_up.vital_status" )])

} else {
  clinical$X_OS_IND <- ifelse(clinical$OS_STATUS == 'LIVING', 0,1)
  clinical$X_OS <- clinical$OS_MONTHS*30.4167

  clinical$X_DFS_IND <- ifelse(clinical$DFS_STATUS == 'DiseaseFree', 0, ifelse(clinical$DFS_STATUS == 'Recurred/Progressed', 1, NA))
  clinical$X_DFS <- as.numeric(clinical$DFS_MONTHS)*30.4167

  clinical$age <- floor(-as.numeric(clinical$DAYS_TO_BIRTH)/365.2422)
}

##### make expression set:

#### take identical patient samples ####
# get patient IDs that are in both clinical and expression data

# colnames(expression_removed_vm_zscore) %in% rownames(clinical)

names_clinical <- rownames(clinical)
names_genomic <- colnames(expressionData)

names_match_genomic <- match(names_clinical,names_genomic)
na_names.match_genomic <- stats::na.omit(names_match_genomic)
na_c_names.match_genomic <- as.numeric(na_names.match_genomic)

expData_ordered <- expressionData[, na_c_names.match_genomic]
names_genomic_ordered <- colnames(expData_ordered)

names_match_clinical <- match(names_genomic_ordered,names_clinical)
na_names.match_clinical <- stats::na.omit(names_match_clinical)
na_c_names.match_clinical <- as.numeric(na_names.match_clinical)

clinical_ordered <- clinical[na_c_names.match_clinical, ]

# Add a column of a gene with self defined cut off values

#### make Expression Set ####
# Combine information in an expression set object



Eset <- Biobase::ExpressionSet(expData_ordered, phenoData = Biobase::AnnotatedDataFrame(clinical_ordered))

Eset


}
