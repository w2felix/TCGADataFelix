#' add a mutation column to the patient data
#'
#'
#' @param path Path to the mutation Files downloaded from Firehose TCGA
#' @param Eset Eset that will be used
#'
#' @return Returns the Expression Set with added mutation data
#' @export
#'
#' @examples fetch_TCGA_mutations(path = "../../../Projects/TCGA-Elisa/LAML.Mutation_Packager_Calls/", Eset = Eset)
#'
#'
#'
fetch_TCGA_mutations <- function(path,
                                 Eset) {

  filenames <- list.files(path = path, pattern = NULL, all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

  filenames <- filenames[grepl("(TCGA)",filenames)]
  patientnames <- substr(patientnames, 1, 12)

  patients <- list()
  mutation_file <- paste(path,filenames,sep="")
  patients = lapply(mutation_file, read.delim)

  names(patients) <- patientnames

  # take only valid mutations:

  for(i in 1:length(filenames)){
    patients[[i]] <- patients[[i]][patients[[i]][25]=="Valid",]
  }

  # read the Expressionset and add information according to patient
  # -> Data is collapsed to a single vector containing all mutations -> needs to be seperated again

  for(i in 1:length(filenames)){
    pData(Eset)[rownames(pData(Eset))==patientnames[i],"mutation"] <- paste(unlist(patients[[i]][1]),collapse=" ")
    pData(Eset)[rownames(pData(Eset))==patientnames[i],"mutation_classification"] <- paste(unlist(patients[[i]][9]),collapse=" ")
    pData(Eset)[rownames(pData(Eset))==patientnames[i],"variant_type"] <- paste(unlist(patients[[i]][10]),collapse=" ")
  }

  Eset
}
