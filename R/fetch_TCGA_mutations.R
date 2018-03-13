#' Add an Mutation Column to an Expression set made by build_TCGA_Eset() from TCGA Mutation Data, downloaded from https://gdac.broadinstitute.org/
#' Please download the Mutation_Packager_Oncotated_Calls
#'
#' @param path Path to the mutation Files downloaded from Firehose TCGA or path to the data_mutations_extended.txt file from cBioportal
#' @param source Use "cBio"/"cBio_z_score" or "firehose" as a source
#'
#' @return Returns the Expression Set with added mutation data
#' @export
#'
#' @examples
#' \dontrun{
#' mutset <- fetch_TCGA_mutations(path = "~/R/Projects/TCGA-Elisa/KIRC.Mutation_Packager_Oncotated_Calls/", Eset = Eset)
#' }

fetch_TCGA_mutations <- function(path,
                                 source) {
  if(source=="firehose"){
    filenames <- list.files(path = path, pattern = NULL, all.files = FALSE,
                            full.names = FALSE, recursive = FALSE,
                            ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

    filenames <- filenames[grepl("(TCGA)",filenames)]
    patientnames <- substr(filenames, 1, 12)

    patients <- list()
    mutation_file <- paste(path,filenames,sep="")
    patients = lapply(mutation_file, utils::read.delim, skip = 3)

    names(patients) <- patientnames

    # take only valid mutations:

    for(i in 1:length(filenames)){
      patients[[i]] <- patients[[i]][patients[[i]][25]=="Valid",]
    }

    mutationdata <- as.data.frame(patients[[1]])
    for(i in 2:length(patients)){
      mutationdata <- rbind(mutationdata,as.data.frame(patients[[i]]))
    }

    # read the Expressionset and add information according to patient
    # -> Data is collapsed to a single vector containing all mutations -> needs to be seperated again
  } else if(source=="cBio" | source=="cBio_z_score"){
    # path <- "../../../Projects/TCGA-Elisa/tcga_KIRC_cBio/data_mutations_extended.txt"

    mutationdata <- as.data.frame(readr::read_delim(path,
                                                  "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1))


  }

  mutationdata <<- mutationdata
}
