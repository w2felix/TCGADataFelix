#' Title
#'
#' @param Eset Expression Set where the mutations will be added as a column in the pData table
#' @param mutations a vector of mutations to be added
#'
#' @return Expression Set where the mutations will are added as a column in the pData table
#' @export
#'
#' @examples
#' \dontrun{
#' add_mutation(Eset = Eset,
#' mutations = c("VHL", "PBRM1")
#' )
#' }
add_mutation <- function (Eset,
                          mutations
                          ) {

  for(i in 1:length(mutations)){
    mutation_patients <- data.frame(unique(substr(mutationdata[mutationdata[,1]==mutations[i],16],1,12)),stringsAsFactors=FALSE)
    colnames(mutation_patients) <- "PATIENT_ID"
    mutation_patients[,mutations[i]] <- "mutated"
    Biobase::pData(Eset) <- dplyr::left_join(pData(Eset),mutation_patients,by = c("PATIENT_ID" = "PATIENT_ID"))
  }

  if(length(mutations)>1){
    multiple <- !is.na(Biobase::pData(Eset)[,mutations])
    Biobase::pData(Eset)[,paste(mutations,  collapse = " & ")] <- ifelse(apply(multiple, 1, sum, na.rm=TRUE)==length(mutations),"all mutated","not all mutated")
  }

  Eset

}
