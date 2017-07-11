#' Title
#'
#' @param Eset Expression Set where the mutations will be added as a column in the pData table
#' @param mutations a vector of mutations to be added
#' @param multiple If more than one mutation was entered, a additional combined column will be added. multiple = "all" -> all mutations present in the patient, multiple = "one" -> mutation in at least one patient
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
                          mutations,
                          multiple="all"
                          ) {

  for(i in 1:length(mutations)){
    mutation_patients <- data.frame(unique(substr(mutationdata[mutationdata[,1]==mutations[i],16],1,12)),stringsAsFactors=FALSE)
    colnames(mutation_patients) <- "PATIENT_ID"
    mutation_patients[,mutations[i]] <- "mutated"
    Biobase::pData(Eset) <- dplyr::left_join(Biobase::pData(Eset),mutation_patients,by = c("PATIENT_ID" = "PATIENT_ID"))
    Biobase::pData(Eset)[,mutations[i]] <- ifelse(!is.na(Biobase::pData(Eset)[,mutations[i]]), "mutated", "not mutated")
  }

  if(length(mutations)>1){
    multiple_mut <- Biobase::pData(Eset)[,mutations]=="mutated"
    if(multiple=="all"){
      Biobase::pData(Eset)[,paste(mutations,  collapse = " & ")] <- ifelse(apply(multiple_mut, 1, sum, na.rm=TRUE)==length(mutations),"all mutated","not all mutated")
    } else if(multiple=="one") {
      Biobase::pData(Eset)[,paste(mutations,  collapse = " & ")] <- ifelse(apply(multiple_mut, 1, sum, na.rm=TRUE)>0,"mutated","not mutated")
    }

  }

  Eset

}
