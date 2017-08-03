#' Small function that lists all the available clinical features of a given Expressionset
#'
#' @param Eset ExpressionSet
#'
#' @return list of all available clinical features
#' @export
#'
#' @examples
#' \dontrun{
#' list_clinical(Eset)
#' }
list_clinical <- function (Eset) {
  cols <- colnames(Biobase::pData(Eset))

 # factors <- list()
#  numbers <- list()

#  for(i in 1:length(colnames(Biobase::pData(Eset)))){
#    factors[[i]] <- levels(factor(Biobase::pData(Eset)[,i]))
#    numbers[[i]] <- unclass(table(factor(Biobase::pData(Eset)[,i])))
#  }

#  df <- data.frame(names(numbers[[1]][1:5]) )
#  for(i in 2:length(factors)){
#    df[,i] <- names(numbers[[i]][1:5])
#  }

#  View(df)
  cols
}
