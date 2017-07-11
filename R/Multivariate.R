#' TCGA Multivariate Analysis
#'
#' @param x A character vector with 1 or more gene names or a column of the clinical patient data
#' @param Eset An Expression Set
#' @param value Defines the value to subdivide the gene expression groups.
#' Numeric: Devided into two groups
#' "q": 25 \% quantile, 25-75 \% quantile and 75 \% quantile
#' @param factor Additional covariate for the multivariate analysis
#' @param factor_list Returns the levels of factors in the analysis
#' @param mutation Vector of mutations, that were added using the add_mutation() function
#' @param exclude Factors to exclude, when no gene is entered but another factor from the clinical patient data, it can also be excluded, the order of the to be exlcuded variables does not matter.
#' @param average for more than one gene, how the value of the averaged z-score is calculated, either median or mean
#' @param optimal calculate the optimal cutpoint, will overide the value value when numeric, does not work when value = "q"
#' @param survival For overall survival = "overall", for disease free survival: survival = "DFS"
#' @param return_df Returns the Dataframe that is used for calculatinf the Cox Regression
#' @param plot_cutpoint Plot the graph how the optimal cutpoint was calculated, normal survival plot will be REPLACED by surv_cutpoint: Determine the optimal cutpoint for each variable using ’maxstat’
#' @param coef Return the coefficients of the multivariate analysis instead of the multivariate object
#' @param ... additional variables that can be added
#'
#' @return The summary of the multivariate analysis
#' @export
#'
#' @examples
#' \dontrun{
#' Multivariate(x = c("FOXA2"), Eset = Eset,
#' value = 0, factor = "ETHNICITY",
#' exclude = c("LATINO"), average = "median",
#' optimal=T, plot_cutpoint=F)
#' }
#'

Multivariate <- function (x, Eset,
                          value = 0,
                          factor,
                          factor_list=F,
                          mutation,
                          exclude,
                          average = "mean",
                          optimal = FALSE,
                          survival="overall",
                          return_df=FALSE,
                          plot_cutpoint=F,
                          coef=F,
                          ...) {

  if(survival=="overall") {
    Biobase::pData(Eset)$time <- Biobase::pData(Eset)$X_OS
    Biobase::pData(Eset)$event <- Biobase::pData(Eset)$X_OS_IND
  } else if(survival=="DFS"){
    Biobase::pData(Eset)$time <- Biobase::pData(Eset)$DFS_MONTH
    Biobase::pData(Eset)$event <- ifelse(Biobase::pData(Eset)$DFS_STATUS=="Recurred/Progressed",1,0)
  } else {
    stop(paste("You cannot use \"", survival, "\" as censoring"), sep="")
  }

  if(x %in% colnames(Biobase::pData(Eset))){
    gene <- ""
  } else {
    gene <- x
  }
  #gene <- c("BCAT1","TP53")

  if(gene!=""){
    if(length(gene)>1){
      if(!gene %in% rownames(Biobase::exprs(Eset))){
        stop(paste(paste(gene[!gene %in% rownames(Biobase::exprs(Eset))], collapse=" & "), "not in gene or patient list"))
      }
      Biobase::pData(Eset)[,gene] <- Biobase::exprs(Eset)[gene,]

    } else {
      if(!gene %in% rownames(Biobase::exprs(Eset))){
        stop(paste(gene, "not in gene or patient list"))
      }
      Biobase::pData(Eset)[,gene] <- Biobase::exprs(Eset)[gene,]
    }
  } else {
    Biobase::pData(Eset)[,x] <- as.factor(Biobase::pData(Eset)[,x])
  }

  # Remove empty rows
  Biobase::pData(Eset) <- Biobase::pData(Eset)[!is.na(Biobase::pData(Eset)$event) & !is.na(Biobase::pData(Eset)$time),  ]
  #Eset <- mutset
  #factor <- "ETHNICITY"
  if (!missing(factor)) {
    for(i in 1:length(factor)){
      Biobase::pData(Eset)[,factor[i]] <- ifelse(Biobase::pData(Eset)[,factor[i]]=="[Not Available]",NA,Biobase::pData(Eset)[,factor[i]])
      Biobase::pData(Eset) <- Biobase::pData(Eset)[!is.na(Biobase::pData(Eset)[,factor[i]]), ]
     }
  }

  if(gene==""){
    Biobase::pData(Eset) <- Biobase::pData(Eset)[!is.na(Biobase::pData(Eset)[,x]),]
    Biobase::pData(Eset) <- droplevels(Biobase::pData(Eset)[,x])
  }


  if(gene!=""){
    if(length(gene)>1){

      ## Z-Score + normalization has already been done for the geneset itself

      #   for(i in 3:(length(z)-1)){
      #     Z_score_1 <- z[,i] - mean(z[,i])
      #     Z_score <- Z_score_1 / stats::sd(z[,i])
      #     z[,i] <- Z_score
      #   }

      df_matrix <- matrix(t(Biobase::pData(Eset)[,gene]), ncol = length(gene), byrow=TRUE)

      Biobase::pData(Eset)$Median <- Biobase::rowMedians(df_matrix)
      Biobase::pData(Eset)$Mean <- base::rowMeans(df_matrix)

      if(value!="q"){
        if(optimal){
          value_cutpoint <- survminer::surv_cutpoint(Biobase::pData(Eset), time="time", event="event", variables = gene)
          value_sum <- summary(value_cutpoint)
          value <- value_sum$cutpoint
          value <- mean(value)
          if(plot_cutpoint){
            graphics::par(mfrow=c(2,1))
            return(graphics::plot(value_cutpoint, gene, palette = "npg"))
          }
        }
        Biobase::pData(Eset)$expression <- ifelse(Biobase::pData(Eset)$Median > value, "High Expression", "Low Expression")
      }
      if(value == "q"){
        if(optimal){
          stop("optimal cutpoint not possible if separation of patients by quantiles")
        }
        quantile_values <- stats::quantile(Biobase::pData(Eset)$Median, c(.25, .50, .75))
        Biobase::pData(Eset)$expression <- ifelse(Biobase::pData(Eset)$Median < quantile_values[1], "Lower Quantile", ifelse(Biobase::pData(Eset)$Median < quantile_values[3], "Intermediate", "Upper Quantile"))
      }
      Biobase::pData(Eset)$expression <- as.factor(Biobase::pData(Eset)$expression)

    } else {
      ## Z-Score + normalization has already been done for the geneset itself
      # Z_score_1 <- z[,3] - mean(z[,3])
      # Z_score <- Z_score_1 / stats::sd(z[,3])
      # z[,3] <- Z_score

      if(value!="q"){
        if(optimal){
          value_cutpoint <- survminer::surv_cutpoint(Biobase::pData(Eset), time="time", event="event", variables = gene)
          value_sum <- summary(value_cutpoint)
          value <- value_sum$cutpoint
          if(plot_cutpoint){
            p <- graphics::plot(value_cutpoint, gene, palette = "npg")
            return(p)
          }
        }
        Biobase::pData(Eset)$expression <- ifelse(Biobase::pData(Eset)[,gene] > value, "High Expression", "Low Expression")
      } else if(value == "q"){
        if(optimal){
          stop("optimal cutpoint not possible if separation of patients by quantiles")
        }
        quantile_values <- stats::quantile(Biobase::pData(Eset)[,gene], c(.25, .50, .75))
        Biobase::pData(Eset)$expression <- ifelse(Biobase::pData(Eset)[, gene] < quantile_values[1], "Lower Quantile", ifelse(Biobase::pData(Eset)[, gene] < quantile_values[3], "Intermediate", "Upper Quantile"))
      }
    }

    Biobase::pData(Eset)$expression <- as.factor(Biobase::pData(Eset)$expression)
  }

  if(factor_list){
    factor_list <- sapply(Biobase::pData(Eset)[,factor],function (x) levels(as.factor(x)))
    return(factor_list)
  }

  if(return_df){
    return(Biobase::pData(Eset))
  }

  multivariate <- survival::coxph(stats::as.formula(paste("Surv(time, event) ~", paste(factor, collapse = " + "), "+ expression")), data = Biobase::pData(Eset))

  if(!coef){
    return(multivariate)
  }

  coefficients_multivariate <- stats::coef(summary(multivariate))

  coefficients_multivariate

}
