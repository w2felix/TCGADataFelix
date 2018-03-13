#' TCGA Multivariate Analysis
#'
#' @param x A character vector with 1 or more gene names and/or columns of the clinical patient data
#' @param Eset An Expression Set
#' @param value Defines the value to subdivide gene expression groups.
#' Numeric: Devided into two groups (for more genes use value = c(cutpoint gene 1, cutpoint gene 2, ...))
#' "q": 25 \% quantile, 25-75 \% quantile and 75 \% quantile
#' @param factor_list Returns the levels of factors in the analysis
#' @param mutation Vector of mutations, that were added using the add_mutation() function
#' @param average for more than one gene, how the value of the averaged z-score is calculated, either median or mean
#' @param optimal calculate the optimal cutpoint, will overide the value value when numeric, does not work when value = "q"
#' @param survival For overall survival = "overall", for disease free survival: survival = "DFS"
#' @param return_df Returns the Dataframe that is used for calculatinf the Cox Regression
#' @param plot_cutpoint Plot the graph how the optimal cutpoint was calculated, normal survival plot will be REPLACED by surv_cutpoint: Determine the optimal cutpoint for each variable using ’maxstat’
#' @param coef Return the coefficients of the multivariate analysis instead of the multivariate object
#' @param stepwise Do a stepwise chosing a model by AIC to get rid of variables that don't play a role in the linear model
#' @param combined includes a comination of the different genes into the analysis. still in development
#' @param ... additional variables that can be added
#'
#' @return The summary of the multivariate analysis
#' @export
#'
#' @import survival
#'
#' @examples
#' \dontrun{
#' Multivariate(x = c("FOXA2", "GENDER"), Eset = Eset,
#' average = "median", stepwise=TRUE, coef=TRUE,
#' optimal=T, plot_cutpoint=F)
#' }
#'
#x <- "SPDEF"
Multivariate <- function (x,
                          Eset,
                          value = 0,
                          factor_list=F,
                          mutation,
                          average = "mean",
                          optimal = FALSE,
                          survival="overall",
                          return_df=FALSE,
                          plot_cutpoint=F,
                          coef=F,
                          combined=F,
                          stepwise=F,
                          ...) {

  if(survival=="overall") {
    Biobase::pData(Eset)$time <- as.numeric(Biobase::pData(Eset)$X_OS)
    Biobase::pData(Eset)$event <- Biobase::pData(Eset)$X_OS_IND
  } else if(survival=="DFS"){
    Biobase::pData(Eset)$time <- as.numeric(Biobase::pData(Eset)$X_DFS)
    Biobase::pData(Eset)$event <- Biobase::pData(Eset)$X_DFS_IND
  } else {
    stop(paste("You cannot use \"", survival, "\" as censoring, use \"overall\" or \"DFS\""), sep="")
  }

  # Check if variables exist in dataframe

  # which exist in pData?
  exist_in_pData <- x[x %in% colnames(Biobase::pData(Eset))]
  exist_in_exprs <- x[x %in% rownames(Biobase::exprs(Eset))]

  if(sum(x %in% c(exist_in_exprs,exist_in_pData))<length(x)){
    not_existing <- x[!x %in% colnames(Biobase::pData(Eset)) & !x %in% rownames(Biobase::exprs(Eset))]
    stop(paste(paste(not_existing,collapse = " & "), "does not exist in clinical and expression data"))
  }

  for(i in 1:length(exist_in_exprs)){
    Biobase::pData(Eset)[,exist_in_exprs] <- Biobase::exprs(Eset)[exist_in_exprs,]
    Biobase::pData(Eset)[,paste(exist_in_exprs,"_value",sep="")] <- Biobase::exprs(Eset)[exist_in_exprs,]
  }

  # remove NA

  Biobase::pData(Eset) <- Biobase::pData(Eset)[!is.na(Biobase::pData(Eset)$event), ]
  Biobase::pData(Eset) <- Biobase::pData(Eset)[!is.na(Biobase::pData(Eset)$time), ]
  if(length(x)>1){
    Biobase::pData(Eset) <-Biobase::pData(Eset)[apply(Biobase::pData(Eset)[,x],1,function (x) sum(!is.na(x)))==length(x), ]
    Biobase::pData(Eset) <-Biobase::pData(Eset)[apply(Biobase::pData(Eset)[,x],1,function (x) sum(x!="[Not Available]"))==length(x), ]
    Biobase::pData(Eset) <-Biobase::pData(Eset)[apply(Biobase::pData(Eset)[,x],1,function (x) sum(x!="[Not Applicable]"))==length(x), ]
  } else {
    Biobase::pData(Eset) <-Biobase::pData(Eset)[!is.na(x), ]
    Biobase::pData(Eset) <-Biobase::pData(Eset)[x!="[Not Available]", ]
    Biobase::pData(Eset) <-Biobase::pData(Eset)[x!="[Not Applicable]", ]
  }

  if(length(exist_in_exprs)>0){
    if(optimal){
      if(!missing(value)){
        warning("Value: When Optimal cutpoint is calculated, the value parameter will be ignored")
      }
      value_cutpoint <- survminer::surv_cutpoint(Biobase::pData(Eset), time="time", event="event", variables = exist_in_exprs)
      value_sum <- summary(value_cutpoint)
      value_single <- value_sum$cutpoint
      value_combined <- mean(value_single)

      if(plot_cutpoint){
        graphics::par(mfrow=c(2,1))
        return(graphics::plot(value_cutpoint, exist_in_exprs, palette = "npg"))
      }
      two_groups <- T
    } else if(missing(value)){
      stop("Please let either an optimal cutpoint be calculated by optimal=TRUE, use a user defined cutpoint by value=X or divide the groups into their quantiles by value=\"q\"")
     } else if(is.numeric(value)){
       if(plot_cutpoint){
         warning("plot_cutpoint: The Optimal cutpoint cannot be plotted when optimal is not TRUE, therefore plot_cutpoint is ignored")
       }
      if(length(value)==length(exist_in_exprs)){
        value_single <- value
      } else if(length(value)==1){
        for(i in 1:length(exist_in_exprs)){
          value_single[i] <- value
        }
      } else {
        stop("The cutpoint value needs to be either as long as genes to be analyzed: \"value = c(gene1, gene2, gene3,...)\"\nor value needs to have one value for all genes \"value = X\"")
      }

      value_combined <- mean(value)
      two_groups <- T
    } else if(value=="q") {
      two_groups <- F
    }

    if(two_groups){
      if(length(exist_in_exprs)>1){
        df_matrix <- matrix(t(Biobase::pData(Eset)[,exist_in_exprs]), ncol = length(exist_in_exprs), byrow=TRUE)
        Biobase::pData(Eset)$median <- Biobase::rowMedians(df_matrix)
        Biobase::pData(Eset)$mean <- base::rowMeans(df_matrix)
        if(average=="mean"){
          Biobase::pData(Eset)$combined <- ifelse(apply(Biobase::pData(Eset)[,exist_in_exprs],1,mean) > value_combined,
                                                  "High Expression",
                                                  "Low Expression")
        } else if(average=="median"){
          Biobase::pData(Eset)$combined <- ifelse(apply(Biobase::pData(Eset)[,exist_in_exprs],1,stats::median) > value_combined,
                                                  "High Expression",
                                                  "Low Expression")
        } else {
          stop("Multiple genes can only be combined using \"median\" or \"mean\" as the parameter for `average`to divide their expression")
        }

      }

      for(i in 1:length(exist_in_exprs)){
        Biobase::pData(Eset)[,exist_in_exprs[i]] <- ifelse(Biobase::pData(Eset)[,exist_in_exprs[i]] > value_single[i],
                                                           "High Expression",
                                                           "Low Expression")
      }
    } else {
      if(length(exist_in_exprs)>1){
        quantile_values <- apply(Biobase::pData(Eset)[,exist_in_exprs],2,function(x) stats::quantile(x,c(.25, .50, .75)))
        if(combined){
          df_matrix <- matrix(t(Biobase::pData(Eset)[,exist_in_exprs]), ncol = length(exist_in_exprs), byrow=TRUE)
          Biobase::pData(Eset)$median <- Biobase::rowMedians(df_matrix)
          Biobase::pData(Eset)$mean <- base::rowMeans(df_matrix)

          Biobase::pData(Eset)$combined <- ifelse(Biobase::pData(Eset)[,average] < mean(quantile_values[1,]), "Lower Quantile", ifelse(Biobase::pData(Eset)$Median < mean(quantile_values[3,]), "Intermediate", "Upper Quantile"))
        }

        for(i in 1:length(exist_in_exprs)){
          Biobase::pData(Eset)[,exist_in_exprs[i]] <- ifelse(Biobase::pData(Eset)[,exist_in_exprs[i]] < quantile_values[1,i], "Lower Quantile", ifelse(Biobase::pData(Eset)[,exist_in_exprs[i]] < quantile_values[3,i], "Intermediate", "Upper Quantile"))
        }

      } else {
        combined=FALSE
        quantile_values <- stats::quantile(Biobase::pData(Eset)[,exist_in_exprs], c(.25, .50, .75))
        Biobase::pData(Eset)[,exist_in_exprs] <- ifelse(Biobase::pData(Eset)[,exist_in_exprs] < quantile_values[1], "Lower Quantile", ifelse(Biobase::pData(Eset)[,exist_in_exprs] < quantile_values[3], "Intermediate", "Upper Quantile"))
      }
    }
  } else {
    if(optimal){
      warning("Optimal & Value: Optimal cutpoint cannot be calculated for clinical data and patients not divided by it or any given value, parameter optimal and value will be ignored")
    }
    if(!missing(value)){
      warning("Value: Clinical data cannot be divided into groups by a user defined gene expression value, parameter value and optimal will be ignored")
    }
    if(combined & length(exist_in_exprs)==1){
      warning("Combined: Combined gene expression cannot be calculated if only one gene there to be analyzed, parameter optimal will be ignored")
    } else if(combined){
      warning("Combined: Combined gene expression cannot be calculated if no genes are there to be analyzed, parameter optimal will be ignored")
    }
  }

  if(factor_list){
    if(length(x)>1){
      return_factor_list <- sapply(Biobase::pData(Eset)[,x],function (x) levels(as.factor(x)))
    } else {
      return_factor_list <- paste(x,": ",paste(levels(factor(Biobase::pData(Eset)[,x])),collapse = ", "), sep="")
    }
    if(return_df){
      warning("return_df: Please use \"factor_list = FALSE\" to return the data frame")
    }
    return(return_factor_list)
  }

  if(return_df){
    return(Biobase::pData(Eset))
  }

  if(length(exist_in_pData)!=0){
    pData_formula <- paste(paste("`", exist_in_pData, "`", sep=""),collapse = " + ")
  } else {
    pData_formula <- NULL
  }

  if(length(exist_in_exprs)!=0){
    exprs_formula <- paste(paste("`", exist_in_exprs, "`", sep=""),collapse = " + ")
    if(length(exist_in_pData!=0)){
      exprs_formula <- paste("+",exprs_formula)
    }
  } else {
    exprs_formula <- NULL
  }

  if(combined){
    additional_combined <- "+ combined"
  } else {
    additional_combined <- NULL
  }

  formula <- paste(pData_formula,exprs_formula,additional_combined)

  Survobject <- stats::as.formula(paste("survival::Surv(time, event) ~", formula))
  multivariate <-  do.call(survival::coxph,
                           list(formula = Survobject,
                                data = Biobase::pData(Eset)
                           )
  )

  if(!coef){
    return(multivariate)
  }

  if(stepwise){
    step_multivariate <- stats::step(multivariate)
    if(!coef){
      return(summary(step_multivariate))
    } else {
      return(coef(summary(step_multivariate)))
    }
  }

  coefficients_multivariate <- stats::coef(summary(multivariate))
  coefficients_multivariate

}
