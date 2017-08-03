#' Create a Kaplan Meier Estimator by using a TCGA Dataset
#'
#' \code{survival_adaptable()} plots a highly modifyable survival estimator using the survminer package as a basis
#'
#' @param Eset An Expression Set
#' @param survival For overall survival = "overall", for disease free survival: survival = "DFS"
#' @param clinical a character vector with one or more factors of the clinical patient data
#' @param mutation a character vector with one or more factors of the patients mutation data added to the ExpressionSet by add_mutation()
#' @param expression a character vector with one or more gene expression data
#' @param optimal calculate the optimal cutpoint using the max_stat ranking method, will overide the value value
#' @param value Defines the value to subdivide the gene expression groups.
#' Numeric: Devided into two groups
#' "q": 25 \% quantile, 25-75 \% quantile and 75 \% quantile
#' @param gene_signature For more than one gene, how the genes for the signature are averaged, either "median" or "mean". If left out, all gene expression combinations will be plotted individually
#' @param groupings You can group factors from the clinical data
#' @param exclude Factors to exclude from any of the given clinical, expression or mutation data. The order of the to be exlcuded variables does not matter.
#' @param show_only Works only if exclude is used. If you want to exclude specific values but only show the regression for either "gene expression", "clinical" data or "mutation data"
#' @param p.val Displays the p-Value on the graph, only possible if a single factor is analyzed
#' @param logrisk Displays how the p-value was calculated
#' @param xlabel User defined x-axis label
#' @param ylabel User defined y-axis label
#' @param plot_title Title of the Plot, if not stated, no title will be shown
#' @param legend_position Where should the legend be? doesn't work with theme_bw
#' @param legend_rows How many rows are used for the legend
#' @param plot_cutpoint Plot the graph how the optimal cutpoint was calculated, normal survival plot will be REPLACED by surv_cutpoint: Determine the optimal cutpoint for each variable using ’maxstat’
#' @param risk_table If the risk table is shown or not, use FALSE or TRUE
#' @param return_df Returns the Dataframe that is used for calculatinf the Cox Regression
#' @param return_fit Returns the fit of the Regression. [[1]] = Fit, [[2]] = p-value, [[3]] = Coxph
#' @param factor_list Returns the list of factors that are in the dataframe. These factors can be used e.g. for exclusen.
#' @param ... Additional variables that can be passed over to the ggsurvplot function of the survminer package
#'
#' @return The modified survival Estimator from the survminer package to be able to handle TCGA Data
#' @export
#'
#' @import survival
#'
#'
#' @examples
#' \dontrun{
#' Survival_adaptable(expression = c("FOXA2"), Eset = Eset,
#' optimal = TRUE, p.val=TRUE,
#' xlabel="Days", legend_position = "top", average = "median",
#' plot_cutpoint=F,risk_table=TRUE)
#' }
Survival_adaptable <- function (Eset,
                                survival="overall",
                                clinical,
                                mutation,
                                expression,
                                optimal = FALSE,
                                value,
                                gene_signature,
                                groupings= FALSE,
                                exclude,
                                show_only=FALSE,
                                p.val = FALSE,
                                logrisk = TRUE,
                                xlabel,
                                ylabel,
                                plot_title = "",
                                legend_position,
                                legend_rows = 1,
                                plot_cutpoint=FALSE,
                                risk_table = TRUE,
                                return_df=FALSE,
                                return_fit=FALSE,
                                factor_list=FALSE,
                                ...) {

  if (missing(xlabel)){
    xlabel <- "time"
  }
  if (missing(ylabel)){
    ylabel <- "Survival"
  }

  if (missing(clinical) & missing(mutation) & missing(expression)) {
    stop("You have to define at least one gene, mutation or clinical feature you want to calculate survival estimators")
  }

  if(survival=="overall") {
    time <- as.numeric(Biobase::pData(Eset)$X_OS)
    event <- Biobase::pData(Eset)$X_OS_IND
  } else if(survival=="DFS"){
    time <- as.numeric(Biobase::pData(Eset)$X_DFS)
    event <- Biobase::pData(Eset)$X_DFS_IND
  } else {
    stop(paste("You cannot use \"", survival, "\" as censoring, use \"overall\" or \"DFS\""), sep="")
  }

  # Check if variables exist in dataframe
  not_existing_exprs <- NULL
  not_existing_clinical <- NULL
  not_existing_mutation <- NULL

  if(!missing(clinical)){
    exist_in_pData <- clinical[clinical %in% colnames(Biobase::pData(Eset))]
    if(length(exist_in_pData)<length(clinical)){
      not_existing_clinical <- clinical[!clinical %in% colnames(Biobase::pData(Eset))]
    }
  } else {
    exist_in_pData <- NULL
  }
  if(!missing(expression)){
    exist_in_exprs <- expression[expression %in% rownames(Biobase::exprs(Eset))]
    if(length(exist_in_exprs)<length(expression)){
      not_existing_exprs <- expression[!expression %in% rownames(Biobase::exprs(Eset))]
    }
  } else {
    exist_in_exprs <- NULL
  }
  if(!missing(mutation)){
    exist_in_mutation <- mutation[mutation %in% colnames(Biobase::pData(Eset))]
    if(length(exist_in_mutation)<length(mutation)){
      not_existing_mutation <- mutation[!mutation %in% colnames(Biobase::pData(Eset))]
    }
  } else {
    exist_in_mutation <- NULL
  }



  if(length(not_existing_mutation)!=0 | length(not_existing_exprs)!=0 | length(not_existing_clinical)!=0){
    if(length(not_existing_mutation)!=0){
      not_existing_mutation_stop <- paste(paste(not_existing_mutation,collapse = " & "), ": Not added to the ExpressionSet by add_mutation(); ",sep="")
    } else {
      not_existing_mutation_stop <- NULL
    }
    stop(paste(not_existing_mutation_stop,paste(c(not_existing_exprs,not_existing_clinical),collapse = " & "), ": not in the clinical and expression data", sep=""))
  }

  if(length(exist_in_exprs)==1){
    expression_data <- Biobase::exprs(Eset)[exist_in_exprs,]
  } else {
    expression_data <- t(Biobase::exprs(Eset)[exist_in_exprs,])
  }

  z <- data.frame(time, event, expression_data, Biobase::pData(Eset)[,exist_in_pData], Biobase::pData(Eset)[,exist_in_mutation])
  colnames(z) <- c("time","event",exist_in_exprs,exist_in_pData,exist_in_mutation)

  # remove NA

  for(i in 1:length(colnames(z))){
    z <- z[z[,i]!="[Not Available]",]
    z <- z[z[,i]!="[Not Applicable]",]
    z <- z[!is.na(z[,i]),]
  }

  if(dim(z)[1]==0){
    stop("No data to do further analysis, maybe your clinical data factor was empty?")
  }

  # Groupings

  if(groupings[1]!=FALSE){
    grouping_factors <- "groupings"
    r <- NULL
    l <- NULL
   # paste(c(exist_in_exprs,exist_in_mutation,exist_in_pData)[!groupings])

      for(j in 3:length(colnames(z))){
        grouping_indicator <- sapply(z[,j],function(x) x==groupings)
        i <- j-2
        r[[i]] <- apply(grouping_indicator,2,any)
        if(i>1){
          if(i>2){
            l <- r[[i]] + l
          } else {
            l <- r[[i]] + r[[i-1]]
          }

        }
      }
      z$groupings <- ifelse(l>0,paste("Group:",paste(groupings,collapse=", ")),"rest")
  } else {
    grouping_factors <- NULL
  }

  if (!missing(exclude)){
    for(i in 1:length(exclude)){
      for(j in 3:length(colnames(z))){
        z <- z[!z[,j]==exclude[i], ]
      }
    }
    if(length(rownames(z))==0){
      stop("All data was excluded with your exclude settings, no further analysis possible.")
    }
  }

  if(length(exist_in_exprs)>0){
    if(optimal){
      if(!missing(value)){
        warning("Value: When Optimal cutpoint is calculated, the value parameter will be ignored")
      }
      value_cutpoint <- survminer::surv_cutpoint(z, time="time", event="event", variables = exist_in_exprs)
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
        if(gene_signature=="mean" | gene_signature=="median"){
          df_matrix <- matrix(t(z[,exist_in_exprs]), ncol = length(exist_in_exprs), byrow=TRUE)
          z$median <- Biobase::rowMedians(df_matrix)
          z$mean <- rowMeans(df_matrix)
          if(gene_signature=="mean"){
            z$gene_signature <- ifelse(apply(z[,exist_in_exprs],1,mean) > value_combined,
                                       "High Expression",
                                       "Low Expression")
          } else if(gene_signature=="median"){
            z$gene_signature <- ifelse(apply(z[,exist_in_exprs],1,stats::median) > value_combined,
                                       "High Expression",
                                       "Low Expression")
          }
          gene_signature_factor <- "gene_signature"
        } else if(!missing(gene_signature)){
          stop("Multiple genes can only be combined using \"median\" or \"mean\" as the parameter for `gene_signature`to use them as a gene expression signature. If you do not want a gene expression signature, leave the parameter out")
        } else {
          gene_signature_factor <- NULL
        }
      } else {
        if(!missing(gene_signature) & length(exist_in_exprs)==1){
          warning("gene_signature: Gene expression signature cannot be calculated if only one gene there to be analyzed, parameter average will be ignored")
        }
        gene_signature_factor <- NULL
      }
      for(i in 1:length(exist_in_exprs)){
        z[,exist_in_exprs[i]] <- ifelse(z[,exist_in_exprs[i]] > value_single[i],
                                        "High Expression",
                                        "Low Expression")
      }
    } else {
      if(length(exist_in_exprs)>1){
        quantile_values <- apply(z[,exist_in_exprs],2,function(x) stats::quantile(x,c(.25, .50, .75)))
        if(gene_signature=="mean" | gene_signature=="median"){
          df_matrix <- matrix(t(z[,exist_in_exprs]), ncol = length(exist_in_exprs), byrow=TRUE)
          z$median <- Biobase::rowMedians(df_matrix)
          z$mean <- rowMeans(df_matrix)

          z$gene_signature <- ifelse(z[,gene_signature] < mean(quantile_values[1,]), "Lower Quantile", ifelse(z[,gene_signature] < mean(quantile_values[3,]), "Intermediate", "Upper Quantile"))
          gene_signature_factor <- "gene_signature"
        } else {
          gene_signature_factor <- NULL
        }

        for(i in 1:length(exist_in_exprs)){
          z[,exist_in_exprs[i]] <- ifelse(z[,exist_in_exprs[i]] < quantile_values[1,i], "Lower Quantile", ifelse(z[,exist_in_exprs[i]] < quantile_values[3,i], "Intermediate", "Upper Quantile"))
        }

      } else {
        gene_signature_factor <- NULL
        quantile_values <- stats::quantile(z[,exist_in_exprs], c(.25, .50, .75))
        z[,exist_in_exprs] <- ifelse(z[,exist_in_exprs] < quantile_values[1], "Lower Quantile", ifelse(z[,exist_in_exprs] < quantile_values[3], "Intermediate", "Upper Quantile"))
      }
    }

  } else {
    gene_signature_factor <- NULL
    if(optimal){
      warning("Optimal: Optimal cutpoint cannot be calculated for clinical data or mutation data, parameter optimal will be ignored")
    }
    if(!missing(value)){
      warning("Value: Clinical and mutational data cannot be divided into groups by a user defined gene expression value, parameter value and optimal will be ignored")
    }
    if(!missing(gene_signature)){
      warning("gene_signature: A gene expression signature cannot be calculated if no genes are there to be analyzed, parameter signature will be ignored")
    }
  }

  if(factor_list){
    if(length(colnames(z))>3){
      return_factor_list <- sapply(z[,c(exist_in_exprs,exist_in_pData,exist_in_mutation,gene_signature_factor,grouping_factors)],function (x) levels(as.factor(x)))
    } else {
      return_factor_list <- paste(colnames(z)[3],": ", paste(levels(factor(z[,3])),collapse = ", "), sep="")
    }
    if(return_df){
      warning("return_df: Please use \"factor_list = FALSE\" to return the data frame")
    }
    if(return_fit){
      warning("return_fit: Please use \"factor_list = FALSE\" to return the fit")
    }
    return(return_factor_list)
  }

  # Factorize
  if(length(c(exist_in_exprs,exist_in_pData,exist_in_mutation,gene_signature_factor,grouping_factors))>1){
    z[c(exist_in_exprs,exist_in_pData,exist_in_mutation,gene_signature_factor,grouping_factors)] <- lapply(z[,c(exist_in_exprs,exist_in_pData,exist_in_mutation,gene_signature_factor,grouping_factors)],as.factor)
  } else if(length(gene_signature_factor)!=0){
    z[,3] <- as.factor(z[,3])
  }

  if(return_df){
    return(z)
  }

  if(length(exist_in_pData)!=0 & show_only!="expression data" & show_only!="mutation data"){
    pData_formula <- paste(paste("`", exist_in_pData, "`", sep=""),collapse = " + ")
  } else {
    pData_formula <- NULL
  }

  if(length(exist_in_exprs)!=0 & show_only!="clinical" & show_only!="mutation data") {
    exprs_formula <- paste(paste("`", exist_in_exprs, "`", sep=""),collapse = " + ")
  #  if(length(exist_in_pData!=0)){
  #    exprs_formula <- paste("+",exprs_formula)
  #  }
  } else {
    exprs_formula <- NULL
  }

  if(length(exist_in_mutation)!=0 & show_only!="clinical" & show_only!="expression data" ){
    mutation_formula <- paste(paste("`", exist_in_mutation, "`", sep=""),collapse = " + ")
  #  if(length(exist_in_exprs!=0)){
  #    mutation_formula <- paste("+",mutation_formula)
  #  }
  } else {
    mutation_formula <- NULL
  }

  # Calculate survfit

  if(length(gene_signature_factor)>0 & length(exist_in_exprs)>1){
    if(groupings[1]!=FALSE){
      formula <- "groupings + gene_signature"
    } else {
      formula <- paste(c(pData_formula,mutation_formula,"gene_signature"),collapse = " + ")
    }

  } else {
    if(groupings[1]!=FALSE){
      formula <- paste(c("groupings",exprs_formula),collapse=" + ")
    } else {
      formula <- paste(c(pData_formula,mutation_formula,exprs_formula),collapse=" + ")
    }
  }

  # p-value:
  if(length(colnames(z))==3){
    z[,3] <- factor(z[,3])
  } else {
    z[,3:length(colnames(z))] <- lapply(z[,3:length(colnames(z))],factor)
  }

  count_exps <- nlevels(z[,exist_in_exprs])/length(exist_in_exprs)
  count_mut <- nlevels(z[,exist_in_mutation])*length(exist_in_mutation)
  count_clin <- nlevels(z[,exist_in_pData])*length(exist_in_pData)

  if(count_exps==0) count_exps=1
  if(count_mut==0) count_mut=1
  if(count_clin==0) count_clin=1

  if((count_exps*count_mut*count_clin)<3){
    if(p.val == TRUE){
      p.val = TRUE
    } else {
      p.val = FALSE
    }
  } else {
    p.val = FALSE
  }

   Survobject <- stats::as.formula(paste("Surv(time = time, event = event) ~", formula))

  fit <- do.call(survival::survfit,
                 list(formula = Survobject,
                      data = z,
                      type="kaplan-meier",
                      conf.type="log"
                      )
                 )


  if(return_fit){
    cox <- do.call(survival::coxph,
                   list(formula = Survobject,
                        data = z)
    )

    sdf <- do.call(survival::survdiff,
                   list(formula = Survobject,
                        data = z)
    )
    p.val <- 1 - stats::pchisq(sdf$chisq, length(sdf$n) - 1)

    quantile_survival <- quantile(fit, probs = c(0.05,0.25,0.5,0.75,0.95))

    return(list(fit,cox,paste("P-Value =",p.val),quantile_survival))
  }

  # Calculate scaling of X-Axis
  if(max(z$time,na.rm=TRUE)>500){
    breaks <- 500
  } else {
    breaks <- 50
  }

  #  if(length(colnames(z))>3){
  #    naming_factors <- sapply(z[,c(exist_in_exprs,exist_in_pData,exist_in_mutation,combined_factor)],function (x) levels(as.factor(x)))
  #  } else {
  #    naming_factors <- levels(factor(z[,3]))
  #  }

  if(length(gene_signature_factor)!=0){
    expression_title <- paste("gene signature of",paste(exist_in_exprs,collapse=" & "))
  } else {
    expression_title <- exist_in_exprs
  }
  legend_title <- paste(c(expression_title,exist_in_pData,exist_in_mutation),collapse=" & ")

  p <- survminer::ggsurvplot(fit,
                             data=z,
                             title = plot_title,
                             surv.scale = c("percent"),
                             #palette = "RdBu",
                             #xscale = "d_y",
                             legend = legend_position,
                             #legend.labs = naming_factors,
                             legend.title = legend_title,
                             # xlim = c(0,2000),        # present narrower X axis, but not affect survival estimates.
                             break.time.by = breaks,     # break X axis in time intervals by 500.

                             pval = p.val,             # show p-value of log-rank test.
                             pval.method = TRUE,
                             log.rank.weights = "n",
                             conf.int = F,         # show confidence intervals for point estimaes of survival curves.
                             xlab = xlabel,
                             ylab = ylabel,
                             risk.table = risk_table,       # show risk table.
                             #  risk.table.y.text.col = T, # colour risk table text annotations.
                             #  risk.table.y.text = F, # show bars instead of names in text annotations in legend of risk table

                             #ncensor.plot = TRUE
                             #cumevents=T
                             #cumcensor=T
                             # tables.theme = theme_cleantable(),
                             # ggtheme = theme_bw() # Change ggplot2 theme

                             surv.median.line="hv",
                             ...

  )
  p <- p + ggplot2::guides(colour = ggplot2::guide_legend(nrow = legend_rows))
  p
}
