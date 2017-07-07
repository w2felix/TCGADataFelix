#' Create a Kaplan Meier Estimator by using a TCGA Dataset
#'
#' \code{survival_adaptable()} plots a highly modifyable survival estimator using the survminer package as a basis
#'
#' @param x a character vector with 1 or more gene names or a column of the clinical patient data
#' @param Eset An Expression Set
#' @param value Defines the value to subdivide the gene expression groups.
#' Numeric: Devided into two groups
#' "q": 25 \% quantile, 25-75 \% quantile and 75 \% quantile
#' @param additional Additional covariate to subset groups
#' @param exclude Factors to exclude, when no gene is entered but another factor from the clinical patient data, it can also be excluded, the order of the to be exlcuded variables does not matter.
#' @param p.val display the p-Value on the graph
#' @param xlabel User defined x-axis label
#' @param ylabel User defined y-axis label
#' @param logrisk logrisk?!
#' @param legend_position where should the legend be? doesn't work with theme_bw
#' @param legend_rows How many rows are used for the legend
#' @param average for more than one gene, how the value of the averaged z-score is calculated, either median or mean
#' @param optimal calculate the optimal cutpoint, will overide the value value when numeric, does not work when value = "q"
#' @param plot_cutpoint Plot the graph how the optimal cutpoint was calculated, normal survival plot will be REPLACED by surv_cutpoint: Determine the optimal cutpoint for each variable using ’maxstat’
#' @param risk_table If the risk table is shown or not, use FALSE or TRUE
#' @param plot_title Title of the Plot, if not stated, no title will be shown
#' @param ... additional variables that can be added
#'
#' @return A survival Estimator
#' @export
#'
#' @examples
#' \dontrun{
#' Survival_adaptable(x = c("FOXA2"), Eset = Eset,
#' value = 0, additional = "pathologic_stage",
#' exclude = c("Stage IV"), p.val=TRUE,
#' xlabel="Days", legend_position = "top", average = "median",
#' optimal=T, plot_cutpoint=F,risk_table=TRUE)
#' }
Survival_adaptable <- function (x, Eset,
                                value = 0,
                                additional,
                                exclude,
                                p.val = FALSE,
                                xlabel,
                                ylabel,
                                logrisk = TRUE,
                                legend_position,
                                legend_rows = 1,
                                average = "mean",
                                optimal = FALSE,
                                plot_cutpoint=FALSE,
                                risk_table = TRUE,
                                plot_title = "",
                                ...) {


  if (missing(x)) {
    stop("You have to define a gene you want to calculate survival")
  }

  time <- Biobase::pData(Eset)$X_OS
  event <- Biobase::pData(Eset)$X_OS_IND

  if(x %in% colnames(Biobase::pData(Eset))){
    gene <- ""
  } else {
    gene <- x
  }

  if (missing(xlabel)){
    xlabel <- "time"
  }
  if (missing(ylabel)){
    ylabel <- "Survival"
  }

  # get the data to the Surv variables
  # Adding a column with gene expression labels to the phenotype, updates if new values

  z <- data.frame(time = time, event = event)
  if(gene!=""){
    if(length(gene)>1){

      z <- data.frame(time = time, event = event)
      z[,gene] <- t(Biobase::exprs(Eset)[gene,])
      rownames(z) <- colnames(Biobase::exprs(Eset)[gene,])

    } else {
      if(!gene %in% rownames(Biobase::exprs(Eset))){
        stop(paste(gene, "not in gene or patient list"))
      }
      z[,gene] <- Biobase::exprs(Eset)[gene,]
      rownames(z) <- colnames(Biobase::exprs(Eset)[gene,])
    }
  } else {
    z[,x] <- as.factor(Biobase::pData(Eset)[,x])
  }


  if (missing(additional)){
    z$additional <- 0
  } else if (is.na(additional)){
    z$additional <- 0
  } else if (additional==""){
    z$additional <- 0
  } else if (!additional %in% colnames(Biobase::pData(Eset))){
    stop(paste(additional, "not in phenotype list"))
  } else {
    z$additional <- as.factor(Biobase::pData(Eset)[,additional])
  }

  # Get rid of empty rows
  z <- z[!is.na(z$event) & !is.na(z$time) & !is.na(z$additional) & !z$additional=="",  ]
  if(gene==""){
    z <- z[!is.na(z[,x]),]
  }


  if (!missing(exclude)){
    if(length(exclude)>1){
      for(i in 1:length(exclude)){
        z <- z[!z$additional==exclude[i], ]
        if(gene==""){
          z <- z[!z[,x]==exclude[i], ]
        }
      }
    } else {
      z <- z[!z$additional==exclude, ]
      if(gene==""){
        z <- z[!z[,x]==exclude, ]
      }
    }
  }

  if (z$additional[1]!=0 && is.factor(z$additional)){
    z$additional <- droplevels(z$additional)
  }

  if(gene==""){
    z[,x] <- droplevels(z[,x])
  }

  if(gene!=""){
    if(length(gene)>1){

    ## Z-Score + normalization has already been done for the geneset itself

 #   for(i in 3:(length(z)-1)){
 #     Z_score_1 <- z[,i] - mean(z[,i])
 #     Z_score <- Z_score_1 / stats::sd(z[,i])
 #     z[,i] <- Z_score
 #   }

      df_matrix <- matrix(t(z[,gene]), ncol = length(gene), byrow=TRUE)

      z$Median <- Biobase::rowMedians(df_matrix)
      z$Mean <- base::rowMeans(df_matrix)

      if(value!="q"){
        if(optimal){
          value_cutpoint <- survminer::surv_cutpoint(z, time="time", event="event", variables = gene)
          value_sum <- summary(value_cutpoint)
          value <- value_sum$cutpoint
          value <- mean(value)
          if(plot_cutpoint){
            graphics::par(mfrow=c(2,1))
            return(graphics::plot(value_cutpoint, gene, palette = "npg"))
          }
        }
        z$gene <- ifelse(z$Median > value, "High Expression", "Low Expression")
      }
      if(value == "q"){
        if(optimal){
          stop("optimal cutpoint not possible if separation of patients by quantiles")
        }
        quantile_values <- stats::quantile(z$Median, c(.25, .50, .75))
        z$gene <- ifelse(z$Median < quantile_values[1], "Lower Quantile", ifelse(z$Median < quantile_values[3], "Intermediate", "Upper Quantile"))
      }
      z$gene <- as.factor(z$gene)

    } else {
    ## Z-Score + normalization has already been done for the geneset itself
   # Z_score_1 <- z[,3] - mean(z[,3])
   # Z_score <- Z_score_1 / stats::sd(z[,3])
   # z[,3] <- Z_score

      if(value!="q"){
        if(optimal){
          value_cutpoint <- survminer::surv_cutpoint(z, time="time", event="event", variables = gene)
          value_sum <- summary(value_cutpoint)
          value <- value_sum$cutpoint
          if(plot_cutpoint){
            p <- graphics::plot(value_cutpoint, gene, palette = "npg")
            return(p)
          }
        }
        z$gene <- ifelse(z[,3] > value, "High Expression", "Low Expression")
      } else if(value == "q"){
        if(optimal){
          stop("optimal cutpoint not possible if separation of patients by quantiles")
        }
        quantile_values <- stats::quantile(z[,3], c(.25, .50, .75))
        z$gene <- ifelse(z[, 3] < quantile_values[1], "Lower Quantile", ifelse(z[, 3] < quantile_values[3], "Intermediate", "Upper Quantile"))
      }
    }

    z$gene <- as.factor(z$gene)
  }


  if(gene!=""){
    if(length(gene)>1){
      for(i in 1:length(gene)){
        z <- z[!is.na(z[,gene[i]]), ]
      }
      if(average=="mean"){
        z <- z[order(z$Mean), ]
        coxph1 <- survival::coxph(survival::Surv(time, event) ~ Mean, data = z)
      } else {
        z <- z[order(z$Median), ]
        coxph1 <- survival::coxph(survival::Surv(time, event) ~ Median, data = z)
      }
    } else {
      z <- z[!is.na(z$gene), ]
      z <- z[order(z$gene), ]
      coxph1 <- survival::coxph(survival::Surv(time, event) ~ gene, data = z)
    }
  } else {
    z <- z[order(z[,x]), ]
    z2 <- z
    colnames(z2)[3] <- "factor"
    coxph1 <- survival::coxph(survival::Surv(time, event) ~ factor, data = z2)
  }

  coeffs <- stats::coef(summary(coxph1))



  if (z$additional[1]!=0){
    if(p.val==TRUE){

      # The second survdiff() argument, rho, designates the weights according to ^S(t)p
      # and may be any numeric  value. The default is rho=0
      # and corresponds to the log-rank test. The Peto& Peto modi cation of the Gehan-Wilcoxon test
      # is computed using rho=1
      if(gene!=""){
        p_value_survdiff <- survival::survdiff(survival::Surv(time, event) ~ gene + additional, data = z, rho=1)
        p.value = 1 - stats::pchisq(p_value_survdiff$chisq, length(p_value_survdiff$n) - 1)
      } else {
        p_value_survdiff <- survival::survdiff(survival::Surv(time, event) ~ factor + additional, data = z2, rho=1)
        p.value = 1 - stats::pchisq(p_value_survdiff$chisq, length(p_value_survdiff$n) - 1)
      }

    }

    # get Labels for legend

    z$additional <- droplevels(z$additional)
    if(gene!=""){
      z$gene <- droplevels(z$gene)
      numbers_quantiles <- data.frame(table(z$additional, z$gene))

      # t 1/2

      median_survival <- survival::survfit(survival::Surv(time, event) ~ gene + additional, data = z)
      median_survival_df <- summary(median_survival)
      #     median_survival_df$table[i,"median"]

      Nameing_vector_all <- c(paste(numbers_quantiles[,2], " - ", additional, ": ", numbers_quantiles[,1], "; n = ", numbers_quantiles[,3], "; t1/2 = ",     median_survival_df$table[,"median"], sep=""))
      Nameing_vector <- c(paste(numbers_quantiles[,2], " Expression - ", additional, ": ", numbers_quantiles[,1], sep=""))

      Nameing_vector_all <- Nameing_vector[numbers_quantiles$Freq!=0]
      Nameing_vector <- Nameing_vector[numbers_quantiles$Freq!=0]

      if(length(gene)>1){
        label_mean <- ifelse(average=="mean","Mean","Median")
        if(length(gene)>5){
          length_gene_rounded <- floor(length(gene)/5)
          gene_names <- c()
          for(i in 1:length_gene_rounded){
            start <- 1+((i-1)*5)
            end <- i*5
            gene_names <- paste(gene_names, paste(gene[start:end],collapse=" & "), "\n")
          }
          if(length(gene)%%5!=0){
            end_length <- length_gene_rounded*5
            gene_names <- paste(gene_names, paste(gene[end_length:length(gene)],collapse=" & "), "\n")
          }
        } else {
          gene_names <- paste(gene,collapse=" & ")
        }
        gene_name <- paste(label_mean, "of", gene_names)
      } else {
        gene_name <- gene[1]
      }

     } else {
        median_survival <- survival::survfit(survival::Surv(time, event) ~ factor + additional, data = z2)
        median_survival_df <- summary(median_survival)
        #     median_survival_df$table[i,"median"]


        numbers_quantiles <- data.frame(table(z$additional, z[,x]))
        Nameing_vector_all <- c(paste(numbers_quantiles[,2], " - ", additional, ": ", numbers_quantiles[,1], "; n = ", numbers_quantiles[,3], "; t1/2 = ", median_survival_df$table[,"median"], sep=""))
        Nameing_vector <- c(paste(numbers_quantiles[,2], "  - ", additional, ": ", numbers_quantiles[,1], sep=""))

        Nameing_vector_all <- Nameing_vector[numbers_quantiles$Freq!=0]
        Nameing_vector <- Nameing_vector[numbers_quantiles$Freq!=0]

        gene_name <- x
      }

    #survival function:
    if(gene!=""){
      fit <- survival::survfit(survival::Surv(time = time,
                                              event = event) ~ gene + additional,
                               data = z,
                               type="kaplan-meier",
                               conf.type="log")
    } else {
      fit <- survival::survfit(survival::Surv(time = time,
                                              event = event) ~ factor + additional,
                               data = z2,
                               type="kaplan-meier",
                               conf.type="log")
      z <- z2
    }

    p <- survminer::ggsurvplot(
      fit,                     # survfit object with calculated statistics.
      data=z,
      title = plot_title,
      surv.scale = c("percent"),
      legend = legend_position,
      legend.labs = Nameing_vector,
      legend.title = gene_name,
      # xlim = c(0,2000),        # present narrower X axis, but not affect survival estimates.
      break.time.by = 500,     # break X axis in time intervals by 500.

      pval = p.val,             # show p-value of log-rank test.
      pval.method = TRUE,
      #log.rank.weights = "n",
      conf.int = F,         # show confidence intervals for point estimaes of survival curves.
      xlab = xlabel,
      ylab = ylabel,
      risk.table = risk_table,       # show risk table.
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = F, # show bars instead of names in text annotations in legend of risk table

      #ncensor.plot = TRUE
      # cumevents=T
      # cumcensor=T

      surv.median.line="hv"
      # tables.theme = theme_cleantable(),
      # ggtheme = theme_bw() # Change ggplot2 theme
    )
    p <- p + ggplot2::guides(colour = ggplot2::guide_legend(nrow = legend_rows))
    p

  } else {
    if(p.val==TRUE){

      # The second survdiff() argument, rho, designates the weights according to ^S(t)p
      # and may be any numeric  value. The default is rho=0
      # and corresponds to the log-rank test. The Peto& Peto modi cation of the Gehan-Wilcoxon test
      # is computed using rho=1
      if(gene!=""){
        p_value_survdiff <- survival::survdiff(survival::Surv(time, event) ~ gene, data = z, rho=1)
        p.value = 1 - stats::pchisq(p_value_survdiff$chisq, length(p_value_survdiff$n) - 1)
      } else {
        p_value_survdiff <- survival::survdiff(survival::Surv(time, event) ~ factor, data = z2, rho=1)
        p.value = 1 - stats::pchisq(p_value_survdiff$chisq, length(p_value_survdiff$n) - 1)
      }

    }

    # t 1/2
    if(gene!=""){
      median_survival <- survival::survfit(survival::Surv(time, event) ~ gene, data = z)
      median_survival_df <- summary(median_survival)
      #     median_survival_df$table[i,"median"]

      z$gene <- as.factor(z$gene)
      gene_levels <- nlevels(z$gene)
      Nameing_vector <- c()
      for(i in 1:gene_levels) {
        Nameing_vector <- c(Nameing_vector, paste(levels(z$gene)[i], "; n = ", nrow(z[z$gene==levels(z$gene)[i], ]), sep=""))
      }

      if(length(gene)>1){
        label_mean <- ifelse(average == "mean","Mean","Median")
        if(length(gene)>5){
          length_gene_rounded <- floor(length(gene)/5)
          gene_names <- c()
          for(i in 1:length_gene_rounded){
            start <- 1+((i-1)*5)
            end <- i*5
            gene_names <- paste(gene_names, paste(gene[start:end],collapse=" & "), "\n")
          }
          if(length(gene)%%5!=0){
            end_length <- length_gene_rounded*5
            gene_names <- paste(gene_names, paste(gene[end_length:length(gene)],collapse=" & "), "\n")
          }
        } else {
          gene_names <- paste(gene,collapse=" & ")
        }
        gene_name <- paste(label_mean, "of", gene_names)
      } else {
        gene_name <- gene[1]
      }
    } else {
      median_survival <- survival::survfit(survival::Surv(time, event) ~ factor, data = z2)
      median_survival_df <- summary(median_survival)
      #     median_survival_df$table[i,"median"]

      gene_levels <- nlevels(z[,x])
      Nameing_vector <- c()
      for(i in 1:gene_levels) {
        Nameing_vector <- c(Nameing_vector, paste(levels(z[,x])[i], "; n = ", nrow(z[z[,x]==levels(z[,x])[i], ]), sep=""))
      }
      gene_name <- x
    }

    # Get p value for the difference of the gene and plot them on the graph

    if(gene!=""){
      fit <- survival::survfit(survival::Surv(time = time,
                                              event = event) ~ gene,
                               data = z,
                               type="kaplan-meier",
                               conf.type="log")
    } else {
      fit <- survival::survfit(survival::Surv(time = time,
                                              event = event) ~ factor,
                               data = z2,
                               type="kaplan-meier",
                               conf.type="log")
      z <- z2
    }

    p <- survminer::ggsurvplot(
      fit,                     # survfit object with calculated statistics.
      data=z,
      surv.scale = c("percent"),
      legend = legend_position,
      title = plot_title,
      legend.labs = Nameing_vector,
      legend.title = gene_name,
      # xlim = c(0,2000),        # present narrower X axis, but not affect survival estimates.
      break.time.by = 500,     # break X axis in time intervals by 500.

      pval = p.val,             # show p-value of log-rank test.
      pval.method = TRUE,
      #log.rank.weights = "n",
      conf.int = F,         # show confidence intervals for point estimaes of survival curves.
      xlab = xlabel,
      ylab = ylabel,
      risk.table = risk_table,       # show risk table.
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = F, # show bars instead of names in text annotations in legend of risk table

      #ncensor.plot = TRUE
      # cumevents=T
      # cumcensor=T

      surv.median.line="hv"
      # tables.theme = theme_cleantable(),
      # ggtheme = theme_bw() # Change ggplot2 theme
    )
    p <- p + ggplot2::guides(colour = ggplot2::guide_legend(nrow = legend_rows))
    p
  }

}
