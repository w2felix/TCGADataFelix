#' Title
#'
#' @param x a character vector with 1 or more gene names
#' @param Eset An Expression Set
#' @param value Defines the value to subdivide the gene expression groups.
#' Numeric: Devided into two groups
#' "q": 25% quantile, 25-75% quantile and 75% quantile
#' @param additional Additional covariate to subset groups
#' @param exclude Factors of the covariates to exclude
#' @param p.val display the p-Value on the graph
#' @param xlabel User defined x-axis label
#' @param ylabel User defined y-axis label
#' @param logrisk logrisk?!
#' @param legend_position where should the legend be? doesn't work with theme_bw
#' @param durchschnitt for more than one gene, how the value of the averaged z-score is calculated, either median or mean
#' @param optimal calculate the optimal cutpoint, will overide the value value when numeric, does not work when value = "q"
#' @param plot_cutpoint Plot the graph how the optimal cutpoint was calculated, normal survival plot will be REPLACED by surv_cutpoint: Determine the optimal cutpoint for each variable using ’maxstat’
#' @param ...
#'
#' @return A survival Estimator
#' @export
#'
#' @examples Survival_adaptable(x = c("FOXA2"), Eset = Eset, value = 0,
#' additional = "pathologic_stage", exclude = c("Stage IV"), zscore=T, p.val=TRUE,
#' xlabel="Days", legend_position = "top", durchschnitt = "median", optimal=T, plot_cutpoint=F)
Survival_adaptable <- function (x, Eset, value = 0,
                                additional, exclude,
                                p.val = FALSE, xlabel, ylabel,
                                logrisk = TRUE, legend_position,
                                durchschnitt = "mean", optimal = FALSE,
                                plot_cutpoint=FALSE, ...) {
  if (missing(x)) {
    stop("You have to define a gene you want to calculate survival")
  }

  time <- pData(Eset)$X_OS
  event <- pData(Eset)$X_OS_IND
  gene <- x

  # time <- pData(kikaEset)$X_OS
  # event <- pData(kikaEset)$X_OS_IND
  # value <- 0
  # additional <- "hemoglobin_result"
  # value seperates the groups

  # p.value only sensible for two groups to compare

  if (missing(xlabel)) xlabel <- "time"
  if (missing(ylabel)) ylabel <- "Survival"

  # get the data to the Surv variables
  # Adding a column with gene expression labels to the phenotype, updates if new values


  z <- data.frame(time = time, event = event)
  if(length(gene)>1){

    z <- data.frame(time = time, event = event)
    z[,gene] <- t(exprs(kikaEset)[gene,])
    rownames(z) <- colnames(exprs(kikaEset)[gene,])

  } else {
    if(!gene %in% rownames(exprs(kikaEset))){
      stop(paste(gene, "not in gene list"))
    }
    z[,gene] <- exprs(kikaEset)[gene,]
    rownames(z) <- colnames(exprs(kikaEset)[gene,])
  }

  if (missing(additional)){
    z$additional <- 0
  } else if (is.na(additional)){
    z$additional <- 0
  } else if (additional==""){
    z$additional <- 0
  } else if (!additional %in% colnames(pData(kikaEset))){
    stop(paste(additional, "not in phenotype list"))
  } else {
    z$additional <- pData(kikaEset)[,additional]
  }

  # Get rid of empty rows
  z <- z[!is.na(z$event) & !is.na(z$time) & !is.na(z$additional) & !z$additional=="",  ]

  if (!missing(exclude)){
    if(length(exclude)>1){
      for(i in 1:length(exclude)){
        z <- z[!z$additional==exclude[i], ]
      }
    } else {
      z <- z[!z$additional==exclude, ]
    }
  }

  if (z$additional[1]!=0 && is.factor(z$additional)){
    z$additional <- droplevels(z$additional)
  }

  if(length(gene)>1){

    ## Z-Score with more than 1 gene:

    for(i in 3:(length(z)-1)){
      Z_score_1 <- z[,i] - mean(z[,i])
      Z_score <- Z_score_1 / sd(z[,i])
      z[,i] <- Z_score
    }

    df_matrix <- matrix(t(z[,gene]), ncol = length(gene), byrow=TRUE)

    z$Median <- rowMedians(df_matrix)
    z$Mean <- rowMeans(df_matrix)

    if(value!="q"){
      if(optimal){
        value_cutpoint <- surv_cutpoint(z, time="time", event="event", variables = gene)
        value_sum <- summary(value_cutpoint)
        value <- value_sum$cutpoint
        value <- mean(value)
        if(plot_cutpoint){
          par(mfrow=c(2,1))
          return(plot(value_cutpoint, gene, palette = "npg"))
        }
      }
      z$gene <- ifelse(z$Median > value, "High Expression", "Low Expression")
    }
    if(value == "q"){
      if(optimal){
        stop("optimal cutpoint not possible if separation of patients by quantiles")
      }
      quantile_values <- quantile(z$Median, c(.25, .50, .75))
      z$gene <- ifelse(z$Median < quantile_values[1], "Lower Quantile", ifelse(z$Median < quantile_values[3], "Intermediate", "Upper Quantile"))
    }
    z$gene <- as.factor(z$gene)


  } else {
    ## Z-Score only 1 gene:
    Z_score_1 <- z[,3] - mean(z[,3])
    Z_score <- Z_score_1 / sd(z[,3])
    z[,3] <- Z_score

    if(value!="q"){
      if(optimal){
        value_cutpoint <- surv_cutpoint(z, time="time", event="event", variables = gene)
        value_sum <- summary(value_cutpoint)
        value <- value_sum$cutpoint
        if(plot_cutpoint){
          p <- plot(value_cutpoint, gene, palette = "npg")
          return(p)
        }
      }
      z$gene <- ifelse(z[,3] > value, "High Expression", "Low Expression")
    } else if(value == "q"){
      if(optimal){
        stop("optimal cutpoint not possible if separation of patients by quantiles")
      }
      quantile_values <- quantile(z[,3], c(.25, .50, .75))
      z$gene <- ifelse(z[, 3] < quantile_values[1], "Lower Quantile", ifelse(z[, 3] < quantile_values[3], "Intermediate", "Upper Quantile"))
    }
  }

  z$gene <- as.factor(z$gene)

  if(length(gene)>1){
    for(i in 1:length(gene)){
      z <- z[!is.na(z[,gene[i]]), ]
    }
    if(durchschnitt=="mean"){
      z <- z[order(z$Mean), ]
      coxph1 <- coxph(Surv(time, event) ~ Mean, data = z)
    } else {
      z <- z[order(z$Median), ]
      coxph1 <- coxph(Surv(time, event) ~ Median, data = z)
    }
  } else {
    z <- z[!is.na(z$gene), ]
    z <- z[order(z$gene), ]
    coxph1 <- coxph(Surv(time, event) ~ gene, data = z)
  }
  coeffs <- coef(summary(coxph1))

  if(p.val==TRUE){


    # The second survdiff() argument, rho, designates the weights according to ^S(t)p
    # and may be any numeric  value. The default is rho=0
    # and corresponds to the log-rank test. The Peto& Peto modi cation of the Gehan-Wilcoxon test
    # is computed using rho=1

    p_value_survdiff <- survdiff(Surv(time, event) ~ gene, data = z, rho=1)
    p.value = 1 - pchisq(p_value_survdiff$chisq, length(p_value_survdiff$n) - 1)
  }

  if (z$additional[1]!=0){

    # get Labels for legend HIER FEHLT DIE NUMMERIERUNG DER EINZELNEN KATEGORIEN; WIE VIELE N pro GRUPPE!

    z$additional <- droplevels(z$additional)
    z$gene <- droplevels(z$gene)
    numbers_quantiles <- data.frame(table(z$additional, z$gene))

    # t 1/2

    median_survival <- survfit(Surv(time, event) ~ gene + additional, data = z)
    median_survival_df <- summary(median_survival)
    #     median_survival_df$table[i,"median"]

    Nameing_vector_all <- c(paste(numbers_quantiles[,2], " - ", additional, ": ", numbers_quantiles[,1], "; n = ", numbers_quantiles[,3], "; t1/2 = ",     median_survival_df$table[,"median"], sep=""))
    Nameing_vector <- c(paste(numbers_quantiles[,2], " Expression - ", additional, ": ", numbers_quantiles[,1], sep=""))

    if(length(gene)>1){
      label_mean <- ifelse(durchschnitt=="mean","Mean","Median")
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

    #survival function:

    #fit <- survfit(Surv(time, event) ~ gene + additional, data = z,
    fit <- survfit(Surv(time, event) ~ gene + additional, data = z,
                   type="kaplan-meier",
                   conf.type="log")

    p <- ggsurvplot(
      fit,                     # survfit object with calculated statistics.
      data=z,
      surv.scale = c("percent"),
      legend = legend_position,
      legend.labs = Nameing_vector,
      legend.title = paste(gene_name, "Expression:"),
      # xlim = c(0,2000),        # present narrower X axis, but not affect survival estimates.
      break.time.by = 500,     # break X axis in time intervals by 500.

      pval = TRUE,             # show p-value of log-rank test.
      pval.method = TRUE,
      #log.rank.weights = "n",
      conf.int = F,         # show confidence intervals for point estimaes of survival curves.

      risk.table = TRUE,       # show risk table.
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = F, # show bars instead of names in text annotations in legend of risk table

      #ncensor.plot = TRUE
      # cumevents=T
      # cumcensor=T

      surv.median.line="hv",
      # tables.theme = theme_cleantable(),
      ggtheme = theme_bw() # Change ggplot2 theme
    )
    p

  } else {


    # t 1/2

    median_survival <- survfit(Surv(time, event) ~ gene, data = z)
    median_survival_df <- summary(median_survival)
    #     median_survival_df$table[i,"median"]

    z$gene <- as.factor(z$gene)
    gene_levels <- nlevels(z$gene)
    Nameing_vector <- c()
    for(i in 1:gene_levels) {
      Nameing_vector <- c(Nameing_vector, paste(levels(z$gene)[i], "; n = ", nrow(z[z$gene==levels(z$gene)[i], ]), sep=""))
    }


    if(length(gene)>1){
      label_mean <- ifelse(durchschnitt == "mean","Mean","Median")
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

    # Get p value for the difference of the gene and plot them on the graph

    fit <- survfit(Surv(time = time,
                        event = event) ~ gene,
                   data = z,
                   type="kaplan-meier",
                   conf.type="log")
    p <- ggsurvplot(
      fit,                     # survfit object with calculated statistics.
      data=z,
      surv.scale = c("percent"),
      legend = legend_position,
      legend.labs = Nameing_vector,
      legend.title = paste(gene_name, "Expression:"),
      # xlim = c(0,2000),        # present narrower X axis, but not affect survival estimates.
      break.time.by = 500,     # break X axis in time intervals by 500.

      pval = TRUE,             # show p-value of log-rank test.
      pval.method = TRUE,
      #log.rank.weights = "n",
      conf.int = F,         # show confidence intervals for point estimaes of survival curves.

      risk.table = TRUE,       # show risk table.
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = F, # show bars instead of names in text annotations in legend of risk table

      #ncensor.plot = TRUE
      # cumevents=T
      # cumcensor=T

      surv.median.line="hv",
      # tables.theme = theme_cleantable(),
      ggtheme = theme_bw() # Change ggplot2 theme
    )
    p
  }

}
