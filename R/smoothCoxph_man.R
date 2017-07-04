#' Create a correlation of gene expression to log hazard ratio
#'
#' @param x 1 or multiple Genes, as a character vector list, to analyze
#' @param Eset The Expression Set, containing the expression data in columns and the survival indicators in the columns X_OS and X_OS_IND
#' @param exclude_values c("where to look for exclusion", "exclude group 1", "exclude group 2", ...)
#' @param xlimit Manually define the x-axis limit
#' @param ylimit Manually define the y-axis limit
#' @param xlabel Manually define the x-axis label
#' @param ylabel manually define the y-axis label
#' @param logrisk don't exactly know... need to check
#' @param average values: `mean` or `median`, defines how >1 genes are averaged
#' @param ... additional variables that can be added
#'
#' @return Plots the log Hazard ratio in correlation to its relative gene expression level using ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' smoothCoxph_man(x = c("BAMBI", "CD44"), Eset=Eset,
#' exclude_values=c("pathologic_stage", "Stage II", "Stage I", "Stage IV"),
#' logrisk=T, average="median")
#' }
#'
smoothCoxph_man <- function (x, Eset, exclude_values,
                               xlimit, ylimit, xlabel, ylabel,
                               logrisk = TRUE, average = "median", ...) {

  time <- Biobase::pData(Eset)$X_OS
  event <- Biobase::pData(Eset)$X_OS_IND
  gene <- x

  # for package reasons: Define variables first:
  hazard_ratio <- NULL
  upper95 <- NULL
  lower95 <- NULL



  if(length(gene)>1){
    z <- data.frame(time = time, event = event)
    z[,gene] <- t(Biobase::exprs(Eset)[gene,])

    # z-score the tcga data

    for(i in 3:length(z)){
      Z_score_1 <- z[,i] - mean(z[,i])
      Z_score <- Z_score_1 / stats::sd(z[,i])
      z[,i] <- Z_score
    }

    df_matrix <- matrix(t(z[,gene]), ncol = length(gene), byrow=TRUE)

    z$Median <- Biobase::rowMedians(df_matrix)
    z$Mean <- base::rowMeans(df_matrix)

  } else {
    x <- Biobase::exprs(Eset)[gene[1],]
    z <- data.frame(time = time, event = event, x = x)

    # z-score the data
    for(i in 3:length(z)){
      Z_score_1 <- z[,i] - mean(z[,i])
      Z_score <- Z_score_1 / stats::sd(z[,i])
      z[,i] <- Z_score
      x <- Z_score
    }
  }

  if (!missing(exclude_values)){
    if (exclude_values[1] %in% colnames(Biobase::pData(Eset))){
      z$exclude_values <- Biobase::pData(Eset)[,exclude_values[1]]
    } else {
      stop(paste(exclude_values[1],"does not exist in the phenoData"))
    }
    if(length(exclude_values)<2){
      stop(paste("Please state values to be excluded in",exclude_values[1]))
    } else {
      for(i in 2:length(exclude_values)){
        z <- z[!z$exclude_values==exclude_values[i], ]
      }
    }
  }

  z <- z[!is.na(z$event) & !is.na(z$time), ]

  if(length(gene)>1){
    for(i in 1:length(gene)){
      z <- z[!is.na(z[,gene[i]]), ]
    }
    if(average=="mean"){
      z <- z[order(z$Mean), ]
      coxph1 <- survival::coxph(survival::Surv(time, event) ~ survival::pspline(Mean, df = 4), data = z)
    } else {
      z <- z[order(z$Median), ]
      coxph1 <- survival::coxph(survival::Surv(time, event) ~ survival::pspline(Median, df = 4), data = z)
    }
  } else {
    z <- z[!is.na(z$x), ]
    z <- z[order(z$x), ]
    coxph1 <- survival::coxph(survival::Surv(time, event) ~ survival::pspline(x, df = 4), data = z)
  }

  if (logrisk) {
    pred <- stats::predict(coxph1, type = "lp", se.fit = TRUE)
    y <- cbind(pred$fit, pred$fit - 1.96 * pred$se.fit, pred$fit +
                 1.96 * pred$se.fit)
  } else {
    pred <- stats::predict(coxph1, type = "risk", se.fit = TRUE)
    y <- cbind(log(pred$fit), log(pred$fit - 1.96 * pred$se.fit),
               log(pred$fit + 1.96 * pred$se.fit))
  }
  if (missing(xlimit)){
    if(length(gene)>1){
      if(average=="mean"){
        xlimit <- range(z$Mean)
      } else {
        xlimit <- range(z$Median)
      }
    } else {
      xlimit <- range(x)
    }
  }
  if (missing(ylimit)) ylimit <- c(min(y[, 1], stats::median(y[, 2]), na.rm = T), max(y[, 1], stats::median(y[, 3]), na.rm = T)) * 1.5
  if (missing(xlabel)) xlabel <- "Relative gene expression level"
  if (missing(ylabel)) ylabel <- "log Hazard Ratio"

  if(length(gene)>1){
    if(average=="mean"){
      coxph2 <- survival::coxph(survival::Surv(time, event) ~ Mean, data = z)
      data_smooth <- data.frame(expression = z$Mean[1:length(pred$fit)],
                                hazard_ratio = y[, 1],
                                lower95 = y[, 2],
                                upper95 = y[, 3])
    } else {
      coxph2 <- survival::coxph(survival::Surv(time, event) ~ Median, data = z)
      data_smooth <- data.frame(expression = z$Median[1:length(pred$fit)],
                                hazard_ratio = y[, 1],
                                lower95 = y[, 2],
                                upper95 = y[, 3])
    }

  } else {
    coxph2 <- survival::coxph(survival::Surv(time, event) ~ x, data = z)
    data_smooth <- data.frame(expression = z$x[1:length(pred$fit)],
                              hazard_ratio = y[, 1],
                              lower95 = y[, 2],
                              upper95 = y[, 3])
  }

  coeffs <- stats::coef(summary(coxph2))
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

  Plot_Title <- paste("Correlation of", gene_name, "expression levels\nto log hazard ratio")

  ggplot2::ggplot(data_smooth, ggplot2::aes(x = expression)) +
    ggplot2::ggtitle(Plot_Title) +
    ggplot2::theme(plot.title = ggplot2::element_text(size=11,
                                    face="bold",
                                    margin = ggplot2::margin(10, 0, 20, 0))) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept=0, size=1, color = "grey") +
    ggplot2::geom_line(ggplot2::aes(y = hazard_ratio), color = "dark blue", lwd = 1) +
    ggplot2::geom_line(ggplot2::aes(y = upper95), color = "dark grey", linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = lower95), color = "dark grey", linetype = "dashed") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower95, ymax = upper95), fill = "steelblue2", color=NA, alpha = 0.1) +
    ggplot2::ylim(ylimit) +
    ggplot2::xlim(xlimit) +
    ggplot2::ylab(label=ylabel) +
    ggplot2::xlab(label=xlabel) +
    ggplot2::annotate("text", x=xlimit[2], y=(ylimit[2]-(ylimit[2]-ylimit[1])*0.95), hjust = 1, label=paste("HR =", round(coeffs[,2], digits=4), "+/-", round(coeffs[,3], digits=4))) +
    ggplot2::annotate("text", x=xlimit[2], y=(ylimit[2]-(ylimit[2]-ylimit[1])*0.99), hjust = 1, label=paste("p <", ifelse(coeffs[,5]>0.05, "NS", ifelse(coeffs[,5]>0.01, "0.05", ifelse(coeffs[,5]>0.005, "0.01", "0.005")))))
}
