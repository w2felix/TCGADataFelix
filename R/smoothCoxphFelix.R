smoothCoxph_Felix <- function (x, Eset, gene, exclude_values,
                               xlimit, ylimit, xlabel, ylabel,
                               logrisk = TRUE, durchschnitt = "median", zscore=T, ...) {

  # x <- dataframe with time as first column and event as second column
  # exclude_values <- c("where to look for exclusion", "exclude group 1", "exclude group 2", ...)
  # gene = c("TSPAN8", "MYC")
  # Eset <- kika_TCGA_Eset
  # time = pData(kika_TCGA_Eset)$X_OS
  # event = pData(kika_TCGA_Eset)$X_OS_IND
  # ylabel <- "log Hazard Ratio"
  # xlabel <- "Gene expression"
  # ylimit <- c(min(y[, 1], median(y[, 2]), na.rm = T), max(y[, 1], median(y[, 3]), na.rm = T)) * 1.5
  # xlimit <- range(x)
  # durchschnitt = "mean", or "median"

  # if more that one gene in genelist -> add them all to the data frame and make an average out of them?
  #gene <- "MYC"
  #gene <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")
  #length(gene)

  time <- x[,1]
  event <- x[,2]
 return(paste(time,event))
}
