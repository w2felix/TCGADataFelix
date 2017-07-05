#' Title
#'
#' @param x a character vector with 1 or more gene names
#' @param Eset An Expression Set
#' @param additional covariate to subset groups
#' @param exclude Factors of the covariates to exclude
#'
#' @return A List of Hazard Ratios of a list of genes
#' @export
#'
#' @examples
#' \dontrun{
#' hazard_list(x = c("FOXA2"), Eset = Eset, additional = "pathologic_stage", exclude = c("Stage IV")
#' }

hazard_list <- function (x,
                         Eset,
                         additional,
                         exclude) {

  # x = gene = e.g.: "TSPAN8"
  # if more that one gene in genelist -> add them all to the data frame and make an average out of them?
  #gene <- c("MYC", "TSPAN8")

  gene <- x

  time <- Biobase::pData(Eset)$X_OS
  event <- Biobase::pData(Eset)$X_OS_IND

  coxph1 <- NULL
  for(i in 1:length(gene)){
    z <- data.frame(time = time, event = event)

    genename <- as.character(gene[i])
    if(genename %in% rownames(Biobase::exprs(Eset))){

      z[,3] <- Biobase::exprs(Eset)[genename,]

      # z-score:
      Z_score_1 <- z[,3] - mean(z[,3])
      Z_score <- Z_score_1 / stats::sd(z[,3])
      z[,3] <- Z_score

      if (missing(additional)){
        z$additional <- 0
      } else if (is.na(additional)){
        z$additional <- 0
      } else if (additional==""){
        z$additional <- 0
      } else if (!additional %in% colnames(Biobase::pData(Eset))){
        stop(paste(additional, "not in phenotype list"))
      } else {
        z$additional <- Biobase::pData(Eset)[,additional]
      }

      z <- z[!is.na(z$event) & !is.na(z$time) & !is.na(z$additional) & !z$additional=="", ]

      if (z$additional[1]!=0 && is.factor(z$additional)){
        z$additional <- droplevels(z$additional)
      }

      if (!missing(exclude)){
        if(length(exclude)>1){
          for(i in 1:length(exclude)){
            z <- z[!z$additional==exclude[i], ]
          }
        } else {
          z <- z[!z$additional==exclude, ]
        }
      }

      z <- z[!is.na(z[,3]), ]
      z <- z[order(z[,3]), ]

      if(nrow(z)!=0){
        coeff <- stats::coef(summary(survival::coxph(survival::Surv(z$time, z$event) ~ z[,3])))
        coxph1 <- rbind(coxph1, coeff)
      } else {
        coxph1 <- rbind(coxph1, NA)
      }
    } else {
      coxph1 <- rbind(coxph1, NA)
    }
  }
  rownames(coxph1) <- gene

  return(coxph1)
}
