Hazard_list <- function (x,
                         time,
                         event,
                         Eset,
                         additional,
                         exclude)
  {
  # x = gene = e.g.: "TSPAN8"
  # time = pData(kikaEset)$X_OS
  # event = pData(kikaEset)$X_OS_IND

  # if more that one gene in genelist -> add them all to the data frame and make an average out of them?
  #gene <- c("MYC", "TSPAN8")
  #gene <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")
  #length(gene)
  #gene = count_exp[,1]

  geneset <- Eset
  gene <- x

  coxph1 <- NULL
  for(i in 1:length(gene)){
    z <- data.frame(time = time, event = event)

    genename <- as.character(gene[i])
    if(genename %in% rownames(exprs(geneset))){

      z[,3] <- exprs(geneset)[genename,]

      # z-score:
      Z_score_1 <- z[,3] - mean(z[,3])
      Z_score <- Z_score_1 / sd(z[,3])
      z[,3] <- Z_score

      if (missing(additional)){
        z$additional <- 0
      } else if (is.na(additional)){
        z$additional <- 0
      } else if (additional==""){
        z$additional <- 0
      } else if (!additional %in% colnames(pData(geneset))){
        stop(paste(additional, "not in phenotype list"))
      } else {
        z$additional <- pData(geneset)[,additional]
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
        coeff <- coef(summary(coxph(Surv(z$time, z$event) ~ z[,3])))
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
