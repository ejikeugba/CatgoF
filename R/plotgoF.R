#' Plot Function for Objects of Class goF
#' @description A Plot function for goodness of fit measures, it returns plots of pseudo R2's & information weights of categorical response models.
#' @usage plotgoF(obj, type)
#' @param obj saved object of class goF
#' @param type respectively selects either 'ps' or "wt" for the pseudo R2 and the information weight.
#' @author Ejike R. Ugba
#' @examples
#'
#' require("ordinal")
#' set.seed(63)
#'
#' x1  <- rnorm(200,0,1); x2 <- rnorm(200,2,3)
#' x3  <- rnorm(200,5,2); x4 <- runif(200,2,5)
#' cnt <- 0.2*x1 - 0.5*x2 + 0.4*x3 - 0.3*x4 + 2 + rnorm(200, 0, 1)
#'
#' qnt <- quantile(cnt)
#' bdv <- cut(cnt, breaks = c(-Inf,qnt[3],Inf), include.lowest=TRUE,labels = c(0,1))
#' pdv <- as.factor(cut(cnt, breaks=quantile(cnt), include.lowest=TRUE,
#'                  labels=c(1, 2, 3, 4), ordered=TRUE))
#'
#' q1 <- clm2(pdv ~ x1 + x3, link = "cloglog")
#' q2 <- clm2(pdv ~ x1 + x2 + x4, link = "cauchit")
#' q3 <- clm2(pdv ~ x1 + x2, link = "logistic")
#' q4 <- clm2(pdv ~ x1 + x2 + x4, link = "probit")
#' q5 <- clm2(pdv ~ x1 + x2, link = "logistic")
#' qlst <- list(q1,q2,q3,q4,q5)
#'
#' ## AIC weights
#' jr <- goFit(qlst, crit = "AIC", display = "all")
#' plotgoF(jr, type = "wt")
#'
#' ## Pseudo R2
#' kr <- glm(bdv ~ x1 + x2, family = binomial())
#' sr <- goFit(kr, crit = "AICc", display = "all")
#' plotgoF(sr)
#' @seealso \code{\link{goFit}}
#' @export
#'
plotgoF <- function(obj, type=NULL){

  mfun <- function(object, type){
    incall <- mget(ls())
    incall
  }
  ss <- mfun(obj, type)
  ob  <- ss$object
  tp <- ss$type

  if(all(tp == c(NULL, "ps")) == TRUE){

    if(class(ob) != "goF") stop("object is not of class goF")
    par(mar=c(5.1, 8.5, 4.5, 2.1))
    ps <- as.matrix((ob$Pseudo.R2)[,1])
    x <- ps[ order(row.names(ps), decreasing = TRUE),]

    plt <- barplot(x, horiz=TRUE, las=1, col= "#0000FFFF", xlab="value", cex.axis = 1,
                   cex.name = 1, cex.lab=1, main=expression("Pseudo"~R^2),
                   font.main=1 ,xlim=c(0,1),space = 1)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
  } else

    if(type=="wt"){

      if(class(ob) != "goF") stop("object not of class goF")
      if (is.null(nrow(ob$Info.Criterion)) ==T) stop("two or more models required for comparison")

      rw <- rownames(ob$Info.Criterion)
      cl <- colnames(ob$Info.Criterion)[4]
      vv <- as.matrix(ob$Info.Criterion[,6])

      if (nrow(ob$Info.Criterion) > 10){
        tp <- (ob$Info.Criterion)[c(1:10), ]
        rw <- rownames(tp)
        cl <- colnames(tp)[4]
        vv  <- as.matrix(tp[,6])
      }
      if (cl=="AIC")  hh <- "AIC weights"
      if (cl=="AICc") hh <- "AICc weights"
      if (cl=="BIC")  hh <- "BIC weights"

      x <- c(vv)
      if (nrow(ob$Info.Criterion) ==2){
        bp <- barplot(x, ylab="value", xaxt="n",xlab="Model",col="#0000FFFF",
                      font.main=1 , ylim = c(0,max(x)), space=3, main=hh)
        axis(1, at = bp, labels = rw, cex = 2,las=2)
      } else {
        bp <- barplot(x, ylab="value", xaxt="n",xlab="Model",col="#0000FFFF",
                      font.main=1 , ylim = c(0,max(x)), space=1, main=hh)
        axis(1, at = bp, labels = rw, cex = 2,las=2)
      }

    } else {
      message("Error: wrong type input! should be either 'ps' or 'wt'")
    }
}
