#' Measures of Fit for Categorical Response Model
#' @description this function performs the goodness-of-fit of categorical outcome models, it returns valuable information measures and pseudo R2s of fitted models.
#' @usage goFit(obj, crit, display, ...)
#' @param obj single, multiple or list of model objects
#' @param crit  information criterion type on which model comparison is based.
#' @param display  chooses to display either 'best', 'top' or 'all' result
#' @param ... .
#' @details \code{goFit} supports objects of class glm, vglm, clm2, polr and multinom. In addition to computing several pseudo-R2 measures relevant to the analysis of categorical outcome models, \code{goFit} also conducts model comparison given two or more models in the function. It provides an automatically ranked information criterion table with the best model on top. The initial object names supplied to the function are replaced with an ordered internally generated numbering as follows: (fit1, fit2, ...) with respect to the ordering of objects in the function.
#' @import MASS, VGAM, ordinal, nnet, graphics, methods, stats, utils
#' @return \item{Categ}{Number of response category}
#' @return \item{logLik}{Log-likelihood of fitted model}
#' @return \item{Param}{Number of parameters in model including intercept}
#' @return \item{AIC}{Akaike Information Criterion}
#' @return \item{AICc}{Akaike Information Criterion for small sample}
#' @return \item{BIC}{Bayesian Information Criterion}
#' @return \item{AIC.diff}{AIC difference, same also for AICc.diff and BIC.diff}
#' @return \item{Akaike.wt}{AIC weight of fitted model}
#' @return \item{Schwarz.wt}{BIC weight of fitted model}
#' @return \item{ER}{Evidence ratio of fitted model}
#' @return \item{Pseudo.R2}{Pseudo R2 of fitted model}
#' @author Ejike R. Ugba
#' @references Long, J.S. (1997). \emph{Regression Models for Categorical and Limited Dependent Variables}. California: Sage Publications.
#' @references Burnham, K. P., & Anderson, D. R. (2002). \emph{Model Selection and Multimodel Inference: A Practical Information-Theoretical Approach}. Springer.
#' @references Ugba, E. R. & Gertheiss, J. (2017). A Generalized Likelihood Ratio Index for Discrete Limited Dependent Variable Models — with Emphasis on the Cumulative Link Model. (preprint)
#' @seealso \code{\link{plotgoF}}
#' @examples
#'
#' require("MASS") ## for polr method
#' set.seed(63)
#'
#' x1  <- rnorm(200,0,1); x2 <- rnorm(200,2,3)
#' x3  <- rnorm(200,5,2); x4 <- runif(200,2,5)
#' cnt <- 0.2*x1 - 0.5*x2 + 0.4*x3 - 0.3*x4 + 2 + rnorm(200, 0, 1)
#'
#' pdv <- as.factor(cut(cnt, breaks=quantile(cnt), include.lowest=TRUE,
#'                      labels=c(1, 2, 3, 4)))
#' qnt <- quantile(cnt)
#' bdv <- cut(cnt, breaks = c(-Inf,qnt[3],Inf), include.lowest=TRUE, labels = c(0,1))
#'
#' ### Binary & polytmous outcome models
#'
#' m1 <-  glm(bdv ~ x1 + x2, family = binomial(link = "probit"))
#' m2 <- polr(pdv ~ x1 * x2 + x3, method = c("logistic"))
#' m3 <- polr(pdv ~ x1 + x2 * x3, method = c("logistic"))
#' m4 <- polr(pdv ~ x1 + x2 + x3, method = c("logistic"))
#' m5 <- polr(pdv ~ x1 * x2 * x3, method = c("logistic"))
#' m6 <- polr(pdv ~ x1 * x2 + x3 * x4, method = c("logistic"))
#'
#' mlst <- list(m2, m3, m4, m5, m6)
#'
#' ### single object
#' goFit(m1)
#'
#' ### multiple object
#' goFit(m2,m3,m4,crit="AICc",display="top")
#'
#' ### list of objects
#' goFit(mlst, crit = "BIC",display="all")
#' @export
#'

goFit <- function(obj,
                  crit = "AICc",
                  display = "best",
                  ...){
  mobj <- function(obj, ...) c(as.list(environment()), list(...))
  fn <- function(obj, ...){
    incall <- mget(ls())
    incall
  }

  ### Compute variance & probabilities
  lginfo  <- function(fit, ...){

    ## glm object
    if (((is(fit)=="glm")[1]) == TRUE){
      ft <- cbind(as.numeric_version(fit$model[,1]))
      tr <- as.numeric(fit$fitted.values)
      lp <- as.numeric(fit$linear.predictors)
      comp <- as.matrix.data.frame(cbind(ft, tr, lp))

      obs   <- length(fit$y)
      npar  <- length(fit$coef)
      grp1  <- subset(comp, comp[,1] < 1)
      grp2  <- subset(comp, comp[,1] > 0)
      vrianc <- var(comp[,3]) / (obs/(obs - 1))

      if (fit$family$link == "logit") err = ((pi)^{2})/3
      if (fit$family$link == "probit") err = 1
      if (fit$family$link == "cloglog") err = ((pi)^{2})/6
      if (all(fit$family$link != (c("logit", "probit","cloglog")))) err = NA

      if (all(class(fit$model[,1]) == "factor" |
              class(fit$model[,1]) == "ordered") == FALSE){
        comp <- diag(NA,3)
      }
    }

    ## vglm object
    if (((is(fit) == "vglm")[1])==TRUE){
      if ((length(table(VGAM::model.frame(fit)[,1]))==2)){


        ft <- cbind(as.numeric_version(VGAM::model.frame(fit)[,1]))
        tr <- as.numeric(VGAM::fitted(fit)[,2])
        lp <- as.numeric(VGAM::predictors(fit))
        comp <- as.matrix.data.frame(cbind(ft, tr, lp))

        obs   <- VGAM::nobs(fit) ;npar <- VGAM::nparam(fit)
        grp1  <- subset(comp, comp[,1] < 1)
        grp2  <- subset(comp, comp[,1] > 0)
        vrianc <- as.numeric(var(comp[,3]) / (obs/(obs - 1)))

        if (VGAM::linkfun(fit)=="logit") err = ((pi)^{2})/3
        if (VGAM::linkfun(fit)=="probit") err = 1
        if (VGAM::linkfun(fit)=="cloglog") err = ((pi)^{2})/6
        if (all(VGAM::linkfun(fit) != (c("logit", "probit","cloglog")))) err = NA
      }
      else obs <- VGAM::nobs(fit)

      if (all(class(VGAM::model.frame(fit)[,1])=="factor" |
              class(VGAM::model.frame(fit)[,1])=="ordered") == FALSE){
        comp <- diag(NA,3)
      }
    }

    ## clm2 object
    if ((is(fit)=="clm2")[1]  == TRUE ){
      if (length(fit$lev) == 2){
        ft <- cbind(as.numeric_version(fit$location[,1]),fit$fitted.values)
        fv <- abs(cbind((1 - unlist(ft[,1])) - (unlist(ft[,2]))))
        tr <- as.numeric_version(fit$location[,1])

        if ((length(fit$beta) == dim((fit$location[,-1]))[2])==TRUE){
          lp <- -(fit$Alpha) + as.matrix(fit$location[,-1])%*% fit$beta
          comp <- as.matrix.data.frame(cbind(tr, fv, lp))
        } else comp <- diag(NA,3)

        obs  <- fit$nobs
        npar <- length(fit$coefficients)
        grp1 <- subset(comp, comp[,1] < 1)
        grp2 <- subset(comp, comp[,1] > 0)
        vrianc <- (var(comp[,3]) / (obs/(obs - 1)))

        if (fit$link == "logistic") err = ((pi)^{2})/3
        if (fit$link == "probit") err = 1
        if (fit$link == "cloglog") err = ((pi)^{2})/6
        if (all(fit$link != (c("logistic", "probit","cloglog")))) err = NA
      }else obs <- fit$nobs
    }

    ## multinom object
    if (((is(fit)=="multinom")[1]) == TRUE){
      if (length(table(fit$model[,1])) == 2){
        ft <- cbind(as.numeric_version(fit$model[,1]),fit$fitted.values)
        fv <- abs(cbind((1 - unlist(ft[,1])) - (unlist(ft[,2]))))
        tr <- as.numeric_version(fit$model[,1])

        if ((length(coef(fit)[-1]) == dim((fit$model[,-1]))[2])==TRUE){
          lp <- -(coef(fit)[1]) + (as.matrix(fit$model[,-1]) %*% coef(fit)[-1])
          comp <- as.matrix.data.frame(cbind(tr, fv, lp))
        } else comp <- diag(NA,3)

        obs  <- nrow(fit$residuals)
        npar <- length(fit$vcoefnames)
        grp1 <- subset(comp, comp[,1] < 1)
        grp2 <- subset(comp, comp[,1] > 0)
        vrianc <- NA
        err  <- NA
      } else obs <- nrow(fit$residuals)
    }

    ## polr object
    if (((is(fit)=="polr")[1]) == TRUE){
      obs  <- fit$nobs
      npar <- (length(fit$zeta) + length(fit$coeff))
    }

    ### Pseudo R2 measures
    bindex <- data.frame(matrix(nrow = 1,ncol = 9))
    names(bindex)=c("Cox & Snell",
                    "Nagelkerke",
                    "McFadden",
                    "Aldrich & Nelson",
                    "Veal & Zimmermann",
                    "Ugba & Gertheiss",
                    "McKlevey & Zavoina",
                    "Efron",
                    "Tjur")
    cindex <- bindex


    LLf  <- LLik(fit)[1] ;LL0  <- LLik(fit)[2] ;LRT  <- (-2) * (LL0 - LLf)

    if (is.na(LLf)==TRUE | is.na(LL0)==TRUE) stop(" model log-likelihood missing")

    if (LLik(fit)[3] == 2){
      ncat <- LLik(fit)[3]
      bindex[1] <- as.vector(1 - exp((2/obs) * (LL0 - LLf)))		                                            # Cox&Snell
      bindex[2] <- as.vector((1 - exp((2/obs) * (LL0 - LLf))) / (1 - exp(LL0)^{2/obs}))                     # Nagelkerke
      bindex[3] <- as.vector(1 - (LLf/LL0))                                                                 # McFadden
      bindex[4] <- as.vector(LRT / (LRT + obs))	                                                            # Aldrich $ Nelson
      bindex[5] <- as.vector((LRT / (LRT + obs)) * ((obs - (2 * LL0)) / (-2 * LL0)))                        # Veal & Zimmermann
      bindex[6] <- as.vector(1 - (LLf / LL0)^{ncat})                                                        # Ugba & Gertheis
      bindex[7] <- as.vector(vrianc / (vrianc + err))                                                       # Mcklevey & Zavoina
      bindex[8] <- as.vector(1 - (sum((comp[,1] - comp[,2])^{2}) / sum((comp[,1] - mean(comp[,1]))^{2})))   # Effron
      bindex[9] <- as.vector(mean(grp2[,2]) - mean(grp1[,2]))                                               # Tjur
    }
    if (LLik(fit)[3] > 2){
      ncat <- LLik(fit)[3]
      cindex[1] <- as.vector(1 - exp((2/obs) * (LL0 - LLf)))
      cindex[2] <- as.vector((1 - exp((2/obs)*(LL0 - LLf)))/(1 - exp(LL0)^(2/obs)))
      cindex[3] <- as.vector(1 - (LLf / LL0))
      cindex[4] <- as.vector(LRT / (LRT + obs))
      cindex[5] <- as.vector((LRT / (LRT + obs)) * ((obs - (2 * LL0)) / (-2 * LL0)))
      cindex[6] <- as.vector(1 - (LLf / LL0)^{sqrt(2 * ncat)})
      cindex[c(7:9)] <- as.numeric(NA)
    }
    if (LLik(fit)[3] == 2) calc <- bindex
    if (LLik(fit)[3]  > 2) calc <- cindex

    calc$ncat <- ncat
    calc$LLf  <- LLf
    calc$npar <- LLik(fit)[4]
    calc$BIC  <- LLik(fit)[7]
    calc$AIC  <- LLik(fit)[5]
    calc$AICc <- LLik(fit)[6]
    core <- t(round(calc,3))
    return(core)
  }

  mobject <- mobj(obj, ...)
  if((length(is(mobject[[1]])=="fit")==3)==T){
    mobject <- mobject
  }
  if((length(is(mobject[[1]])=="fit")==2)==T){
    mobject <- unlist(mobj(obj, ...), recursive = FALSE)
  }
  hr <- fn(obj, ...)
  if(is.list(hr$crit)|is.list(hr$display) == TRUE) stop("arguments criterion and output missing")

  collect <- list()
  if ((class(mobject)=="list")==TRUE){
    if (all(class(mobject[[1]])[1] != (c("glm", "vglm", "clm2", "polr", "multinom"))) == TRUE ){
      stop("Object of class not supported")
    } else {
      collect <- sapply(mobject,lginfo)
    }
  }
  if ((class(mobject) !="list") ==TRUE){
    if (all(class(mobject)[1] != (c("glm", "vglm", "clm2", "polr", "multinom"))) == TRUE ){
      stop("Object of class not supported")
    } else {
      collect <- sapply(mobject,lginfo)
    }
  }

  if ((length(mobject) >= 2)==TRUE){
    hh <- matrix(NA, nrow = length(mobject), ncol = 1)
    for(r in 1:length(mobject)){
      hh[r] <- class(mobject[[r]])[1]
    }
    if (all(hh == hh[1])==FALSE) warning("different objects of class compared")
  }
  if (is.matrix(collect) == TRUE  ){
    hr <- collect
    r  <- ncol(hr)
    rnames <- c("Cox & Snell",
                "Nagelkerke",
                "McFadden",
                "Aldrich & Nelson",
                "Veal & Zimmermann",
                "Ugba & Gertheiss",
                "McKlevey & Zavoina",
                "Efron",
                "Tjur",
                "Categ",
                "logLik",
                "Param",
                "BIC",
                "AIC",
                "AICc")
    fun <- function (r){
      v <- sprintf('fit%d', 1:r,r)
      return(v)
    }
    cnames <- fun(r)
    colnames(hr) <- c(cnames)
    rownames(hr) <- c(rnames)
  }
  if (is.list(collect) == TRUE){
    count  <- sapply(collect, nrow)
    rnames <- lapply(collect[which.max(count)], rownames)
    hr  <- sapply(collect, '[', seq(max(sapply(collect, length))))
    rownames(hr) <- c(rnames[[1]])
    r   <- length(collect)
    fun <- function (r){
      v <- sprintf('fit:%d', 1:r,r)
      return(v)
    }
    cnames <- fun(r)
    colnames(hr) <- c(cnames)
  }

  ### Akaike weights, Schwarz weights & Evidence ratio

  Info <- t(hr[c(10:15), ])
  if(length(crit) > 1){
    Info = Info[ ,-c(4,5)]
  } else {
    Info <- switch(crit,
                   'AICc'= Info[ ,-c(4,5)],
                   'AIC' = Info[ ,-c(4,6)],
                   'BIC' = Info[ ,-c(5,6)])
  }
  if (ncol(as.matrix(Info)) > 1){
    delta <- matrix(NA, nrow = 1, ncol = nrow(Info))
    for(i in 1:length(delta)){
      delta[i] <- Info[,4][i] - min(Info[,4])
    }

    weight <- matrix(NA, nrow = 1, ncol = nrow(Info))
    for(i in 1:length(weight)){
      weight[i] <- exp(-delta[i]/2) / sum(exp(-delta/2))
    }
    weight <- round(weight,3)
    criterion   <- cbind(Info,t(delta),t(weight))

    if (crit == "AICc") colnames(criterion) = c("categ","logLik","Param","AICc","AICc.df","Akaike.wt")
    if (crit == "AIC")  colnames(criterion) = c("categ","logLik","Param","AIC","AIC.df","Akaike.wt")
    if (crit == "BIC")  colnames(criterion) = c("categ","logLik","Param","BIC","BIC.df","Schwarz.wt")

    criterion <- criterion[ order(-criterion[,5], decreasing = TRUE), ]
    criterion <- round(criterion,2)

    EvdR <- matrix(NA, nrow = 1, ncol = length(criterion[,5]) - 1)
    for(i in 1:length(EvdR)){
      EvdR[i] <- exp(-criterion[,5][1]/2) / exp(-criterion[,5][i + 1]/2)
    }
    EvdR <- round(t(EvdR),2)
    rnam <- rownames(criterion)
    h    <- cbind(rnam[1],rnam[-1])

    v <- matrix(NA,nrow = 1, ncol = nrow(h))
    for (r in 1:nrow(h)){
      v[r] <- sprintf('%s|%s', h[r,1],h[r,2])
    }
    v <- t(v)
    EvdR <- data.frame(Model = v, ER = EvdR)

    Pseu <- as.matrix(hr[-c(10:15), ])
    best  <- rownames(criterion)[1]
    selc  <- Pseu[ , best, drop = FALSE]
    Pseu <- cbind(selc[complete.cases(selc), ])
    colnames(Pseu) <- colnames( selc )
  } else EvdR <- NULL

  Ps <- as.matrix(hr[-c(10:15), 1])
  Ps <- cbind(Ps[complete.cases(Ps), ])
  colnames(Ps) <- ("value")

  if (length(mobject) == 1){
    infocrt <- hr[c(10:15), ]
    core1 <- list(Info.Criterion = infocrt, Pseudo.R2 = Ps)
  }
  if (length(mobject) >1){
    jr <- list(Info.Criterion = criterion[1,],
               Evidence.Ratio = EvdR[1,],
               Pseudo.R2 = Pseu)

    jk <- list(Info.Criterion = head(criterion,2),
               Evidence.Ratio = head(EvdR,1),
               Pseudo.R2 = Pseu)

    jw <- list(Info.Criterion = criterion,
               Evidence.Ratio = EvdR,
               Pseudo.R2 = Pseu)

    core2 <- switch(display,
                    'best'= jr,
                    'top' = jk,
                    'all' = jw)
  }
  if (length(mobject)==1) res <- core1
  if (length(mobject) >1) res <- core2

  class(res) <- c("goF")
  return(res)
}
