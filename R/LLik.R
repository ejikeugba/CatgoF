#' Loglikelihoods of Categorical Response Model
#' @param fit model objects
#'
LLik <- function(fit){

  ## clm2 LogL
  if (((is(fit)=="clm2")[1]) == TRUE){
    Lclm2 <- function(fit){
      LLf <- fit$logLik
      LL0 <- (clm2(as.factor(fit$y) ~ 1))$logLik
      return(c(LLf, LL0))
    }
    f1   <- Lclm2(fit)
    ncat <- length(fit$lev)
    npar <- length(fit$coefficients)
    AC   <- (-2)*fit$logLik + (2*npar)
    ACc  <- (-2)*fit$logLik + (2)*(ncol(fit$location) - 1) + ((2*npar*(npar + 1))/(fit$nobs - npar - 1))
    BC   <- (-2)*fit$logLik + (length(fit$coefficients) * log(fit$nobs))
  }

  ## polr LogL
  if (((is(fit) == "polr")[1]) == TRUE){
    Lpolr <- function(fit){
      LLf <- (fit$deviance)/(-2)
      LL0 <- ((polr(as.factor(fit$model[,1]) ~ 1,method = c("logistic")))$deviance)/(-2)
      return(c(LLf, LL0))
    }
    f2   <- Lpolr(fit)
    ncat <- length(fit$lev)
    npar <- (length(fit$zeta) + length(fit$coeff))
    AC   <- fit$deviance + (2)*npar
    ACc  <- fit$deviance + (2)*(ncol(fit$model) - 1) + ((2*npar*(npar + 1))/(fit$nobs - npar - 1))
    BC   <- fit$deviance + ((length(fit$zeta) + length(fit$coeff)) * log(fit$nobs))
  }

  ## multinom LogL
  if (((is(fit)=="multinom")[1]) == TRUE){
    Lmultinom <- function(fit){
      if (is.null(fit$model)== TRUE){
        stop("argument model=T missing in multinom")
      }
      LLf  <- logLik(fit)
      fnt  <- function(fit){
        y  <- as.matrix(fit$model[,1])
        if(ncol(y) > 1) {y <- as.ordered(apply(y, 1, function(x) which(x==1)))}
        h <- multinom(y ~ 1, trace = FALSE)
        LL0  <- logLik(h)
        return(LL0)
      }
      LL0  <- fnt(fit)
      return(c(LLf, LL0))
    }
    f3   <- Lmultinom(fit)
    ncat <- length(fit$lev)
    npar <- length(fit$vcoefnames) - 1
    AC   <- fit$AIC
    ACc  <- fit$AIC + ((2*npar*(npar + 1))/(nrow(fit$residuals) - npar - 1))
    BC   <- fit$deviance + (length(fit$vcoefnames) * log(nrow(fit$residuals)))
  }

  ## vglm LogL
  if (((is(fit)=="vglm")[1]) == TRUE){

    Lvglm <- function(fit){
      if (all(familyname(fit) != (c("acat", "cumulative", "propodds", "cratio",  "sratio",
                                    "multinomial","brat","bratt","ordpoisson"))) == TRUE ){
        stop("family not supported")
      }
      dvar <- VGAM::model.frame(fit)[,1]

      if ((familyname(fit)=="acat")==TRUE)        LL0 <- logLik(vglm(dvar ~ 1, family=acat()))
      if ((familyname(fit)=="cumulative")==TRUE)  LL0 <- logLik(vglm(dvar ~ 1, family=cumulative()))
      if ((familyname(fit)=="propodds")==TRUE)    LL0 <- logLik(vglm(dvar ~ 1, family=propodds()))
      if ((familyname(fit)=="cratio")==TRUE)      LL0 <- logLik(vglm(dvar ~ 1, family=cratio()))
      if ((familyname(fit)=="sratio")==TRUE)      LL0 <- logLik(vglm(dvar ~ 1, family=sratio()))
      if ((familyname(fit)=="multinomial")==TRUE) LL0 <- logLik(vglm(dvar ~ 1, family=multinomial()))
      if ((familyname(fit)=="brat")==TRUE)        LL0 <- logLik(vglm(dvar ~ 1, family=brat()))
      if ((familyname(fit)=="bratt")==TRUE)       LL0 <- logLik(vglm(dvar ~ 1, family=bratt()))
      if ((familyname(fit)=="ordpoisson")==TRUE)  LL0 <- logLik(vglm(dvar ~ 1, family=ordpoisson()))

      LLf <- logLik(fit)
      return(c(LLf, LL0))
    }
    f4    <- Lvglm(fit)
    ncat  <- length(table(VGAM::model.frame(fit)[,1]))
    npar  <- VGAM::nparam(fit)
    AC    <- VGAM::AIC(fit)
    ACc   <- VGAM::AIC(fit) + ((2*npar*(npar + 1))/(VGAM::nobs(fit) - npar - 1))
    BC    <- VGAM::deviance(fit) + (VGAM::nparam(fit)) * log(VGAM::nobs(fit))
  }

  ## glm LogL
  if (((is(fit)=="glm")[1]) == TRUE){
    if (length(table(fit$model[ ,1])) != 2){
      stop("binary response not used in glm")
    }
    Lglm <- function(fit){
      if (all(family(fit)[1] != (c("gaussian", "binomial", "Gamma", "poisson",
                                   "inverse.gaussian","quasi"))) == TRUE ){
        stop("family not supported")
      }
      if ((family(fit)[1]=="gaussian")==TRUE) LL0 <- as.numeric(logLik(glm((fit$y) ~ 1, family = gaussian())))
      if ((family(fit)[1]=="binomial")==TRUE) LL0 <- as.numeric(logLik(glm((fit$y) ~ 1, family = binomial())))
      if ((family(fit)[1]=="Gamma")==TRUE)    LL0 <- as.numeric(logLik(glm((fit$y) ~ 1, family = Gamma())))
      if ((family(fit)[1]=="poisson")==TRUE)  LL0 <- as.numeric(logLik(glm((fit$y) ~ 1, family = poisson())))
      if ((family(fit)[1]=="inverse.gaussian")==TRUE) LL0 <- as.numeric(logLik(glm((fit$y) ~ 1, family = inverse.gaussian())))
      if ((family(fit)[1]=="quasi")==TRUE)    LL0 <- as.numeric(logLik(glm((fit$y) ~ 1, family = quasi())))

      LLf <- logLik(fit)
      return(c(LLf, LL0))
    }
    f5   <- Lglm(fit)
    ncat <- length(table(fit$y))
    npar <- length(fit$coef)
    AC   <- fit$aic
    ACc  <- fit$aic + ((2*npar*(npar + 1))/(length(fit$y) - npar - 1))
    BC   <- fit$deviance + (length(fit$coef) * log(length(fit$y)))
  }

  if (((is(fit)=="clm2")[1])==T) return(c(f1, ncat, npar, AC, ACc, BC))
  if (((is(fit)=="polr")[1])==T) return(c(f2, ncat, npar, AC, ACc, BC))
  if (((is(fit)=="multinom")[1])==T) return(c(f3, ncat, npar, AC, ACc, BC))
  if (((is(fit)=="vglm")[1])==T) return(c(f4, ncat, npar, AC, ACc, BC))
  if (((is(fit)== "glm")[1])==T) return(c(f5, ncat, npar, AC, ACc, BC))
}
