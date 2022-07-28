#' Conditional inference  
#'
#' Conditional inference for lm and glm models
#' @param object An lm() or glm() object
#' @param cond.data Optional, a dataframe for the conditioning set; set as all covariates in the lm() or glm() object formula if not provided
#' @param param Optional, a vector of coefficients to conduct conditional inference; fit all coefficients if not provided; can be a mixture of string name and index
#' @param alg Optinal, a string for name of algorithm, current options are 'loess' and 'grf'
#' @param random.seed Optional, random seed for sample splitting
#' @param other.params Optional, other parameters for the regression algorithm; can include span and degree for loess
#' @param folds Optional, a list of two folds of indices for sample splitting; can be useful to control sample splitting
#' @param verbose Optional, whether or not to print summary of inference; TRUE by default
#' @return Standard error for conditional parameter, super-population parameter, fitted empirical parameter, confidence interval for conditional parameter
#' @examples 
#' X = matrix(rnorm(1000*10), nrow=1000)
#' Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
#' Z = data.frame(X[,1:2])
#' lm.mdl = lm(Y~., data = data.frame(X))
#' cond.inf(lm.mdl, Z, c("X2",2)) 
#' 
#' X = matrix(rnorm(1000*10), nrow=1000)
#' logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
#' Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x)))
#' Z = data.frame(X[,1:2])
#' glm.mdl = glm(Y~., data = data.frame(X), family='binomial')
#' cond.inf(glm.mdl, cond.data=Z)
#' cond.inf(glm.mdl, cond.data=Z, c(1, "X1", "X2"), alg='grf')
#' 
#' @export
cond.inf <- function(object,cond.data=NULL,param=NULL,alg="loess",
                     random.seed=NULL,other.params=NULL,
                     folds=NULL,verbose=TRUE){
  infl = influence(object)
  coefs = coef(object)
  if (!is.null(random.seed)){
    set.seed(random.seed)
  }
  # if cond.data is not given, then use the covariates in the formula
  if (is.null(cond.data)){
    all.coef = names(object$coefficients)[2:length(object$coefficients)] 
    cond.data = object$model[all.coef]
  }
  
  # conditional inference for all coefficients by default
  if (is.null(param)){param = 1:length(object$coefficients)}
  names = param
  
  for (i.par in 1:length(param)){
    if (!is.na(suppressWarnings(as.integer(param[i.par])))){
      names[i.par] = names(object$coefficients)[as.integer(param[i.par])]
    }
  }
  
  n = length(object$fitted.values)
  infl.vals = matrix(n * infl$coefficients[,names], ncol=length(param))
  fitted.coef = coefs[names]
  hat.sigmas = rep(0, length(param))
  sup.sds = rep(0, length(param))
  ci.lows = rep(0, length(param))
  ci.his = rep(0, length(param))
  
  
  # sample splitting for cross fitting
  if (is.null(folds)){
    fold1 = sample(1:n, floor(n/2)) 
  }else{
    fold1 = folds[[1]]
  } 
  fold2 = setdiff(1:n, fold1)
  cond.data1 = data.frame(cond.data[fold1,])
  colnames(cond.data1) = colnames(cond.data)
  cond.data2 = data.frame(cond.data[fold2,])
  colnames(cond.data2) = colnames(cond.data)
  
  for (i.par in 1:length(param)){
    fit.Zr = infl.vals
    
    # LOESS regression
    if (alg == 'loess'){
      loess.span = ifelse(is.null(other.params$span), 0.75, other.params$span)
      loess.deg = ifelse(is.null(other.params$degree), 2, other.params$degree) 
      # cross-fit nonparametric regression models
      Zr.1 = loess(infl.vals[fold2,i.par]~., data=cond.data2,
                   span=loess.span, degree=loess.deg) 
      Zr.2 = loess(infl.vals[fold1,i.par]~., data=cond.data1,
                   span=loess.span, degree=loess.deg) 
      fit.Zr[fold1,i.par] = predict(Zr.1, cond.data1)
      fit.Zr[fold2,i.par] = predict(Zr.2, cond.data2)
    }
    
    # random forest regression
    if (alg == 'grf'){
      # cross-fit nonparametric regression models
      Zr.1 = regression_forest(cond.data2, infl.vals[fold2,i.par], num.threads=1)
      Zr.2 = regression_forest(cond.data1, infl.vals[fold1,i.par], num.threads=1)
      fit.Zr[fold1,i.par] = predict(Zr.1, cond.data1)$predictions
      fit.Zr[fold2,i.par] = predict(Zr.2, cond.data2)$predictions
    }
    
    hat.sigmas[i.par] = min(sd(infl.vals[,i.par]), 
                            sd(infl.vals[,i.par] - fit.Zr[,i.par], na.rm=TRUE))
    sup.sds[i.par] = sd(infl.vals[,i.par])
    ci.lows[i.par] = fitted.coef[i.par] - qnorm(0.975) * hat.sigmas[i.par] / sqrt(n)
    ci.his[i.par] = fitted.coef[i.par] + qnorm(0.975) * hat.sigmas[i.par] / sqrt(n)
  }
  
  # print the results 
  sup.p_vals = 2 * pnorm(abs(fitted.coef) / (sup.sds/sqrt(n)), lower.tail=FALSE)
  cond.p_vals = 2 * pnorm(abs(fitted.coef) / (hat.sigmas/sqrt(n)), lower.tail=FALSE)
  ret_table = cbind(fitted.coef, hat.sigmas, cond.p_vals, sup.sds, sup.p_vals)
  colnames(ret_table) = c("Estimate", "Cond. Std. Error", "Cond. Pr(>|z|)",
                          "Sup. Std. Error", "Sup. Pr(>|z|)")
  rownames(ret_table) = names
  
  if (verbose){
    cat("\n")
    cat("Summary of conditional inference")
    cat("\n\n")
    print(ret_table)
    cat("\n")
  }
  
  
  invisible(list("cond.std.err" = hat.sigmas, "std.err" = sup.sds,
              "fitted.coef" = fitted.coef, 
              "cond.ci.low" = ci.lows, "cond.ci.upp" = ci.his,
              "summary" = ret_table))
}