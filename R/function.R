suppressMessages(require("grf"))

#' Conditional inference  
#'
#' Conditional inference for lm and glm models
#' @param model An object returned from lm() or glm() functions
#' @param Z A dataframe for the conditioning set
#' @param alg Optinal, a string for name of algorithm, current options are 'loess' and 'grf'
#' @param random.seed Optional, random seed for sample splitting
#' @param other.params Optional, other parameters for the regression algorithm; can include span and degree for loess
#' @param folds Optional, a list of two folds of indices for sample splitting; can be useful to control sample splitting
#' @return Standard error for conditional parameter, super-population parameter, fitted empirical parameter, confidence interval for conditional parameter
#' @examples 
#' X = matrix(rnorm(1000*10), nrow=1000)
#' Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
#' Z = data.frame(X[,1:2])
#' lm.mdl = lm(Y~., data = data.frame(X))
#' cond.inf(lm.mdl, Z, 2) 
#' 
#' X = matrix(rnorm(1000*10), nrow=1000)
#' logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
#' Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x)))
#' Z = data.frame(X[,1:2])
#' glm.mdl = glm(Y~., data = data.frame(X), family='binomial')
#' cond.inf(glm.mdl, Z, 2)
#' 
#' @export
cond.inf <- function(model,Z,param,alg="loess",
                     random.seed=NULL,other.params=NULL,folds=NULL){
  infl = influence(model)
  coefs = coef(model)
  if (!is.null(random.seed)){
    set.seed(random.seed)
  }
  
  n = length(model$fitted.values)
  infl.vals = n * infl$coefficients[,param]
  fitted.coef = coefs[param]
  
  # sample splitting for cross fitting
  fold1 = ifelse(is.null(folds[[1]]), sample(1:n, floor(n/2)), folds[[1]]) 
  fold2 = setdiff(1:n, fold1)
  fit.Zr = rep(0, n)
  
  # LOESS regression
  if (alg == 'loess'){
    loess.span = ifelse(is.null(other.params$span), 0.75, other.params$span)
    loess.deg = ifelse(is.null(other.params$degree), 2, other.params$degree) 
    # cross-fit nonparametric regression models
    Zr.1 = loess(infl.vals[fold2]~., data=data.frame(Z[fold2,]),
                 span=loess.span, degree=loess.deg) 
    Zr.2 = loess(infl.vals[fold1]~., data=data.frame(Z[fold1,]),
                 span=loess.span, degree=loess.deg) 
    fit.Zr[fold1] = predict(Zr.1, data.frame(Z[fold1,]))
    fit.Zr[fold2] = predict(Zr.2, data.frame(Z[fold2,]))
  }
  
  # random forest regression
  if (alg == 'grf'){
    # cross-fit nonparametric regression models
    Zr.1 = regression_forest(data.frame(Z[fold2,]), infl.vals[fold2], num.threads=1)
    Zr.2 = regression_forest(data.frame(Z[fold1,]), infl.vals[fold1], num.threads=1)
    fit.Zr[fold1] = predict(Zr.1, data.frame(Z[fold1,]))$predictions
    fit.Zr[fold2] = predict(Zr.2, data.frame(Z[fold2,]))$predictions
  }
  
  hat.sigma = min(sd(infl.vals), sd(infl.vals - fit.Zr, na.rm=TRUE))
  ci.low = fitted.coef - qnorm(0.975) * hat.sigma / sqrt(n)
  ci.hi = fitted.coef + qnorm(0.975) * hat.sigma / sqrt(n)
  
  return(list("cond.std.err" = hat.sigma, "std.err" = sd(infl.vals),
              "fitted.coef" = fitted.coef, 
              "cond.ci.low" = ci.low, "cond.ci.upp" = ci.hi))
}

# trans.inf <- function(model,Z,new.Z,alg="loess"){
#   
# }



X = matrix(rnorm(1000*10), nrow=1000)
Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
Z = data.frame(X[,1:2])
lm.mdl = lm(Y~., data = data.frame(X))
cond.inf(lm.mdl, Z, 2)

X = matrix(rnorm(1000*10), nrow=1000)
Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,4]**2 + rnorm(1000) * 0.1
Z = data.frame(X[,1:2])
lm.mdl = lm(Y~., data = data.frame(X))
cond.inf(lm.mdl, Z, 2, alg='grf')


X = matrix(rnorm(1000*10), nrow=1000)
logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x)))
Z = data.frame(X[,1:2])
glm.mdl = glm(Y~., data = data.frame(X), family='binomial')
cond.inf(glm.mdl, Z, 2)


X = matrix(rnorm(1000*10), nrow=1000)
logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x)))
Z = data.frame(X[,1:2])
glm.mdl = glm(Y~., data = data.frame(X), family='binomial')
cond.inf(glm.mdl, Z, 2)
