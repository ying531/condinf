# condinf
R package that implements conditional inference of parameters from lm and glm models. Details of conditional inference is in [this paper](https://arxiv.org/abs/2104.04565). The current implementation uses the influence function from `lm()` and `glm()` in R and two-fold cross-fitting for conditional inference standard error estimation. 

## Installation

1. The [devtools](https://github.com/hadley/devtools) package has to be installed. You can install it using `install.packages("devtools")`.
2. The latest development version can then be installied using `devtools::install_github("ying531/condinf")`.



### Helper

| Arguments      | Description                                                  |
| -------------- | ------------------------------------------------------------ |
| `model`        | An object returned from lm() or glm() functions              |
| `Z`            | A dataframe for the conditioning set                         |
| `alg`          | Optinal, a string for name of algorithm, current options are 'loess' and 'grf' |
| `random.seed`  | Optional, random seed for sample splitting                   |
| `other.params` | Optional, other parameters for the regression algorithm; can include span and degree for loess |
| `folds`        | Optional, a list of two folds of indices for sample splitting; can be useful to control sample splitting |

To use `alg = 'grf'` as the regressor, the R package `grf` is required to be installed.

## Usage

```
cond.inf(
  model,
  Z,
  param,
  alg = "loess",
  random.seed = NULL,
  other.params = NULL,
  folds = NULL
)
```

### Conditional inference

#### Conditional inference for linear models

The following example works out conditional inference of linear regression coefficients for a well-specified linear model. In this case, the super-population parameter is the same as conditional parameter, and the standard errors are the same. 

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
> Z = data.frame(X[,1:2])
> lm.mdl = lm(Y~., data = data.frame(X)) 
> cond.inf(lm.mdl, Z, param=1)
$cond.std.err
[1] 0.1016477

$std.err
[1] 0.1016477

$fitted.coef
 (Intercept) 
0.0002876785 

$cond.ci.low
 (Intercept) 
-0.006012397 

$cond.ci.upp
(Intercept) 
0.006587754 
```



The following example conducts conditional inference for a misspecified linear model. In this case, the inference for conditional parameter can be different from that for super-population parameter (see the `cond.std.err` for conditional inference and `std.err` for super-population inference). The regression algorithis `grf`.

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
> Z = data.frame(X[,1:2])
> lm.mdl = lm(Y~., data = data.frame(X))
> cond.inf(lm.mdl, Z, 2, alg='grf')
$cond.std.err
[1] 1.488649

$std.err
[1] 3.287264

$fitted.coef
      X1 
1.228476 

$cond.ci.low
     X1 
1.13621 

$cond.ci.upp
      X1 
1.320741 

```



#### Conditional inference for generalized linear models 

The following example works out conditional inference for parameters from a well-specified logistic model. The regression algorithm is `loess`. 

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
> Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x)))
> Z = data.frame(X[,1:2])
> glm.mdl = glm(Y~., data = data.frame(X), family='binomial')
> cond.inf(glm.mdl, Z, 2)
$cond.std.err
[1] 3.911343

$std.err
[1] 3.935879

$fitted.coef
      X1 
1.338414 

$cond.ci.low
      X1 
1.095991 

$cond.ci.upp
      X1 
1.580837 
```

This is an example of conditional inference for parameters from a misspecified logistic regression model. The regression algorithm is `loess`. 

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
> Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x)))
> Z = data.frame(X[,1:2])
> glm.mdl = glm(Y~., data = data.frame(X), family='binomial')
> cond.inf(glm.mdl, Z, 2)
$cond.std.err
[1] 3.184585

$std.err
[1] 3.492766

$fitted.coef
       X1 
0.7752392 

$cond.ci.low
       X1 
0.5778602 

$cond.ci.upp
       X1 
0.9726182 
```

