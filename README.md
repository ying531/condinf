# condinf
R package that implements conditional inference of parameters from lm and glm models. Details of conditional inference is in [this paper](https://arxiv.org/abs/2104.04565). The current implementation uses the influence function from `lm()` and `glm()` in R and two-fold cross-fitting for conditional inference standard error estimation. 

## Installation

1. The [devtools](https://github.com/hadley/devtools) package has to be installed. You can install it using `install.packages("devtools")`.
2. The latest development version can then be installied using `devtools::install_github("ying531/condinf")`.

## Helper

##### Usage


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

This function wraps around any `lm()` or `glm()` model. It prints the summary of conditional inference on the specified coefficients, and returns a list of results as described in the following. 

| Arguments      | Description                                                  |
| -------------- | ------------------------------------------------------------ |
| `model`        | An object returned from lm() or glm() functions              |
| `Z`            | A dataframe for the conditioning set                         |
| `param`        | A vector of coefficients to conduct conditional inference; can be a mixture of string name and index |
| `alg`          | Optinal, a string for name of algorithm, current options are 'loess' and 'grf' |
| `random.seed`  | Optional, random seed for sample splitting                   |
| `other.params` | Optional, other parameters for the regression algorithm; can include span and degree for loess |
| `folds`        | Optional, a list of two folds of indices for sample splitting; can be useful to control sample splitting |

To use `alg = 'grf'` as the regressor, the R package `grf` is required to be installed.

| Output         | Description                                                  |
| -------------- | ------------------------------------------------------------ |
| `cond.std.err` | Estimated standard error for inference of conditional parameters |
| `std.err`      | Estimated standard error for inference of super-population parameters |
| `fitted.coef`  | Fitted (empirical) coefficient from the model                |
| `cond.ci.low`  | Lower 0.95-confidence bound for conditional parameter        |
| `cond.ci.upp`  | Upper 0.95-confidence bound for conditional parameter        |



## Examples


### Conditional inference

#### Conditional inference for linear models

The following example works out conditional inference of linear regression coefficients (setting `param=1` selects the intercept) for a well-specified linear model. In this case, the super-population parameter is the same as conditional parameter, and the standard errors are the same. 

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
> Z = data.frame(X[,1:2])
> lm.mdl = lm(Y~., data = data.frame(X)) 
> cond.inf(lm.mdl, Z, param=1)

Summary of conditional inference

                Estimate Cond. Std. Error Cond. Pr(>|z|) Sup. Std. Error Sup. Pr(>|z|)
(Intercept) 0.0001085111        0.1081838      0.9746966       0.1081838     0.9746966
```

In the above summary, `Estimate` is the original estimator, `Cond. Std. Err` is the estimated standard error for inferring the conditional parameter, and `Cond. Pr(>|z|)` is the p-value for testing whether the conditional parameter is zero. `Sup. Std. Error` and `Sup. Pr(>|z|)` are the standard error and p-value for standard super-population inference. 



The following example conducts conditional inference for a misspecified linear model. In this case, the inference for the conditional parameter can be different from that for the super-population parameter (see the `cond.std.err` for conditional inference and `std.err` for super-population inference). The regression algorithm is `grf` and we focus on the coefficients for `"X1"` and `"X2"`.

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
> Z = data.frame(X[,1:2])
> lm.mdl = lm(Y~., data = data.frame(X))
> cond.inf(lm.mdl, Z, c("X1", "X2"), alg='grf')

Summary of conditional inference

    Estimate Cond. Std. Error Cond. Pr(>|z|) Sup. Std. Error Sup. Pr(>|z|)
X1 0.9618222         1.893261   4.478929e-58        3.589425  2.378448e-17
X2 2.0099977         1.137924   0.000000e+00        1.423715  0.000000e+00
```



#### Conditional inference for generalized linear models 

The following example works out conditional inference for parameters from a well-specified logistic model. The regression algorithm is `loess`. We focus on the coefficients for the second variable (`"X1"`).

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
> Y = rbinom(1000, 1, exp(logit.x)/(1+exp(logit.x)))
> Z = data.frame(X[,1:2])
> glm.mdl = glm(Y~., data = data.frame(X), family='binomial')
> cond.inf(glm.mdl, Z, 2)

Summary of conditional inference

  Estimate Cond. Std. Error Cond. Pr(>|z|) Sup. Std. Error Sup. Pr(>|z|)
2  1.10893         3.710156   3.332376e-21        3.710156  3.332376e-21
```

This is an example of conditional inference for parameters from a misspecified logistic regression model. The regression algorithm is `loess`. We focus on the coefficients for variable `"X3"` and the third variable (`"X2"`).

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
> Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x)))
> Z = data.frame(X[,1:2])
> glm.mdl = glm(Y~., data = data.frame(X), family='binomial')
> cond.inf(glm.mdl, Z, c("X3", 3))

Summary of conditional inference

   Estimate Cond. Std. Error Cond. Pr(>|z|) Sup. Std. Error Sup. Pr(>|z|)
X3 2.631113         4.492922   1.459031e-76        4.500887  2.680111e-76
X2 2.005930         4.116081   1.379739e-53        4.193597  1.088415e-51
```

