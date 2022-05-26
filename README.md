# condinf
R package that implements conditional inference of parameters from lm and glm models. Details of conditional inference is in [this paper](https://arxiv.org/abs/2104.04565). The current implementation uses the influence function from `lm()` and `glm()` in R and two-fold cross-fitting for conditional inference standard error estimation. 

## Installation

1. The [devtools](https://github.com/hadley/devtools) package has to be installed. You can install it using `install.packages("devtools")`.
2. The latest development version can then be installied using `devtools::install_github("ying531/condinf")`.

## Usage

This package computes the standard error and confidence intervals for conditional inference

### Conditional inference

The following example works out conditional inference of linear regression coefficients for a well-specified linear model. 

```R
X = matrix(rnorm(1000*10), nrow=1000)
Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
Z = data.frame(X[,1:2])
lm.mdl = lm(Y~., data = data.frame(X))
cond.inf(lm.mdl, Z, 2)
```

