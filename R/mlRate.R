# %% ####################################################
#' Flexible covariate adjustment using ML predictions
#' @description
#' First fits a predictive model of covariates on outcome \eqn{Y} to construct fitted values \eqn{g(X)}.
#' Then, it runs ML-based version of the Lin (2013) regression with predictions \eqn{g(X)}.
#' \eqn{Y = \beta_0 + \beta_1 a + \beta_2 g(X) + \beta_4 (a \times (g(X) - \bar{g}(X) )}
#' where the treatment is interacted with predictions g(). This a dimension reduced veresion of \eqn{X}s instead of using the full matrix.
#' Treatment effect is \eqn{\beta_1} with a heteroskedasticity robust variance estimate that is asymptotically bounded from above by the naive Neyman variance but typically lower if covariates explain substantial variation in Y.
#' @param y outcome vector
#' @param a treatment dummy vector
#' @param X covariate matrix
#' @param nuisMod ML algorithm to fit nuisance function. Defaults to glmnet, can also use generalized random forests.
#' @param glmnet_lamchoice choice of lambda (shrinkage parameter) for regularized linear regressions. Only relevant when nuisMod == "rlm"
#' @param glmnet_alpha in [0, 1], choice of alpha in glmnet. 1 (default) corresponds with L1 regularization (LASSO) and 0 corresponds with L2 regularization (ridge), while intermediate values correspond with a mix of the two (elastic net)
#' @param glmnet_mu_family  GLM family for outcome model. Gaussian by default.
#' @param glmnet_parl Boolean for parallelization in glmnet. Need to enable parallelized cluster beforehand.
#' @param tuneRf boolean for whether to tune RF
#' @param noi Boolean for noisy
#' @return treatment effect and SE
#' @export
#' @references Guo, Y., Coey, D., Konutgan, M., Li, W., Schoener, C., & Goldman, M. (2021). Machine learning for variance reduction in online experiments. Advances in Neural Information Processing Systems, 34, 8637-8648.
#' @references Lin, Winston. "Agnostic notes on regression adjustments to experimental data: Reexamining Freedman’s critique." The Annals of Applied Statistics 7.1 (2013): 295-318.
mlRate = function(y, a, X,
                  nuisMod = c("rlm", "rf"), # model to fit double robust score - can add more
                  glmnet_lamchoice = "lambda.min",
                  glmnet_alpha = 1,
                  glmnet_mu_family = "gaussian",
                  glmnet_parl = FALSE,
                  # rf choices
                  tuneRf = "none",
                  noi = FALSE) {
  # housekeeping
  nuisMod = match.arg(nuisMod)
  n = dim(X)[1]; avals = names(table(a))
  n.avals = length(avals)
  ############################################################
  # fit nuisance functions
  ############################################################
  if (nuisMod == "rlm") { # calls glmnet
    gFit = cv.glmnet(X, y,
      standardize = FALSE,
      keep = TRUE,
      family = glmnet_mu_family,
      alpha = glmnet_alpha,
      trace.it = noi,
      parallel = glmnet_parl
    )
    g = fitGet(gFit, glmnet_lamchoice)
  } else if (nuisMod == "rf") {
    gFit = regression_forest(X, y, tune.parameters = tuneRf)
    g = predict(gFit)[, 1]
  }
  X = cbind(1, a, g, a * (g - mean(g)))
  # fit regression
  beta = lm.fit(X, y)$coefficients
  # variance
  phat = mean(a)
  sighat = (
    # vanilla terms
    var(y[a == 1]) / (phat) + var(y[a == 0]) / (1 - phat) -
      # reductions from ml
      var(g) / (phat * (1 - phat)) *
        (beta[3] * phat + (beta[3] + beta[4]) * (1 - phat))
  )
  c(beta[2], sqrt(sighat / n))
}
