# %% ####################################################
#' AIPW estimation under unconfoundedness
#' @description Simplified version of `ateGT` with no selection indicator. Imputes counterfactual outcome for each observation i under each treament a as
#' \deqn{Y^a = \frac{A = a}{\pi^a (X)} (Y - \mu^a(X)) + \mu^a(X) }
#' Average treatment effects are defined as averages of differences between counterfactual outcomes \eqn{Y^a - Y^{a'}}.
#' and returns marginal means, causal contrasts, and their influence-function based variance estimates. Uses same underlying nuisance function fitter as `ateGT` with no missing data or surrogates.
#' @param y outcome vector (may contain missings ; missings must correspond with s = 0)
#' @param a treatment vector (no missings; can be relaxed with some tinkering)
#' @param X covariate matrix (no missings)
#' @param treatProb propensity score vector (of length n_treatment) or matrix (n_treatment X n_obs), where latter is for covariate adaptive designs; must sum to 1. NULL by default, so pscore is fitted. When provided, no propensity score is fit. With discrete covariates, estimated propensity score is advisable even if treatment was randomized.
#' @param nuisMod one of c("rlm", "rf") : choose how to fit nuisance functions (cross-fit).
#' @param estimator one of c("AIPW", "IPW", "OM"). The default is AIPW.
#' @param separateMus boolean for whether to fit separate outcome models for each treatment group or a single pooled model. The former is recommended and is the default, but a pooled model may be fit when data is scarce / computation is burdensome.
#' @param glmnet_lamchoice choice of lambda (shrinkage parameter) for regularized linear regressions. Only relevant when nuisMod == "rlm"
#' @param glmnet_alpha in [0, 1], choice of alpha in glmnet. 1 (default) corresponds with L1 regularization (LASSO) and 0 corresponds with L2 regularization (ridge), while intermediate values correspond with a mix of the two (elastic net)
#' @param glmnet_rho_family  GLM family for selection model. "binomial" by default but can be safely switched to "gaussian" for linear probability models with discrete covariates for faster compute
#' @param glmnet_pi_family  GLM family for propensity model. "binomial" by default but can be safely switched to "gaussian" for linear probability models with discrete covariates for faster compute
#' @param glmnet_mu_family  GLM family for outcome model. Gaussian by default.
#' @param glmnet_parl Boolean for parallelization in glmnet. Need to enable parallelized cluster beforehand.
#' @param hajekize boolean for whether to divide the inverse probability weights term for each treatment level by the sum of weights in that treatment level. This guards against instability from very large weights from extremely small selection or propensity scores.
#' @param grf_tunerf Tune rf hyperparameters? Passed to grf's regression forest. Use 'all' for hyperparameter tuning.
#' @param noi boolean for printing marginal means and causal contrasts table (it gets returned anyway). off by default.
#' @return list containing treatment effects table , nuisance function estimates, and influence function values
#' @export
#' @examples
#' # simulation with no selection bias (generate by passing null function for selection)
#' df2 = selDGP(selF = NULL)
#' df2$tau |> mean() # true effect : sinh(1)
#' # lasso
#' aipw(df2$y, df2$a, df2$X, noi = FALSE, nuisMod = "rlm")$res
#' # random forest
#' aipw(df2$y, df2$a, df2$X, noi = FALSE, nuisMod = "rf")$res
aipw = function(y,
                a,
                X,
                treatProb = NULL, # treatment probability - if supplied no propensity score fit
                nuisMod = c("rlm", "rf"), # model to fit double robust score - can add more
                estimator = c("AIPW", "IPW", "OM"),
                hajekize = FALSE,
                separateMus = TRUE,
                # glmnet choices
                glmnet_lamchoice = "lambda.min",
                glmnet_alpha = 1,
                glmnet_rho_family = "binomial",
                glmnet_pi_family = "binomial",
                glmnet_mu_family = "gaussian",
                glmnet_parl = FALSE,
                # rf choices
                grf_tunerf = "none",
                # misc
                noi = FALSE) {
  # housekeeping
  nuisMod = match.arg(nuisMod)
  estimator = match.arg(estimator)

  n = dim(X)[1]
  avals = names(table(a)); n.avals = length(avals)

  # nuisance function fitting - pass ACW to avoid fitting selection score
  if (nuisMod == "rlm") { # calls glmnet
    nuis = nuisRLM(
      y, a, X, n, avals, n.avals,
      estimator = "ACW",
      treatProb = treatProb,
      s = NULL, Z = NULL, # passing nulls means no selection score fit
      SepMu = separateMus,
      # glmnet parameters
      lamChoice = glmnet_lamchoice,
      alph = glmnet_alpha,
      rhoF = glmnet_rho_family,
      piF = glmnet_pi_family,
      muF = glmnet_mu_family,
      parl = glmnet_parl,
      # others
      noi = noi
    )
  } else if (nuisMod == "rf") {
    nuis = nuisGRF(
      y, a, X, n, avals, n.avals,
      s = NULL, Z = NULL, # passing nulls means no selection score fit
      estimator = "ACW",
      treatProb = treatProb,
      SepMu = separateMus,
      noi = noi,
      tune = grf_tunerf
    )
  } else {
    stop("Nuisance function model not supported.")
  }
  if (sum(is.na(y)) > 0) stop("outcome contains missings; fit with sample correction")

  # compute influence function
  if (estimator == "AIPW") {
    # classic infunc
    ipw_piece = ((nuis$amat == nuis$alevel) * (nuis$ymat - nuis$muhat)) / (nuis$pihat)
    if (hajekize) ipw_piece = ipw_piece / colSums(nuis$pihat)
    # add in outcome model piece
    ifvals = ipw_piece + nuis$muhat
  } else if (estimator == "ISW") {
    ipw_piece = ((nuis$amat == nuis$alevel) * nuis$ymat) / (nuis$pihat)
    if (hajekize) ipw_piece = ipw_piece / colSums(nuis$pihat)
    ifvals = ipw_piece
  } else if (estimator == "OM") {
    ifvals = nuis$muhat
  }

  # estimate marginal means and treatment effects and their SEs
  res = estQois(ifvals, avals, noi = noi)

  outobj = list(
    res       = res,
    nuis      = nuis,
    ifvals    = as.data.frame(ifvals),
    nuisMod   = nuisMod,
    estimator = estimator,
    n         = n,
    nA        = n.avals
  )
  class(outobj) = "aipw"

  return(invisible(outobj))
}

#' Summary method - prints treatment effects and marginal means
#' @param fitted aipw object
#' @export
summary.aipw = function(obj) print(obj$res)
