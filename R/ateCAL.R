# %% ####################################################
#' ateCAL: reweight experimental results to match external sample
#' @description Generalization or transportation estimator with method of
#'  moments weighting. Imputes counterfactual outcome for each observation i under each treament a as
#' \deqn{Y^a = \omega \frac{A = a}{\pi^a (X)} (Y - \hat{\mu}^a(X)) + \mu^a(X) }
#' Where \eqn{\omega} is solved for using entropy balancing for balance with target moments.
#' @param y outcome vector [length n]
#' @param a treatment vector [length n]
#' @param X covariate matrix [n X k]
#' @param targetMoments k-vector of target moments (means in overall sample for
#'  generalization, means in target sample for transportation
#' @param treatProb treatment probability, non-null results in no pscore fit
#' @param nuismod one of c("rlm", "rf") : choose how to fit nuisance functions
#'  (cross-fit) for  outcome models
#' @param estimator one of c("CW", "ACW"). Calibration weights (CW), which fit a set of balancing weights that reweights the sample to match target sample moments. ACW augments this with an outcome model (default)
#' @param glmnet_lamchoice choice of lambda (shrinkage parameter) for regularized linear regressions. Only relevant when nuismod == "rlm"
#' @param glmnet_alpha in [0, 1], choice of alpha in glmnet. 1 (default) corresponds with L1 regularization (LASSO) and 0 corresponds with L2 regularization (ridge), while intermediate values correspond with a mix of the two (elastic net)
#' @param separateMus boolean for whether to fit separate outcome models for each treatment group or a single pooled model. The former is recommended and is the default, but a pooled model may be fit when data is scarce / computation is burdensome.
#' @param noi boolean for printing marginal means and causal contrasts table (it gets returned anyway). On by default.
#' @return list containing treatment effects table , nuisance function estimates, and influence function values
#' @examples
#' df = selDGP(n = 10000)
#' id = with(df, s == 1)
#' X = df$X[id == 1, ]; y = df$y[id == 1]; a = df$a[id == 1]
#' Xtar = colMeans(df$X[!id, ])
#' ateCAL(y, a, X, targetMoments = Xtar, estimator = "ACW") %>% summary()
#' @export

ateCAL = function(y,
                  a,
                  X,
                  targetMoments,
                  treatProb = NULL,
                  nuismod = c("rlm", "rf"),
                  estimator = c("ACW", "CW"),
                  glmnet_lamchoice = "lambda.min",
                  glmnet_alpha = 1,
                  separateMus = TRUE,
                  noi = FALSE) {
  # housekeeping
  nuismod = match.arg(nuismod); estimator = match.arg(estimator)
  n = dim(X)[1]; avals = names(table(a)); n.avals = length(avals)

  if (!is.null(treatProb)) {
    if (noi) cat("Propensity scores passed manually. No pscore model will be fit. \n")
    # turn it into a vector if only p(a = 1) is passed
    if (length(treatProb) == 1) treatProb = c(treatProb, 1 - treatProb)
  }

  ############################################################
  # fit nuisance functions
  if (nuismod == "rlm") { # calls glmnet
    nuis = nuisRLM(y, a,
      s = NULL, X, n, avals, n.avals,
      estimator = estimator, treatProb = treatProb, SepMu = separateMus,
      lamChoice = glmnet_lamchoice, alph = glmnet_alpha, noi = noi
    )
  } else if (nuismod == "rf") {
    nuis = nuisGRF(y, a,
      s = NULL, X, n, avals, n.avals,
      estimator = estimator, treatProb = treatProb, SepMu = separateMus, noi = noi
    )
  } else {
    stop("Nuisance function model not supported.")
  }

  ############################################################
  # estimation - reweighting to make the data look like target
  ############################################################

  # solve for weights that set covariate means to target means in whole sample
  ebwts = eb_solve_dual(c(1, targetMoments), cbind(1, X))$Weights.ebal
  # need to scale each weight by n1 so that averaging works
  wtmat = matrix(rep(ebwts * n, n.avals), nrow = n, byrow = FALSE)
  if (estimator == "CW") {
    ipw_piece = wtmat * ((nuis$amat == nuis$alevel) / nuis$pihat) * nuis$ymat
    ifvals = ipw_piece
  } else if (estimator == "ACW") {
    balfit = wtmat * ((nuis$amat == nuis$alevel) / nuis$pihat) * (nuis$ymat - nuis$muhat)
    ifvals = as.matrix(balfit) + nuis$muhat
  } else {
    stop("Estimator not supported.")
  }
  # causal effects and marginal means
  res = estQois(ifvals, avals, noi = noi)
  outobj = list(
    res = res,
    nuis = nuis,
    ifvals = as.data.frame(ifvals),
    estimator = estimator,
    targetMoments = targetMoments,
    n = dim(X)[1],
    nA = length(avals)
  )
  class(outobj) = "atecal"
  return(invisible(outobj))
}

#' summary method - prints treatment effects and marginal means
#' @param fitted atecal object
#' @export
summary.atecal = function(obj) print(obj$res)
