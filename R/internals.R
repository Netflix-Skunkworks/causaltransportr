# %%
#' prepare X matrix with flexible interactions and polynomials
#' @description Prepare matrix with polynomial basis expansion and/or interactions.
#' @param dat              data table / dataframe
#' @param dummies          Names of dummy vars (null by default - auto-detected)
#' @param continuouses     Names of continuous vars to modify fn form (null by default - auto detected)
#' @param corr_cut         cutoff for correlation threshold to drop one of the vars (default = 0.9)
#' @param k                order of interactions: defaults to pairwise interactions
#' @param m                order of functions: defaults to quadratic functions
#' @param raw              raw or orthogonal bases
#' @return                 data.frame with base terms and interactions + basis
#' @export

polySieveM = function(dat,
                      dummies = NULL,
                      continuouses = NULL,
                      corr_cut = 0.90,
                      k = 2,
                      m = 2,
                      raw = FALSE) {
  # coerce to df
  dat = as.data.frame(dat)
  nuniqs = apply(dat, 2, function(x) length(unique(x)))
  # populate lists if missing
  if (is.null(dummies)) dummies = colnames(dat)[which(nuniqs == 2)]
  if (is.null(continuouses)) continuouses = colnames(dat)[which(nuniqs >= 5)]
  # concat names
  controls = c(dummies, continuouses)
  data = dat[, controls]
  if (!is.null(continuouses)) { # if there are any continous variables
    # functional form changes for continuous vars
    # polynomials (orthogonal by default)
    polyfml = paste(paste0(
      "poly(", continuouses, ",", m,
      ", raw=", raw, ")"
    ), collapse = " + ")
    # model matrix with sieve polynomial basis
    powmat = model.matrix(as.formula(paste0("~ -1 + ", polyfml)), dat[, continuouses])
    # log variables with all positives
    loggable = apply(dat[, continuouses], 2, function(x) sum(x < 0) == 0)
    if (sum(loggable) > 0) {
      logmat = log1p(dat[, loggable]) %>% as.matrix()
      colnames(logmat) = paste0("log1p_", names(which(loggable == TRUE)))
      powmat = cbind(powmat, logmat)
    }
    # cbind base terms with smooth fns
    data = cbind(data, powmat)
  }
  if (k > 1) {
    # all n-way interactions with polynomials and dummies
    X = model.matrix(as.formula(glue("~.^{k} - 1")), data)
  } else {
    X = as.matrix(data)
  }
  # drop non-varying Xs (e.g. interactions bw mutually exclusive dummies)
  varyvar = apply(X, 2, function(col) length(unique(col)) > 1)
  X = X[, varyvar]
  # drop highly correlated Xs - this drops polynomials of binary variables
  corm = cor(X)
  hc = findCorr(corm, cutoff = corr_cut) # put any value as a "cutoff"
  hc = sort(hc)
  return(X[, -c(hc)])
}

#' prepare sparse X matrix with interactions
#' @description Prepare sparse matrix with  interactions.
#' @param data              data table / dataframe
#' @param k                order of interactions: defaults to pairwise interactions
#' @param corr_cut         cutoff for correlation threshold to drop one of the vars (default = 0.9)
#' @return                matrix with interactions
#' @export

interSparseM = function(data, k, corr_cut = 0.9) {
  X = Matrix::sparse.model.matrix(
    as.formula(glue("~.^{k} - 1")),
    data
  )
  # drop non-varying Xs (e.g. interactions bw mutually exclusive dummies)
  varyvar = apply(X, 2, function(col) length(unique(col)) > 1)
  X = X[, varyvar]
  # drop highly correlated Xs - this drops polynomials of binary variables
  corm = cor(X)
  hc = findCorr(corm, cutoff = corr_cut) # put any value as a "cutoff"
  hc = sort(hc)
  X[, -c(hc)]
}


######################################################################
# internals
######################################################################
#' min-max scale (maps continuous variable to [0, 1])
#' @param X vector
#' @export
mMscale = function(X) {
  X = as.matrix(X)
  mins = apply(X, 2, min)
  maxs = apply(X, 2, max)
  return(scale(X, center = mins, scale = maxs - mins))
}

# %% fast correlation matrix filtration
findCorr = function(x, cutoff = .90) {
  if (any(!complete.cases(x)))
    stop("The correlation matrix has some missing values.")
  averageCorr = colMeans(abs(x))
  averageCorr = as.numeric(as.factor(averageCorr))
  x[lower.tri(x, diag = TRUE)] = NA
  combsAboveCutoff = which(abs(x) > cutoff)
  colsToCheck = ceiling(combsAboveCutoff / nrow(x))
  rowsToCheck = combsAboveCutoff %% nrow(x)
  colsToDiscard = averageCorr[colsToCheck] > averageCorr[rowsToCheck]
  rowsToDiscard = !colsToDiscard
  deletecol = c(colsToCheck[colsToDiscard], rowsToCheck[rowsToDiscard])
  deletecol = unique(deletecol)
  deletecol
}

# %% glmnet
# internal to get cross-fit predictions for binomial or gaussian fits
fitGet = function(m, lam = "lambda.min") {
  # column index of best fitting model
  predcol = (m$lambda == m[[lam]])
  m$fit.preval[, predcol]
}

elapExtract = function(t) unname(round(t[3]))


clip = function(x, lower = 0, upper = 1) pmax(pmin(x, upper), lower)
