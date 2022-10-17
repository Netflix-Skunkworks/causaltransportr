#' Simulate an experimental dataset with selection bias
#' @param n         sample size (default 10000)
#' @param p         number of covariates (default 10)
#' @param xb        bounds of uniform random variables for X (default -1, 1)
#' @param piF       propensity score function (fixed number indicates simple randomization - default 1/2)
#' @param tauF      treatment heterogeneity function (nonlinear by default)
#' @param y0F       baseline outcome function (nonlinear by default)
#' @param selF      selection bias function (nonlinear by default). When null, no missingness is introduced.
#' @return          list with outcome (with missing values corresponding with s = 0), treatment, covariates, selection, and true treatment effect
#' @export
selDGP = function(n = 10000, p = 10,
                  treatProb = 0.5, xb = c(-1, 1),
                  # propensity score function
                  piF = function(x) 1 / 2,
                  # nonlinear heterogeneity
                  tauF = function(x) 1 / exp(-x[3]),
                  # nonlinear y0
                  y0F = function(x) 3 * pmax(x[1] + x[2], 0) + 5 * sin(x[5]) * 2 * pmax(x[7], 0.5),
                  # nonlinear selection
                  selF = function(x) x[1] - 5 * x[3] + pmax(x[4], 0)) {
  X = matrix(runif(n * p, xb[1], xb[2]), n, p)
  # treatment selection
  pScore = apply(X, 1, tauF)
  a = rbinom(n, 1, plogis(pScore))
  # outcomes
  TAU = apply(X, 1, tauF)
  Y0 = apply(X, 1, y0F)
  # outcome switching
  y = (a * TAU + Y0 + rnorm(n))
  # selection - no selection bias if no selection function passed
  if (!is.null(selF)) {
    selscore = apply(X, 1, selF)
    s = rbinom(n, 1, plogis(selscore)) |> as.logical()
    # set outcomes for s = 0 as missing
    y[s == 0] = NA
  } else {
    s = NULL
  }

  list(
    y = y, a = a, X = X,
    s = s, tau = TAU
  )
}
