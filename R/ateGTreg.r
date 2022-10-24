# %% ####################################################
#' Generalization AISW using fixest fixed effects regressions
#' @description Augmented IPW, generalization, or transport estimator with fast fixest regressions and computation in `data.table` (1.14.5), no regularization, and no sample splitting. Imputes counterfactual outcome for each observation i under each treament a as
#' \eqn{Y^a = \frac{S}{\rho(X)} \frac{A = a}{\pi^a (X)} (Y - \mu^a(X)) + \mu^a(X) }
#' Where the first term is 1 for all observations under no sample selection, and therefore this is the doubly-robust Augmented Inverse Propensity Weighting (AIPW) estimator. When S is supplied, the argument in 'target' is used to fit either the generalization or transportation estimator.
#' Recommended with the nonparametric/bayesian bootstrap for inference.
#' @param d data.table
#' @param yn outcome name
#' @param an treatment indicator name
#' @param sn selection indicator name (null by default - this fits the AIPW regression)
#' @param xn list of covariates (default intercept)
#' @param fe list of fixed effects (default 0)
#' @param target estimand (generalization / transportation/ insample)
#' @return generalization effect and influence functions
#' @export

ateGTreg = function(d,
                    yn = "y",
                    an = "a",
                    sn = NULL,
                    xn = "1",
                    fe = "0",
                    target = c("generalize", "transport", "insample")) {
  # internal fn for LPM
  clip01 = function(x, lower = 0, upper = 1) pmax(pmin(x, upper), lower)
  # data prep
  target = match.arg(target)
  data = copy(d)
  # dummies for treatment
  avals = names(table(data[[an]]))
  data[, paste0("a_", 1:length(avals)) := lapply(
    avals,
    function(x) 1 * (get(an) == x)
  )]
  # selection model business
  if (is.null(sn)) {
    fitrho = FALSE
    data[, s := 1]
    onames = c(yn, an);      nnames = c("y", "a")
  } else {
    fitrho = TRUE
    onames = c(yn, an, sn);  nnames = c('y', 'a', 's')
  }
  # rename for ease later
  setnames(data, onames, nnames)
  # container df for all future computations
  df2 = data.table(data[, .(y, a, s)])
  # need dummies for counterfactual estimation too
  df2[, paste0("a_", 1:length(avals)) := lapply(
    avals,
    function(x) 1 * (a == x)
  )]
  ######################################################################
  # selection model
  ######################################################################
  if (fitrho) {
    selmod = feols(s ~ .[xn] | .[fe], data = data)
    df2$rho = clip01(predict(selmod, data))
  } else {
    df2$rho = 1
  }
  ######################################################################
  # treatment pscore
  ######################################################################
  # j - 1 propensity models
  for (j in 1:(length(avals) - 1)) {
    aname = paste0("a_", j)
    psmod = feols(.[aname] ~ .[xn] | .[fe], data = data[s == 1])
    df2[[paste0("pihat_", j)]] = clip01(predict(psmod, data))
  }
  # 1 - sum of others for last pscore
  pscols = grep("^pihat_*", colnames(df2), value = T)
  df2[[paste0("pihat_", length(avals))]] = df2[, 1 - rowSums(.SD),
    .SDcols = pscols
  ]
  ######################################################################
  # outcome model
  ######################################################################
  for (j in 1:length(avals)) {
    t = avals[j]
    # outcome model fit
    mumod = feols(y ~ .[xn] | .[fe], data = data[a == t & s == 1])
    df2[[paste0("muhat_", j)]] = predict(mumod, data)
    # ipw piece
    ipwlist = list(
      ipwpiece = paste0("ipw_", j),
      s = "s",
      rho = "rho",
      treatdum = paste0("a_", j),
      pscore = paste0("pihat_", j),
      y = "y",
      yhat = paste0("muhat_", j)
    )
    if (target == "generalize") {
      df2[, ipwpiece := (s / rho) * (treatdum) / (pscore) * (y - yhat), env = ipwlist]
    } else if (target == "transport") {
      df2[, ipwpiece := (s * (1 - rho) / rho) * (treatdum) / (pscore) * (y - yhat), env = ipwlist]
    } else if (target == "insample") {
      df2[, ipwpiece := (treatdum) / (pscore) * (y - yhat), env = ipwlist]
    }
    # this is missing for s=0
    df2[is.na(ipw_piece), ipw_piece := 0,
      env = list(ipw_piece = paste0("ipw_", j))
    ]
    # finally impute counterfactual
    df2[, yhat := ipw_piece + omod_piece,
      env = list(
        yhat       = paste0("yhat_", j),
        ipw_piece  = paste0("ipw_", j),
        omod_piece = paste0("muhat_", j)
      )
    ]
  }
  if (df2[rho == 0, .N] > 0) warning("Selection scores = 0 exist and will be dropped")
  df2 = df2[!(rho == 0)]
  # compute ATE for 2 treatment case
  if (length(avals) == 2) {
    df2[, esti := y1 - y0,
      env = list(esti = "gamma_i", y1 = "muhat_2", y0 = "muhat_1")
    ]
    ate_est = mean(df2[, gamma_i])
    return(list(
      if_nuis = df2,
      ate = ate_est
    ))
  } else {
    return(list(if_nuis = df2))
  }
}
