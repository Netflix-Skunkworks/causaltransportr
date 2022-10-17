# %% LASSO
nuisRLM = function(y,
                   a,
                   s,
                   X,
                   n,
                   avals,
                   n.avals,
                   estimator = "AISW",
                   treatProb = NULL,
                   Z = NULL,
                   SepMu = TRUE,
                   lamChoice = "lambda.min",
                   alph = 1,
                   rhoF = "binomial",
                   piF = "binomial",
                   muF = "gaussian",
                   noi = TRUE,
                   parl = FALSE) {
  # for standard matrices - put on [0,1] scale; do this manually for sparse matrices
  if (inherits(X, "matrix")) X = mMscale(X)
  # progress bar
  tracer = 0
  if (noi) tracer = 1
  # sample splits for cross-fitting
  # ROT for nfolds
  k_folds = floor(max(3, min(10, length(a) / 4)))
  foldid = sample(rep(seq(k_folds), length = length(a)))
  if (parl) cat("Parallel requires registered parallel backend. Use doParallel::registerDoParallel() before call. \n")
  # switches - fit only necessary nuisance functions
  if (estimator == "AISW") {
    fitmu = fitrho = fitpi = TRUE
  } else if (estimator == "ISW") {
    fitrho = fitpi = TRUE
    fitmu = FALSE
  } else if (estimator == "OM") {
    fitrho = fitpi = FALSE
    fitmu = TRUE
  } else if (estimator == "CW") {
    fitrho = fitmu = FALSE
    fitpi = TRUE
  } else if (estimator == "ACW") {
    fitrho = FALSE
    fitmu = fitpi = TRUE
  }
  if (fitpi == FALSE) treatProb = 0
  # post-hoc mods conditional on inputs
  if (is.null(s)) {
    fitrho = FALSE; s = rep(1, n)
  }
  if (!is.null(treatProb)) fitpi = FALSE
  if (!is.null(Z)) {
    surrogates = TRUE
  } else {
    surrogates = FALSE
  }
  # containers
  muhat = as.data.frame(matrix(NA, nrow = n, ncol = n.avals))
  colnames(muhat) = paste("a", avals, sep = "")
  pihat = nuhat = muhat
  ######################################################################
  # selection model
  if (fitrho) {
    rhohat = array(NA, dim = n)
    if (surrogates) { # include surrogates in selection score
      if (noi) cat("fitting surrogate selection model\n")
      time1 = system.time({
        rhoMod = cv.glmnet(cbind(X, Z), s,
          standardize = FALSE,
          keep = TRUE,
          foldid = foldid,
          family = rhoF,
          alpha = alph,
          trace.it = tracer,
          parallel = parl
        )
      })
      if (noi) cat(glue("Selection model fitting took {elapExtract(time1)} seconds"), "\n")
    } else {
      if (noi) cat("fitting selection model\n")
      time1 = system.time({
        rhoMod = cv.glmnet(X, s,
          standardize = FALSE,
          keep = TRUE,
          family = rhoF,
          foldid = foldid,
          alpha = alph,
          trace.it = tracer,
          parallel = parl
        )
      })
      if (noi) cat(glue("Selection model fitting took {elapExtract(time1)} seconds"), "\n")
    }
    # selection score for all obs
    rhohat = fitGet(rhoMod, lamChoice)
    if (rhoF == "binomial") rhohat = plogis(rhohat)
  } else { # default selection score is 1 (for OM/ACW etc)
    rhohat = rep(1, n)
  }
  ######################################################################
  # propensity model
  if (fitpi) {
    # first k-1 treatments
    for (j in 1:(n.avals - 1)) {
      ids = s == 1
      # separate model for E[A  = a | X = x, S = 1] for each a
      if (noi) cat(glue("fitting propensity for treatment {avals[j]}\n"))
      time3 = system.time({
        pimod = cv.glmnet(
          X[ids, ],
          as.numeric(a == avals[j])[ids],
          keep = TRUE,
          standardize = FALSE,
          foldid = foldid[ids],
          family = piF,
          parallel = parl,
          trace.it = tracer,
          alpha = alph
        )
      })
      if (noi) cat(glue("Propensity model {j} fitting took {elapExtract(time3)} seconds"), "\n")
      pihat[, j] = predict(pimod, X, s = lamChoice, type = 'response')
      # oob for obs used to train model
      pihat[ids, j] = fitGet(pimod, lamChoice)
      if (piF == "binomial") pihat[ids, j] = plogis(fitGet(pimod, lamChoice))
    }
    # last pscore is 1- sum(others)
    pihat[, n.avals] = 1 - rowSums2(as.matrix(pihat), na.rm = T)
  } else {
    if (inherits(treatProb, "matrix")) { # pscore is a matrix (already n X k)
      pihat = treatProb
    } else { # repeat each pscore vector k times
      pihat = matrix(rep(treatProb, n), nrow = n, byrow = FALSE)
    }
  }
  ######################################################################
  # outcome model
  if (fitmu) {
    for (j in 1:n.avals) {
      ids = a == avals[j] & s == 1
      if (noi) cat(glue("fitting outcome model for treatment {avals[j]}"), "\n")
      # separate model for E[Y | A = a, X = x, S = 1] for each a
      time3 = system.time({
        muMod = cv.glmnet(X[ids, ], y[ids],
          keep = TRUE,
          standardize = FALSE,
          family = muF,
          foldid = foldid[ids],
          alpha = alph,
          trace.it = tracer,
          parallel = parl
        )
      })
      if (noi) cat(glue("Outcome model {j} fitting took {elapExtract(time3)} seconds"), "\n")
      muhat[, j] = predict(muMod, X, s = lamChoice)
      # oob for obs used to train model
      muhat[ids, j] = fitGet(muMod, lamChoice)
    }
    # additional model  for surrogates
    if (surrogates) { # additional outcome model with surrogates
      for (j in 1:n.avals) {
        ids = a == avals[j] & s == 1
        if (noi) cat(glue("fitting surrogate outcome model for treatment {avals[j]}"), "\n")
        # separate model for E[Y | A = a, X = x, Z = z, S = 1] for each a
        time3 = system.time({
          nuMod = cv.glmnet(cbind(X[ids, ], Z[ids, ]), y[ids],
            keep = T,
            standardize = FALSE,
            family = muF,
            foldid = foldid[ids],
            parallel = parl,
            trace.it = tracer,
            alpha = alph
          )
        })
        if (noi) cat(glue("Outcome model {j} fitting took {elapExtract(time3)} seconds"), "\n")
        nuhat[, j] = predict(nuMod, cbind(X, Z), s = lamChoice)
        # oob for obs used to train model
        nuhat[ids, j] = fitGet(nuMod, lamChoice)
      }
    }
  }
  # all matrices will be N X n.avals
  amat = matrix(rep(a, n.avals), nrow = n, byrow = FALSE)
  smat = matrix(rep(s, n.avals), nrow = n, byrow = FALSE)
  rhomat = matrix(rep(rhohat, n.avals), nrow = n, byrow = FALSE)
  alevel = matrix(rep(as.numeric(avals), rep(n, n.avals)), nrow = n, byrow = FALSE)
  ymat = matrix(rep(y, n.avals), nrow = n, byrow = FALSE)
  return(list(
    # nuis fn fits
    muhat = muhat, pihat = pihat, rhomat = rhomat, nuhat = nuhat,
    # numbers
    amat = amat, smat = smat, alevel = alevel, ymat = ymat
  ))
}

# %% RANDOMFOREST
nuisGRF = function(y,
                   a,
                   s,
                   X,
                   n,
                   avals,
                   n.avals,
                   estimator = "AISW",
                   treatProb = NULL,
                   Z = NULL,
                   SepMu = TRUE,
                   noi = TRUE,
                   tune = "none") {
  # switches - fit only necessary nuisance functions
  if (estimator == "AISW") {
    fitmu = fitrho = fitpi = TRUE
  } else if (estimator == "ISW") {
    fitrho = fitpi = TRUE
    fitmu = FALSE
  } else if (estimator == "OM") {
    fitrho = fitpi = FALSE
    fitmu = TRUE
  } else if (estimator == "CW") {
    fitrho = fitmu = FALSE
    fitpi = TRUE
  } else if (estimator == "ACW") {
    fitrho = FALSE
    fitmu = fitpi = TRUE
  }
  # post-hoc mods conditional on inputs
  if (is.null(s)) {
    fitrho = FALSE; s = rep(1, n)
  }
  if (!is.null(treatProb)) fitpi = FALSE
  if (!is.null(Z)) {
    surrogates = TRUE
  } else {
    surrogates = FALSE
  }
  # containers
  muhat = as.data.frame(matrix(NA, nrow = n, ncol = n.avals))
  colnames(muhat) = paste("a", avals, sep = "")
  pihat = muhat
  nuhat = muhat
  # selection model
  if (fitrho) {
    rhohat = array(NA, dim = n)
    if (!surrogates) { # estimate selection score on X alone
      time3 = system.time({
        rhoMod = regression_forest(
          X,
          s == 1,
          tune.parameters = tune
        )
      })
      if (noi) cat(glue("RF: Selection model fitting took {elapExtract(time3)} seconds"), "\n")
    } else { # incorporate surrogates into selection score
      time3 = system.time({
        rhoMod = regression_forest(cbind(X, Z), s == 1)
      })
      if (noi) cat(glue("RF: Selection model fitting took {elapExtract(time3)} seconds \n"), "\n")
    }
    # out of bag predictions for s == 1
    # rhohat = predict(rhoMod, estimate.variance = FALSE)$predictions[, 2]
    rhohat = predict(rhoMod, X)$predictions
    rhohat = predict(rhoMod, estimate.variance = FALSE)$predictions # oob for obs used to train model
  } else {
    # all selection indicators and probabilities set to 1 - fits into AIPW
    rhohat = rep(1, n)
  }
  ######################################################################
  # propensity model
  if (fitpi) { # fit pscore
    for (j in 1:(n.avals - 1)) {
      if (noi) cat(glue("RF: Starting propensity model fit for treatment {j} seconds"), "\n")
      ids = s == 1
      time3 = system.time({
        piMod = regression_forest(
          X[ids, ],
          as.numeric(a == avals[j])[ids],
          tune.parameters = tune
        )
      })
      if (noi) cat(glue("RF: Propensity model fitting took {elapExtract(time3)} seconds"), "\n")
      # out of bag predictions of class probabilities
      pihat[, j] = predict(piMod, X)$predictions
      pihat[ids, j] = predict(piMod, estimate.variance = FALSE)$predictions # oob for obs used to train model
    }
    # last pscore is 1- sum(others)
    pihat[, n.avals] = 1 - rowSums2(as.matrix(pihat), na.rm = T)
  } else {
    if (inherits(treatProb, "matrix")) { # pscore is a matrix (already n X k)
      pihat = treatProb
    } else { # repeat each pscore vector k times
      pihat = matrix(rep(treatProb, n), nrow = n, byrow = FALSE)
    }
  }
  ######################################################################
  # outcome model
  if (fitmu) {
    if (SepMu) { # separate outcome model for each value of treatment
      for (j in 1:n.avals) {
        ids = a == avals[j] & s == 1
        # separate model for E[Y | A = a, X = x, S = 1] for each a
        time3 = system.time({
          muMod = regression_forest(
            X[ids, ],
            y[ids],
            tune.parameters = tune
          )
        })
        if (noi) cat(glue("RF: Outcome model fitting took {elapExtract(time3)} seconds"), "\n")
        muhat[, j] = predict(muMod, X)$predictions
        muhat[ids, j] = predict(muMod, estimate.variance = FALSE)$predictions # oob for obs used to train model
      }
    } else { # pooled model
      # single model for E[Y | A = a, X = x, S = 1]
      ids = s == 1
      time3 = system.time({
        sf = regression_forest(cbind(X, a)[ids, ], y[ids], tune.parameters = tune)
      })
      if (noi) cat(glue("RF: Outcome model fitting took {elapExtract(time3)} seconds"), "\n")
      preds.sf.oob = predict(sf, estimate.variance = FALSE)$predictions
      for (j in 1:n.avals) {
        muhat[, j] = predict(sf, cbind(X, as.numeric(avals[j])))$predictions
        muhat[a == avals[j], j] = preds.sf.oob[a == avals[j]]
      }
    }
    # outcome model with surrogates
    if (surrogates) {
      if (SepMu) { # separate outcome model for each value of treatment
        for (j in 1:n.avals) {
          ids = a == avals[j] & s == 1
          # separate model for E[Y | A = a, X = x, Z = z, S = 1] for each a
          time3 = system.time({
            nuMod = regression_forest(
              cbind(Z[ids, ], X[ids, ]),
              y[ids],
              tune.parameters = tune
            )
          })
          if (noi) cat(glue("RF: Surrogate outcome model fitting took {elapExtract(time3)} seconds"), "\n")
          nuhat[, j] = predict(nuMod, cbind(Z, X))$predictions
          nuhat[ids, j] = predict(nuMod, estimate.variance = FALSE)$predictions # oob for obs used to train model
        }
      } else { # pooled model
        # single model for E[Y | A = a, X = x, Z = z, S = 1]
        ids = s == 1
        time3 = system.time({
          sf = regression_forest(
            cbind(X, Z, a)[ids, ],
            y[ids],
            tune.parameters = tune
          )
        })
        if (noi) cat(glue("RF: Surrogate outcome model fitting took {elapExtract(time3)} seconds"), "\n")
        preds.sf.oob = predict(sf, estimate.variance = FALSE)$predictions
        for (j in 1:n.avals) {
          nuhat[, j] = predict(sf, cbind(X, Z, as.numeric(avals[j])))$predictions
          nuhat[a == avals[j], j] = preds.sf.oob[a == avals[j]]
        }
      }
    }
  }
  # all matrices will be N X n.avals
  amat = matrix(rep(a, n.avals), nrow = n, byrow = FALSE)
  smat = matrix(rep(s, n.avals), nrow = n, byrow = FALSE)
  rhomat = matrix(rep(rhohat, n.avals), nrow = n, byrow = FALSE)
  alevel = matrix(rep(as.numeric(avals), rep(n, n.avals)), nrow = n, byrow = FALSE)
  ymat = matrix(rep(y, n.avals), nrow = n, byrow = FALSE)
  return(list(
    # nuis fns
    muhat = muhat, pihat = pihat, rhomat = rhomat, nuhat = nuhat,
    # numbers
    amat = amat, smat = smat, alevel = alevel, ymat = ymat
  ))
}

# %% !TODO - nuisSL - extract doubly robust scores from sherlock (needs some tinkering because of missing data + selection score
