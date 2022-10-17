# %% calibration functions : cloned from apoorvalal/renyiWeights
#' Entropy balancing by solving dual
#' @param X1m K-vector of target means
#' @param X0 NxK matrix of covariates
#' @param base.weight = [NULL] n-vector of baseline weights
#' @param coefs starting coefs for solution
#' @param maxIterations [200] stopping rule
#' @param constraint.tolerance [1] value for constraint threshold
#' @param printLevel [0, 1, 2, 3] 0 is silent, 1 reports success, 2 and 3 are noisy (for debugging)
#' @param sparsify [T/F] (in progress) run Newton-Raphson with sparse matrix classes from Matrix package
#' @export

eb_solve_dual = function(X1m, X0,
                         coefs = NULL, base.weight = NULL,
                         maxIterations = 200L,
                         constraint.tolerance = 1,
                         printLevel = 0) {
  if (is.null(coefs)) coefs = rep(0, ncol(X0))
  if (is.null(base.weight)) base.weight = rep(1, nrow(X0))
  converged = FALSE
  for (iter in 1:maxIterations) {
    # Z is a R-vector of solution coefficients (coefs)
    #########################################################
    # unnormalised weights
    weights.temp = c(exp(X0 %*% coefs))
    weights.ebal = weights.temp * base.weight
    ### gradient construction
    X0.agg = c(weights.ebal %*% X0)
    # ∇_Z = M - CW
    gradient = X0.agg - X1m
    # minimum reached
    if (max(abs(gradient)) < constraint.tolerance) {
      converged = TRUE
      break
    }
    # noisy
    if (printLevel >= 2) {
      cat(
        "Iteration", iter, "maximum deviation is =",
        format(max(abs(gradient)), digits = 4), "\n"
      )
    }
    hessian = t(X0) %*% (weights.ebal * X0)
    Coefs = coefs
    newton = solve(hessian, gradient)
    coefs = coefs - newton

    # step length - optimal step length section
    loss.new = line.searcher(
      Base.weight = base.weight, X0 = X0,
      X1m = X1m, coefs = coefs, Newton = newton, ss = 1
    )
    loss.old = line.searcher(
      Base.weight = base.weight, X0 = X0,
      X1m = X1m, coefs = Coefs, Newton = newton, ss = 0
    )

    if (is.na(loss.new) | is.na(loss.old)) {
      stop(
        "Optimization ran into problems. Loss function values are",
        "\n New:", loss.new,
        "\n Old:", loss.old,
        "\n Constraint tolerance is", constraint.tolerance,
        "\n Increase constraint tolerance or trim moment conditions problem"
      )
    }

    if (printLevel >= 3) cat("new loss", loss.new, "old loss=", loss.old, "\n")
    if (loss.old <= loss.new) {
      ss.out = optimize(line.searcher,
        lower = .00001, upper = 1, maximum = FALSE,
        Base.weight = base.weight, X0 = X0, X1m = X1m,
        coefs = Coefs, Newton = newton
      )
      if (printLevel >= 3) {
        cat("LS Step Length is ", ss.out$minimum, "\n")
      }
      if (printLevel >= 3) {
        cat("Loss is", ss.out$objective, "\n")
      }
      coefs = Coefs - ss.out$minimum * solve(hessian, gradient)
    }
  }
  # step out of loop
  if (printLevel >= 1 && converged) cat("Converged within tolerance \n")

  return(
    list(
      maxdiff = max(abs(gradient)),
      coefs = coefs,
      Weights.ebal = weights.ebal,
      converged = converged
    )
  )
}

# %% internal for step length
line.searcher = function(Base.weight, X0, X1m, coefs, Newton, ss) {
  weights.temp = c(exp(X0 %*% (coefs - (ss * Newton))))
  weights.temp = weights.temp * Base.weight
  X0.agg = c(weights.temp %*% X0)
  maxdiff = max(abs(X0.agg - X1m))
  return(maxdiff)
}
