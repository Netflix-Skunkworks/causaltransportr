# %% ####################################################
#' Omnibus function for ATE estimation for generalization and transportation
#' @description Augmented IPW, generalization, or transport estimator with ML and cross-fitting for nuisance functions. Imputes counterfactual outcome for each observation i under each treament a as
#' \deqn{Y^a = \omega \frac{A = a}{\pi^a (X)} (Y - \hat{\mu}^a(X)) + \mu^a(X) }
#' Where \eqn{\omega} is 1 for all observations under no sample selection, and therefore this is the doubly-robust Augmented Inverse Propensity Weighting (AIPW) estimator. When S is supplied, the argument in 'target' is used to fit either the generalization or transportation estimator, which corresponds with \eqn{\omega = S/\rho(X)} and \eqn{\omega = (S (1-\rho(X)) /\rho(X)} respectively. When a surrogate vector \eqn{Z} is supplied, an additional residual piece \eqn{a /\pi(X)(\hat{\nu}(X, Z) -  \hat{\mu}(X))}  is added to the influence function.
#' Average treatment effects are defined as averages of differences between counterfactual outcomes \eqn{Y^a - Y^{a'}}.
#' @param y outcome vector (may contain missings ; missings must correspond with s = 0)
#' @param a treatment vector (no missings; can be relaxed with some tinkering)
#' @param X covariate matrix (no missings)
#' @param s selection vector, NULL by default (no missings, 1 corresponds with nonmissing y; 0 corresponds with missing y). May be omitted when the target is "insample" .
#' @param treatProb propensity score vector (of length n_treatment) or matrix (n_treatment X n_obs), where latter is for covariate adaptive designs; must sum to 1. NULL by default, so pscore is fitted. When provided, no propensity score is fit. With discrete covariates, estimated propensity score is advisable even if treatment was randomized.
#' @param Z surrogate matrix, NULL by default (no missings). When nonmissing, the surrogate influence function (Kallus and Mao 2020) is used to compute treatment effects.

#' @param nuisMod one of c("rlm", "rf") : choose how to fit nuisance functions (cross-fit).
#' @param target one of c("generalize", "transport", "insample") estimand to target. "generalize" generalizes (quasi)experimental estimates from the complete data (S == 1) to the overall sample (S == 0 or S == 1). "transport" transports estimates from the S == 1 sample to the S == 0 sample. "insample" estimates causal effects in the S == 1 sample (i.e. conventional quasi/experimental estimation).
#' @param estimator one of c("AISW", "ISW", "OM", "CW", "ACW"). The default is the augmented inverse selection weighting estimator, which augments the inverse selection weighting estimator (ISW) with an outcome model (OM). ACW does the same with calibration weights (CW), which fit a set of entropy balancing weights that reweights the sample to match target sample moments.
#' @param hajekize boolean for whether to divide the inverse probability weights term for each treatment level by the sum of weights in that treatment level. This guards against instability from very large weights from extremely small selection or propensity scores.
#' @param separateMus boolean for whether to fit separate outcome models for each treatment group or a single pooled model. The former is recommended and is the default, but a pooled model may be fit when data is scarce / computation is burdensome.

#' @param glmnet_lamchoice choice of lambda (shrinkage parameter) for regularized linear regressions. Only relevant when nuisMod == "rlm"
#' @param glmnet_alpha in [0, 1], choice of alpha in glmnet. 1 (default) corresponds with L1 regularization (LASSO) and 0 corresponds with L2 regularization (ridge), while intermediate values correspond with a mix of the two (elastic net)
#' @param glmnet_rho_family  GLM family for selection model. "binomial" by default but can be safely switched to "gaussian" for linear probability models with discrete covariates for faster compute
#' @param glmnet_pi_family  GLM family for propensity model. "binomial" by default but can be safely switched to "gaussian" for linear probability models with discrete covariates for faster compute
#' @param glmnet_mu_family  GLM family for outcome model. Gaussian by default.
#' @param glmnet_parl Boolean for parallelization in glmnet. Need to enable parallelized cluster beforehand.

#' @param grf_tuneRf Tune rf hyperparameters? Passed to grf's regression forest. Use 'all' for hyperparameter tuning.
#' @param noi boolean for printing marginal means and causal contrasts table (it gets returned anyway). Off by default.

#' @return list containing treatment effects table , nuisance function estimates, and influence function values

#' @export
#' @examples
#' \dontrun{
#' df = selDGP(n = 1000)
#' df$tau |> mean() # true treatment effect
#' estimation with lasso
#' ateGT(df$y, df$a, df$X, df$s, noi = FALSE, nuisMod = 'rlm')$res
#' estimation with rf
#' ateGT(df$y, df$a, df$X, df$s, noi = FALSE, nuisMod = 'rf')$res
#' }
#' @references Bia, M., M. Huber, and L. Lafférs. (2020): “Double Machine Learning for Sample Selection Models,” arXiv [econ.EM],.
#' @references Dahabreh, I. J., S. E. Robertson, E. J. Tchetgen, E. A. Stuart, and M. A. Hernán. (2019): “Generalizing causal inferences from individuals in randomized trials to all trial-eligible individuals,” Biometrics, 75, 685–94.
#' @references Hirshberg, D. A., A. Maleki, and J. R. Zubizarreta. (2019): “Minimax Linear Estimation of the Retargeted Mean,” arXiv [math.ST],.
#' @references Kallus, N., and X. Mao. (2020): “On the Role of Surrogates in the Efficient Estimation of Treatment Effects with Limited Outcome Data,” arXiv [stat.ML],.

ateGT = function( y,
                  a,
                  X,
                  s = NULL, # selection indicator - null by default
                  treatProb = NULL, # treatment probability - if supplied no propensity score fit
                  Z = NULL, # null by default; if non-null, incorporates surrogates automatically
                  nuisMod = c("rlm", "rf"), # model to fit double robust score - can add more
                  target = c("generalize", "transport", "insample"),
                  estimator = c("AISW", "ISW", "OM", "CW", "ACW"),
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
                  grf_tuneRf = "none",
                  # misc
                  noi = FALSE) {
    # housekeeping
    nuisMod = match.arg(nuisMod)
    target = match.arg(target)
    estimator = match.arg(estimator)
    nosurrogates = is.null(Z)

    n = dim(X)[1]
    avals = names(table(a))
    n.avals = length(avals)

    if (target %in% c("generalize", "transport") & is.null(s)) {
        stop("Cannot transport without selection indicator")
    }

    if(!is.null(treatProb)){
        if(noi) cat("Propensity scores passed manually. No pscore model will be fit. Note that using the estimated propensity score is more efficient even in randomized experiments. \n")
        # turn it into a vector if only p(a = 1) is passed
        if(length(treatProb) == 1) treatProb = c(treatProb, 1 - treatProb)
    }

    ############################################################
    # fit nuisance functions
    ############################################################
    if (nuisMod == "rlm") { # calls glmnet
        nuis = nuisRLM(y, a, s, X, n, avals, n.avals,
            estimator = estimator,
            treatProb = treatProb,
            Z = Z,
            SepMu = separateMus,
            # glmnet parameters
            lamChoice = glmnet_lamchoice,
            alph = glmnet_alpha,
            rhoF = glmnet_rho_family,
            piF  = glmnet_pi_family,
            muF  = glmnet_mu_family,
            parl = glmnet_parl,
            # others
            noi = noi
        )
    } else if (nuisMod == "rf") {
        nuis = nuisGRF(y, a, s, X, n, avals, n.avals,
            estimator = estimator,
            treatProb = treatProb,
            Z = Z,
            SepMu = separateMus,
            noi = noi,
            tune = grf_tuneRf
        )
    } else {
        stop("Nuisance function model not supported.")
    }
    ##################################################################################
    # estimation: construct influence functions to average later for causal parameters
    ##################################################################################
    if(nosurrogates){ # no surrogates - default
        if (target == "generalize") {
            ############################################################
            # generalization
            ############################################################
            # (marginalise over s, s = 0, s = 1)
            if (estimator == "AISW") {
                # selection indic in numerator and selection proba in denominator
                ipw_piece = (
                    ((nuis$amat == nuis$alevel) * nuis$smat) / (nuis$rhomat * nuis$pihat)
                ) * (nuis$ymat - nuis$muhat)
                # indices of observations with missing ipw because y is missing
                ymissing_indices = apply(ipw_piece, 1, function(x) any(is.na(x)))
                # set them to 0 (rely on outcome modeling for those obs)
                ipw_piece[ymissing_indices, ] = 0
                if (hajekize) ipw_piece = ipw_piece / colSums(nuis$pihat)
                # add in outcome model piece
                ifvals = ipw_piece + nuis$muhat
            } else if (estimator == "OM") {
                ifvals = nuis$muhat
            } else if (estimator == "ISW") {
                # selection indic in numerator and selection proba in denominator
                ipw_piece = (
                    ((nuis$amat == nuis$alevel) * nuis$smat) / (nuis$rhomat * nuis$pihat)
                ) * nuis$ymat
                # indices of observations with missing ipw because y is missing
                ymissing_indices = apply(ipw_piece, 1, function(x) any(is.na(x)))
                # set them to 0 (rely on outcome modeling for those obs)
                ipw_piece[ymissing_indices, ] = 0
                if (hajekize) ipw_piece = ipw_piece / colSums(nuis$pihat)
                ifvals = ipw_piece
            } else if (estimator == "CW") {
                s1 = (s == 1)
                n1 = sum(s1)
                # solve for weights that set covariate means to target means in whole sample
                ebwts = eb_solve_dual(
                    c(1, colMeans2(X)), # target is overall sample mean
                    cbind(1, X[s1, ]) # reweighting sample is source (s == 1)
                )$Weights.ebal
                # need to scale each weight by n1 so that averaging works
                wtmat = matrix(rep(ebwts * n1, n.avals), nrow = n1, byrow = FALSE)
                # ipw
                ipw_piece = wtmat * (
                    (nuis$amat[s1, ] == nuis$alevel[s1, ]) / nuis$pihat[s1, ]) *
                    nuis$ymat[s1, ]
                ifvals = ipw_piece
            } else if (estimator == "ACW") {
                # first do calibration for s = 1 group
                s1 = (s == 1)
                n1 = sum(s1)
                # solve for weights that set covariate means to target means in whole sample
                ebwts = eb_solve_dual(
                    c(1, colMeans(X)), # target is overall sample mean
                    cbind(1, X[s1, ]) # reweighting sample is source (s == 1)
                )$Weights.ebal
                # need to scale each weight by n1 so that averaging works
                wtmat = matrix(rep(ebwts * n1, n.avals), nrow = n1, byrow = FALSE)
                # balwt fit
                balfit = wtmat * (
                    (nuis$amat[s1, ] == nuis$alevel[s1, ]) / nuis$pihat[s1, ]) *
                    (nuis$ymat[s1, ] - nuis$muhat[s1, ])
                # initialise container of all obs
                ipw_piece = matrix(0, nrow = nrow(X), ncol = n.avals)
                # fill in obs only for s = 1
                ipw_piece[s1, ] = as.matrix(balfit)
                # add in outcome model piece
                ifvals = ipw_piece + nuis$muhat
            } else {
                stop("Estimator not supported.")
            }
            # estQOIs called once for all  estimators (different ifs)
            res = estQois(ifvals, avals, noi = noi)
        } else if (target == "transport") {
            ############################################################
            # transportation
            ############################################################
            # estQOIs called inside each estimator block because of potentially different denominators
            if (estimator == "AISW") {
                # selection indic and 1 - selection probability in numerator and selection proba in denominator
                ipw_piece = (((nuis$amat == nuis$alevel) * nuis$smat * nuis$rhomat) /
                    ((1 - nuis$rhomat) * nuis$pihat)
                ) * (nuis$ymat - nuis$muhat)
                # indices of observations with missing ipw because y is missing
                ymissing_indices = apply(ipw_piece, 1, function(x) any(is.na(x)))
                # set them to 0 (rely on outcome modeling for those obs)
                ipw_piece[ymissing_indices, ] = 0
                if (hajekize) ipw_piece = ipw_piece / colSums(nuis$pihat)
                # add in outcome model piece (only kicks in for s = 0 group)
                ifvals = ipw_piece + (1 - nuis$smat) * nuis$muhat
                # estimation: pass different denominator
                res = estQois(ifvals, avals, denom = sum(s == 0), noi = noi)
            } else if (estimator == "OM") {
                # only average over s = 0 group
                ifvals = nuis$muhat[s == 0, ]
                res = estQois(ifvals, avals, noi = noi)
            } else if (estimator == "ISW") {
                ipw_piece = (
                    ((nuis$amat == nuis$alevel) * nuis$smat * (1 - nuis$rhomat)) /
                        (nuis$rhomat * nuis$pihat)
                ) * nuis$ymat
                # indices of observations with missing ipw because y is missing
                ymissing_indices = apply(ipw_piece, 1, function(x) any(is.na(x)))
                # set them to 0 (rely on outcome modeling for those obs)
                ipw_piece[ymissing_indices, ] = 0
                if (hajekize) ipw_piece = ipw_piece / colSums(nuis$pihat)
                ifvals = ipw_piece
                res = estQois(ifvals, avals, denom = sum(s == 0), noi = noi)
            } else if (estimator == "CW") {
                s1 = (s == 1)
                n1 = sum(s1)
                # solve for weights that set covariate means to target means in target sample
                ebwts = eb_solve_dual(
                    c(1, colMeans(X[s == 0, ])), # target mean vector
                    cbind(1, X[s1, ]) # reweighting sample is source (s == 1)
                )$Weights.ebal
                # need to scale each weight by n1 so that averaging works
                wtmat = matrix(rep(ebwts * n1, n.avals), nrow = n1, byrow = FALSE)
                # ipw piece is truncated to be just n1 rows, so no need for normalisation
                ipw_piece = wtmat * (
                    (nuis$amat[s1, ] == nuis$alevel[s1, ]) / nuis$pihat[s1, ]) *
                    nuis$ymat[s1, ]
                ifvals = ipw_piece
                res = estQois(ifvals, avals, noi = noi)
            } else if (estimator == "ACW") {
                # first do calibration for s = 1 group
                s1 = (s == 1)
                n1 = sum(s1)
                # solve for weights that set covariate means to target means in target sample
                ebwts = eb_solve_dual(
                    c(1, colMeans(X[s == 0, ])), # target mean vector
                    cbind(1, X[s1, ]) # reweighting sample is source (s == 1)
                )$Weights.ebal
                # need to scale each weight by n1 so that averaging works
                wtmat = matrix(rep(ebwts * n1, n.avals), nrow = n1, byrow = FALSE)
                # balwt fit
                balfit = wtmat * (
                    (nuis$amat[s1, ] == nuis$alevel[s1, ]) / nuis$pihat[s1, ]) *
                    (nuis$ymat[s1, ] - nuis$muhat[s1, ])
                # initialise container of all obs
                ipw_piece = matrix(0, nrow = nrow(X), ncol = n.avals)
                # fill in obs only for s = 1
                ipw_piece[s1, ] = as.matrix(balfit)
                # add in outcome model piece
                ifvals = ipw_piece + nuis$muhat
                ifvals = ipw_piece + (1 - nuis$smat) * nuis$muhat
                # estimation: pass different denominator
                res = estQois(ifvals, avals, denom = sum(s == 0), noi = noi)
            } else {
                stop("Estimator not supported.")
            }
        } else if (target == "insample") { # just fit AIPW
            if (sum(is.na(y)) > 0) stop("outcome contains missings; fit with sample correction")
            if (estimator == "AISW") {
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
            res = estQois(ifvals, avals, noi = noi)
        } else {
            stop("Target not supported. Must be in (generalize, transport, insample)")
        }
    } else { # use surrogates
        if (target == "generalize") {
            if (estimator == "AISW"){
                # Do AIPWS
                ipw_piece = ((nuis$amat == nuis$alevel) * nuis$smat * (nuis$ymat - nuis$muhat)) /
                    (nuis$rhomat * nuis$pihat)
                # indices of observations with missing ipw because y is missing
                ymissing_indices = apply(ipw_piece, 1, function(x) any(is.na(x)))
                # set them to 0 (rely on outcome modeling for those obs)
                ipw_piece[ymissing_indices, ] = 0
                if (hajekize) ipw_piece = ipw_piece / colSums(nuis$pihat)
                # ipw surrogate residual
                ipw_piece2 = ((nuis$amat == nuis$alevel) * (nuis$nuhat - nuis$muhat)) / nuis$pihat
                # add in outcome model piece *and* additional residual
                ifvals = nuis$muhat + ipw_piece + ipw_piece2
            } else if(estimator == "OM"){
                ifvals = nuis$nuhat
            } else {
                stop(glue("{estimator} not supported with surrogates. Use AISW or OM"))
            }
        } else if (target == "transport") {
            if (estimator == "AISW"){
                # Do AIPWS
                ipw_piece = ((nuis$amat == nuis$alevel) * nuis$smat * (1 - nuis$rhomat) * (nuis$ymat - nuis$muhat)) /
                    (nuis$rhomat * nuis$pihat)
                # indices of observations with missing ipw because y is missing
                ymissing_indices = apply(ipw_piece, 1, function(x) any(is.na(x)))
                # set them to 0 (rely on outcome modeling for those obs)
                ipw_piece[ymissing_indices, ] = 0
                if (hajekize) ipw_piece = ipw_piece / colSums(nuis$pihat)
                # ipw surrogate residual
                ipw_piece2 = ((nuis$amat == nuis$alevel) * (nuis$nuhat - nuis$muhat)) / nuis$pihat
                # add in outcome model piece *and* additional residual
                ifvals = nuis$muhat + ipw_piece + ipw_piece2
            } else if (estimator == "OM") {
                ifvals = nuis$nuhat[s == 0,]
            } else {
                stop(glue("{estimator} not supported with surrogates. Use AISW or OM"))
            }
        } else if (target == "insample"){
            stop("Surrogates not supported for insample estimation. Call the function without s and Z.")
        } else {
            stop("Target not supported. Must be in (generalize, transport, insample)")
        }
        res = estQois(ifvals, avals, noi = noi)
    }

    outobj = list(
        res       = res,
        nuis      = nuis,
        ifvals    = as.data.frame(ifvals),
        nuisMod   = nuisMod,
        target    = target,
        estimator = estimator,
        n         = dim(X)[1],
        nA        = length(avals)
    )
    class(outobj) = "ategt"
    return(invisible(outobj))
}

#' summary method - prints treatment effects and marginal means
#' @param fitted ategt object
#' @export
summary.ategt = function(obj) print(obj$res)

#' SMD balance plot method - plots balance across covariates and gains from reweighting
#' @param fit ategt object
#' @param covariate matrix
#' @examples
#' sim = selDGP(selF = \(x) -3 * x[5] + x[2] + 3 * x[10])
#' fit = with(sim, ateGT(y = y, a = a, X = X, s = s))
#' fit %>% summary
#' mean(sim$tau)
#' plot(fit, sim$X)
#' @export
plot.ategt = function(fit, X){
    require(ggplot2)
    # check for names
    if(is.null(colnames(X))) colnames(X) = paste0("X", 1:ncol(X))
    Xn = colnames(X)
    # do SMD calculations in datatable
    baltab = data.table(a = fit$nuis$smat[, 2], pscore = fit$nuis$rhomat[,2], X)
    # raw mean
    groupMeans = baltab[, lapply(.SD, mean), a, .SDcols = Xn]
    # pooled SD
    sd1 = baltab[, lapply(.SD, var), a, .SDcols = Xn][, -1] |> colMeans() |> sqrt()
    std_mean_diffs = abs(-colDiffs(as.matrix(groupMeans[, -1]))/ sd1)
    # weighted mean
    groupMeans2 = baltab[, lapply(.SD, weightedMean, w = pscore), a, .SDcols = Xn]
    # pooled SD
    sd2 = baltab[, lapply(.SD, weightedVar, pscore), a, .SDcols = Xn][, -1] |> colMeans() |> sqrt()
    std_mean_diffs2 = abs(-colDiffs(as.matrix(groupMeans2[, -1]))/ sd2)
    # plot it
    plotdf = data.table(Xn, t(std_mean_diffs), t(std_mean_diffs2));
    colnames(plotdf) = c("X", "Raw", "Wtd")
    melted = melt(plotdf, id = "X")
    # plot
    ggplot(melted, aes(X, value, colour = variable, group = variable)) +
        geom_point(alpha = 0.4, size = 2) + coord_flip() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), legend.pos = "top") +
        geom_hline(aes(yintercept = 0.2), linetype = 'dotted') +
        geom_hline(aes(yintercept = 0.1), linetype = 'dashed') +
        ylim(c(0, max(melted$value) + 0.01)) +
        labs(title = "Standardized Mean Difference in Covariates",
                subtitle = "abs gap between Complete and Missing Samples")
}
