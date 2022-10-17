# %% estimate treatment effects and marginal means
# internal
estQois = function(ifvals, avals, denom = NULL, noi = TRUE) {
  n = dim(ifvals)[1]
  n.avals = dim(ifvals)[2]
  if (!is.null(denom)) n = denom
  # check for missingness in marginal means
  if (is.na(sum(apply(ifvals, 2, sum))))
    warning("Missing values present in marginal means. This may happen due to extreme propensity scores. Subset ifvals and recompute.")
  # Counterfactual means
  est = apply(ifvals, 2, sum, na.rm = T)
  est = est / n # this n can be different from nrow(ifvals)
  se = apply(ifvals, 2, sd, na.rm = T) / sqrt(n)
  ci.ll = est - 1.96 * se
  ci.ul = est + 1.96 * se
  pval = round(2 * (1 - pnorm(abs(est / se))), 3)
  paste("E{Y(", avals, ")}")
  res1 = data.frame(
    parameter = paste("E{Y(", avals, ")}", sep = ""),
    est, se, ci.ll, ci.ul, pval
  )
  # distance  between all possible Y^d
  signdist = function(x) c(as.dist(outer(x, x, "-")))
  ifvals2 = t(apply(ifvals, 1, signdist))
  if (n.avals == 2) ifvals2 = t(ifvals2)
  tmp = expand.grid(1:n.avals, 1:n.avals)
  tmp2 = tmp[tmp[, 1] > tmp[, 2], ]
  contrasts = apply(cbind(avals[tmp2[, 1]], avals[tmp2[, 2]]), 1,
    paste,
    collapse = ")-Y("
  )
  # all pairwise contrasts
  contrasts = paste("E{Y(", contrasts, ")}", sep = "")
  est2 = apply(ifvals2, 2, sum, na.rm = T)
  est2 = est2 / n # this n can be different from nrow(ifvals)
  se2 = apply(ifvals2, 2, sd, na.rm = T) / sqrt(n)
  ci.ll2 = est2 - 1.96 * se2
  ci.ul2 = est2 + 1.96 * se2
  pval2 = round(2 * (1 - pnorm(abs(est2 / se2))), 3)
  res2 = data.frame(
    parameter = contrasts, est = est2, se = se2,
    ci.ll = ci.ll2, ci.ul = ci.ul2, pval = pval2
  )
  # concat marginal means and
  res = rbind(res1, res2)
  rownames(res) = NULL
  if (noi) print(res)
  return(res)
}
