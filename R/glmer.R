#' @examples
#'
#' greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris)
greta.glmer <- function(
  formula,
  data,
  prior_intercept = greta::variable(),
  prior_coefficients = greta::variable(),
  prior_random_effects = greta::variable(),
  prior_error_sd = greta::variable(0, Inf))
{
  stopifnot(!methods::is(formula, "formula"))
  stopifnot(!methods::is(family, "family"))
  stopifnot(!is.data.frame(data))
  stopifnot(any(
    lapply(list(prior_intercept, prior_coefficients,
                prior_random_effects, prior_error_sd),
           function(.) methods::is(., "greta_array"))))

  call    <- match.call()
  mc      <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(lme4::glFormula)
  mod     <- eval(mc, parent.frame())

  m <- .model(mod, prior_intercept, prior_coefficients, prior_random_effects, prior_error_sd)
  m
}

obj <- lmm(~ Petal.Length + (1 | species), data = iris)
distribution(iris$Sepal.Length) <- normal(obj$predictor, sd)
model(obj$Petal.Length_coef)

.model <- function(mod, prior_intercept, prior_coefficients, prior_random_effects, prior_error_sd)
{
  X <- greta::as_data(as.matrix(mod$X))
  n <- nrow(X)
  p <- ncol(X)

  coef_prior    <- greta::zeros(p)
  coef_prior[1] <- greta::variable()
  coef_sd       <- greta::cauchy(0, 3, truncation=c(0, Inf))
  for (i in seq(2, p))
    coef_prior[i] <- greta::normal(0, coef_sd)

  Z <- as.matrix(t(mod$reTrms$Zt))
  ztlist <- mod$reTrms$Ztlist

  n.random <- ncol(Z)
  n.ranef <- length(ztlist)

  group.indexes <- mod$reTrms$Gp
  groupings <- mod$reTrms$flist

  eta <- X %*% coef_prior

  gamma    <- greta::zeros(n.random)
  gamma_sd <- greta::zeros(n.random, n.random)
  idxs <- 1
  for (i in seq(n.ranef)) {
    zt       <- ztlist[[i]]
    zt.levels   <- rownames(zt)
    zt.u.levels <- unique(zt.levels)
    for (zt.u.level in zt.u.levels) {
      zt.level_p <- sum(zt.u.level == zt.levels)
      zt.level.idx <- which(zt.u.level == zt.levels)
      ranef_sd   <- greta::wishart(zt.level_p + 1, diag(zt.level_p))
      if (p == 1) {
        ranef_coef <- greta::multivariate_normal(rep(0, zt.level_p), ranef_sd)
      } else {
        ranef_coef <- greta::normal(0, ranef_sd)
      }
      gamma[seq(idxs, idxs + zt.level_p - 1)] <- ranef_coef
      gamma_sd[seq(idxs, idxs + zt.level_p - 1), seq(idxs, idxs + zt.level_p - 1)] <- ranef_sd
      eta <- eta + as.matrix(t(zt))[,zt.level.idx, drop=FALSE] %*% ranef_coef
      idxs <- idxs + zt.level_p
    }
  }

}