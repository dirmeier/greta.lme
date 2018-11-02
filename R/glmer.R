#' @title Initialize the linear predictor for generalized linear effect models
#'
#' @description \code{greta.glmer} creates the linear predictor $eta$ required
#' for fitting generalized linear mixed models. The linear predictor has the
#' form: \deqn{\eta = X \beta + Z \gamma,}
#' where \eqn{X} is the fixed effects design matrix, \eqn{\beta} the fixed
#' effects, \eqn{Z} the random effects design matrix and \eqn{\gamma} the
#' random effects.
#'
#' @export
#' @import greta
#' @importFrom lme4 glFormula
#' @importFrom methods is
#'
#' @param formula  an \code{lme4}-style formula
#' @param data  the data set containing the variables described in
#'  \code{formula}
#' @param  prior_intercept  a \code{\link{greta_array}} random variable
#'  that specifies the prior for the intercept \eqn{\beta_0} for the
#'  fixed effects design
#'  matrix. If \code{NULL}, the distribution for the intercept is chosen as
#'  \deqn{p(\beta) ∝ const.} The intercept prior
#'  needs to have dimension 1.
#' @param  prior_coefficients  a \code{\link{greta_array}} random
#'  variable that specifies the prior for the coefficients\eqn{\beta} for the
#'  fixed effects
#'  design matrix. If \code{NULL}, the distribution for the coefficients is
#'  chosen as \deqn{p(\beta_i) ~ N(0, \sigma),} and
#'   \deqn{p(\sigma) ~ Half-Cauchy(0, ∞).}
#'  The coefficients prior
#'  need to have the same dimensionalty as the number of columns of
#'  your fixed effects design matrix (without the intercept).
#' @param  prior_random_effects  a \code{greta_array} random
#'  variable that specifies the prior for the random effects \eqn{\gamma} for
#'  the random effects design matrix. If \code{NULL}, the distribution for the
#'   random effects is chosen as
#'   \deqn{p(\gamma) ~ N(0, \tau),} and
#'   \deqn{p(\tau) ~ Inv-Wishart}
#'  If \code{prior_random_effects} is provided
#'  it needs to have the same dimensionalty as the
#'  number of columns of \eqn{Z}
#'  random effect design matrix. These are generally not trivial to construct.
#'  Calling \code{glFormula(formula, data)} creates the specified
#'  random effects matrix \eqn{Z}. From this the dimensionality of
#'  \code{prior_random_effects} can be inferred.
#'
#'
#' @return Returns a list of class \code{greta.glmer} with the following
#'  elements:
#' \item{predictor }{ the linear predictor \eqn{\eta}}
#' \item{X }{ the fixed effects design matrix}
#' \item{x.eta }{ the linear predictor \eqn{X\beta}}
#' \item{coef }{ the prior for the coefficients \eqn{\beta}}
#' \item{Z }{ the random effects design matrix}
#' \item{z.eta }{ the linear predictor \eqn{Z\gamma}}
#' \item{gamma }{ the prior for the random effects \eqn{\gamma}}
#' \item{Ztlist }{ a list of random effect terms used for initializing \eqn{\gamma}}
#' \item{grp.vars }{ the variables used for grouping}
#' \item{call }{ the function call used}
#' \item{formula }{ the formula used}
#'
#' @examples
#' # create a random intercept model
#' greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris)
#'
#' # creates a random slope model
#' greta.glmer(Sepal.Length ~ Sepal.Width + (Sepal.Width | Species), iris)
#'
#' # creates a random slope model with multiple random effect terms
#' greta.glmer(Sepal.Length ~ Sepal.Width + (Sepal.Width | Species) + (Petal.Width | Species),
#'             iris)
#'
#' # creates a random slope model with strong regularizing prior
#' greta.glmer(Sepal.Length ~ Sepal.Width + (Sepal.Width | Species),
#'            iris,
#'            prior_random_effects = greta::normal(0, 1, dim=6))
#'
#' # creates a random slope model with flat coefficient prior
#' greta.glmer(Sepal.Length ~ Sepal.Width + (Sepal.Width | Species),
#'            iris,
#'            prior_coefficients = greta::variable(dim=1))
#'
#' # creates a random slope model with normal intercept prior
#' greta.glmer(Sepal.Length ~ Sepal.Width + (Sepal.Width | Species),
#'            iris,
#'            prior_intercept = greta::normal(0, 1))
#'
greta.glmer <- function(
  formula,
  data,
  prior_intercept = NULL,
  prior_coefficients = NULL,
  prior_random_effects = NULL)
{
  stopifnot(methods::is(formula, "formula"))
  stopifnot(is.data.frame(data))

  mc      <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(lme4::glFormula)
  mc$prior_intercept <- mc$prior_coefficients <- mc$prior_random_effects <- NULL
  mod     <- eval(mc, parent.frame())

  ret <- .model(mod, prior_intercept, prior_coefficients, prior_random_effects)
  ret$call    <- match.call()
  ret$formula <- mod$formula
  class(ret) <- "greta.glmer"

  ret
}


#' @noRd
.model <- function(mod, prior_intercept, prior_coefficients, prior_random_effects)
{
  coef_list  <- .coef_priors(mod, prior_intercept, prior_coefficients)
  ranef_list <- .ranef_priors(mod, prior_random_effects)
  eta        <- coef_list$x.eta + ranef_list$z.eta

  res           <- c(coef_list, ranef_list)
  res$formula   <- mod$mormula

  ret <- c(
    list(predictor = eta),
    res
  )

  ret
}


#' @noRd
#' @importFrom methods is
#' @importFrom greta normal variable zeros multivariate_normal wishart as_data
.ranef_priors <- function(mod, prior_random_effects)
{
  Z        <- greta::as_data(t(as.matrix(mod$reTrms$Zt)))
  Ztlist   <- mod$reTrms$Ztlist
  Zp       <- ncol(Z)
  grp.vars <- mod$reTrms$flist

  if (!is.null(prior_random_effects))
  {
    stopifnot(methods::is(prior_random_effects, "greta_array"))
    desc <- attr(prior_random_effects, "node")$description()
    if (!length(grep("normal", desc)))
      stop("Random effects prior does not follow a normal distribution",
           call.=FALSE)
    check.dimensionality(dim(prior_random_effects), c(Zp, 1))
  }

  gamma_sd <- greta::zeros(Zp, Zp)
  gamma <- if (is.null(prior_random_effects)) {
    greta::zeros(Zp)
  } else {
    prior_random_effects
  }

  eta <- 0
  idxs  <- 1
  for (i in seq_len(length(Ztlist)))
  {
    # current random term
    zt          <- Ztlist[[i]]
    # grouping levels for that term
    zt.levels   <- rownames(zt)
    # unique grouping levels, i.e. how which groups for this term are there
    # TODO: can this be done nicer using `lme4::glFormula` elements?
    zt.u.levels <- unique(zt.levels)

    # iterate over the different groups
    for (zt.u.level in zt.u.levels)
    {
      # dimensionality of the random effects for the current group
      zt.level_p   <- sum(zt.u.level == zt.levels)
      zt.level.idx <- which(zt.u.level == zt.levels)

      # we just use a wishart, because it creates an psd matrix
      ## changing to ranef_sd here throws an error
      # ranef_sd <- solve(greta::wishart(zt.level_p + 1, diag(zt.level_p)))
      ranef_sd <- base::diag(zt.level_p)
      ranef_coef <- if (zt.level_p > 1) {
        greta::multivariate_normal(t(rep(0, zt.level_p)), ranef_sd)
      } else {
        greta::normal(0, ranef_sd)
      }

      # if random effects prior is not provided, create a random variable and
      # an rv for the covariance
      if (is.null(prior_random_effects)) {
        gamma_sd[seq(idxs, idxs + zt.level_p - 1),
                 seq(idxs, idxs + zt.level_p - 1)] <- ranef_sd
        gamma[seq(idxs, idxs + zt.level_p - 1)] <- ranef_coef
        stopifnot(ncol(Z) == nrow(gamma))
      }

      # add design matrix * random effect to the linear predictor
      z    <- greta::as_data(t(as.matrix(zt))[, zt.level.idx, drop=FALSE])
      eta  <- eta + z %*% gamma[seq(idxs, idxs + zt.level_p - 1)]
      idxs <- idxs + zt.level_p
    }
  }

  ret <- list(
    Z = Z,
    z.eta = eta,
    gamma = gamma)
  if (is.null(prior_random_effects))
    ret$gamma_sd <- gamma_sd
  ret$Ztlist <- Ztlist
  ret$grp.vars <- grp.vars

  ret
}


#' @noRd
#' @importFrom methods is
#' @importFrom greta normal variable zeros cauchy as_data
.coef_priors <- function(mod, prior_intercept, prior_coefficients)
{

  X <- greta::as_data(as.matrix(mod$X))
  p <- ncol(X)
  el <- list()

  coef_prior <- greta::zeros(p)

  if (is.null(prior_intercept))
  {
    coef_prior[1] <- greta::variable()
  }
  else
  {
    stopifnot(methods::is(prior_intercept, "greta_array"))
    check.dimensionality(dim(prior_intercept), c(1, 1))
    coef_prior[1] <- prior_intercept
  }

  if (is.null(prior_coefficients))
  {
    coef_sd <- greta::cauchy(0, 3, truncation=c(0, Inf))
    for (i in seq(2, p))
      coef_prior[i] <- greta::normal(0, coef_sd)
    el$coef_sd <- coef_sd
  }
  else
  {
    stopifnot(methods::is(prior_coefficients, "greta_array"))
    check.dimensionality(dim(prior_coefficients), c(p - 1, 1))
    coef_prior[seq(2, p)] <- prior_coefficients
  }

  eta <- X %*% coef_prior

  ret <- c(
    list(
      X = X,
      x.eta = eta,
      coef = coef_prior),
    el)

  ret
}