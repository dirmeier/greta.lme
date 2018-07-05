
greta.glmer <- function(
  formula,
  data,
  family = gaussian(),
  prior_intercept = NULL,
  prior_coefficients = NULL,
  prior_random_effects = NULL,
  prior_error_sd = NULL)
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

  m <- .model(mod, family, prior_intercept, prior_coefficients, prior_random_effects, prior_error_sd)
  m
}

.model <- function(mod, family, prior_intercept, prior_coefficients, prior_random_effects, prior_error_sd)
{
  X <- greta::as_data(mod$X)
  n <- nrow(X)
  p <- ncol(X)

  int <- if (is.null(prior_intercept)) {
    greta::uniform(-1000, 1000, dim=1)
  } else {
      prior_intercept
  }
  coef <- if (is.null(prior_coefficients)) {
    greta::uniform(-1000, 1000, dim=p)
  } else {
    prior_coefficients
  }
  serr <- if (is.null(prior_error_sd)) {
    greta::inverse_gamma(1, 1, dim=1)
  } else {
    prior_error_sd
  }

  n.species  <- length(unique(iris$Species))
  species_id <- as.numeric(iris$Species)
  mu         <- zeros(n.species)

  Z     <- matrix(0, nrow(iris), n.species * 2)
  gamma <- zeros(n.species * 2)

  for (s in unique(species_id)) {
    Z[species_id == s, s * 2  - 1] <- 1
    Z[species_id == s, s * 2]  <- iris$Petal.Length[species_id == s]
    ranef_sd <- wishart(5, diag(2))
    gamma[(s * 2 - 1):(s * 2)] <- multivariate_normal(rep(0, 2), ranef_sd)
  }

  wi <- as_data(iris$Sepal.Width)
  Z  <- as_data(Z)

  mu <- int + coef * wi + Z %*% gamma


  sep <- as_data(iris$Sepal.Length)

  distribution(sep) <- normal(mu, sd)
  m <- model(coef, gamma, sd, int)

}