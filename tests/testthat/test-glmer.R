context("glmer")

data(iris)

testthat::test_that("glmer contains every list item", {
  mod <- greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris)

  testthat::expect_equal(dim(mod$predictor)[1], nrow(iris))

  testthat::expect_equal(ncol(mod$X), 2)
  testthat::expect_equal(dim(mod$x.eta)[1], nrow(iris))
  testthat::expect_equal(dim(mod$coef)[1], 2)
  testthat::expect_equal(dim(mod$coef_sd)[1], 1)

  testthat::expect_equal(ncol(mod$Z), 3)
  testthat::expect_equal(nrow(mod$Z), nrow(iris))
  testthat::expect_equal(dim(mod$z.eta)[1], nrow(iris))
  testthat::expect_equal(dim(mod$gamma)[1], length(unique(iris$Species)))
  testthat::expect_equal(dim(mod$gamma_sd)[1], length(unique(iris$Species)))
  testthat::expect_equal(length(mod$Ztlist), 1)
  testthat::expect_equal(length(mod$grp.vars), 1)

  testthat::expect_equal(as.character(mod$call[[1]]), "greta.glmer")
  testthat::expect_equal(as.character(mod$call[[2]]),
                         c("~", "Sepal.Length", "Sepal.Width + (1 | Species)"))
  testthat::expect_equal(as.character(mod$call[[3]]), "iris")
})


testthat::test_that("glmer builds intercept model with custom coefficients", {
  m <- greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                   prior_coefficients = greta::variable(dim=1))
  testthat::expect_true(is.null(m$coef_sd))
})


testthat::test_that("glmer builds intercept model with custom ranef", {
  m <- greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                   prior_random_effects = greta::normal(0, 3, dim=3))
  testthat::expect_true(is.null(m$gamma_sd))
  testthat::expect_true(dim(m$gamma)[1] == 3)
  testthat::expect_true(dim(m$Z)[1] == nrow(iris))
  testthat::expect_true(dim(m$Z)[2] == 3)
})


testthat::test_that("glmer builds slope model with custom coefficients", {
  m <- greta.glmer(Sepal.Length ~ Sepal.Width + (Sepal.Width | Species), iris,
                   prior_random_effects = greta::normal(0, 1, dim=6))
  testthat::expect_true(is.null(m$gamma_sd))
  testthat::expect_true(dim(m$gamma)[1] == 6)
  testthat::expect_true(dim(m$Z)[1] == nrow(iris))
  testthat::expect_true(dim(m$Z)[2] == 6)
})

testthat::test_that("glmer throws when no raneds", {
  testthat::expect_error(
    greta.glmer(Sepal.Length ~ Sepal.Width, iris)
  )
})


testthat::test_that("glmer fails with incorrect intercept prior dim", {
  testthat::expect_error(
    greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                prior_intercept =  greta::uniform(0, 1, dim=2))
  )
})


testthat::test_that("glmer fails with incorrect coefficient prior dim", {
  testthat::expect_error(
    greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                prior_coefficients =  greta::uniform(0, 1, dim=2))
  )
})


testthat::test_that("glmer fails with incorrect ranef prior dim", {
  testthat::expect_error(
    greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                prior_random_effects =  greta::uniform(0, 1, dim=2))
  )
})


testthat::test_that("glmer fails with ranefs are not normal", {
  testthat::expect_error(
    greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                prior_random_effects =  greta::uniform(0, 1, dim=3))
  )
})


testthat::test_that("glmer fails with classes ", {
  testthat::expect_error(
    greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                prior_intercept =  3)
  )
  testthat::expect_error(
    greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                prior_coefficients =  3)
  )
  testthat::expect_error(
    greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                prior_random_effects =  3)
  )
  testthat::expect_error(
    greta.glmer("Sepal.Length ~ Sepal.Width + (1 | Species)", iris)
  )
  testthat::expect_error(
    greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), 2)
  )
})