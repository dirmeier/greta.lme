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