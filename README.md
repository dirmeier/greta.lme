# greta.lme

[![Project Status](http://www.repostatus.org/badges/latest/concept.svg)](http://www.repostatus.org/#concept)
[![Project Life](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Build travis](https://travis-ci.org/dirmeier/greta.lme.svg?branch=master)](https://travis-ci.org/dirmeier/greta.lme)
[![Build status](https://ci.appveyor.com/api/projects/status/y1q7ynuvc7eo0pt3?svg=true)](https://ci.appveyor.com/project/dirmeier/greta-lme)
[![codecov](https://codecov.io/gh/dirmeier/greta.lme/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/greta.lme)
[![CRAN](http://www.r-pkg.org/badges/version/greta.lme?color=white)](https://cran.r-project.org/package=greta.lme)

a greta extension for hierarchical models 
 
## About

`greta.lme` is an extension package for `greta` to allow easier building of hierarchical models, such as mixed effect models.
For example, a normal hierarchical model using random intercepts for the iris data set in `lme4` formula style:

```R
mod <- greta.glmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris)
sd <- inverse_gamma(1, 1)
 
distribution(iris$Sepal.Length) <- normal(mod$predictor, sd)
m <- model(mod$coef, mod$gamma)
d <- mcmc(m)
```

## Installation

```r
devtools::install_github("greta-dev/greta")
devtools::install_github("dirmeier/greta.lme")
```

The current `greta.lme` release uses the same version numbers as the current `greta` release, such that it's easier to see what works together, e.g.

```r
devtools::install_github("greta-dev/greta@0.3.0")
devtools::install_github("dirmeier/greta.lme@0.3.0")
```


## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@web.de">simon.dirmeier@web.de</a>
