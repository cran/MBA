# MBA

[![CRAN status](https://www.r-pkg.org/badges/version/MBA)](https://CRAN.R-project.org/package=MBA)
[![R-CMD-check](https://github.com/finleya/MBA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/finleya/MBA/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/finleya/MBA/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/finleya/MBA/actions/workflows/pkgdown.yaml)

<img src="man/figures/logo.png" alt="MBA package hex logo" align="left" width="130" style="margin: 0 1.2rem 0.6rem 0;">

`MBA` interpolates irregularly and regularly spaced data using
multilevel B-spline approximation. The package provides functions for
estimating smooth surfaces on a regular grid and for predicting surface
values at arbitrary point locations.

<br clear="left"/>

## Installation

Install the CRAN release with:

```r
install.packages("MBA")
```

Install the development version from GitHub with:

```r
remotes::install_github("finleya/MBA")
```

## Basic Use

```r
library(MBA)

data(LIDAR)
set.seed(1)

train <- sample(seq_len(nrow(LIDAR)), 500)
xyz <- LIDAR[train, ]

surf <- mba.surf(xyz, no.X = 100, no.Y = 100, extend = TRUE)
image(surf$xyz.est, xaxs = "r", yaxs = "r")

pts <- mba.points(
  xyz = xyz,
  xy.est = LIDAR[-train, c("x", "y")],
  h = 8,
  verbose = FALSE
)

head(pts$xyz.est)
```

## Functionality

The package provides:

- Surface interpolation from bivariate scattered data with `mba.surf()`.
- Point prediction from fitted multilevel B-spline approximations with
  `mba.points()`.
- Optional convex-hull masking for grid surfaces.
- Optional `sp` output for gridded surfaces.
- The `LIDAR` example dataset for interpolation examples.

## Citation

The `MBA` implementation calls portions of the SINTEF Multilevel
B-spline Library written by Oyvind Hjelle and implements methods
developed in:

Lee, S., Wolberg, G., and Shin, S. Y. (1997). Scattered data
interpolation with multilevel B-splines. *IEEE Transactions on
Visualization and Computer Graphics*, 3(3), 229-244.
doi:10.1109/2945.620490.
