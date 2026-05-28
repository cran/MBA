test_that("mba.surf returns an image-style surface", {
  set.seed(1)
  xyz <- cbind(
    x = runif(40),
    y = runif(40),
    z = sin(runif(40) * pi) + cos(runif(40) * pi)
  )

  fit <- mba.surf(xyz, no.X = 12, no.Y = 10, h = 3, extend = TRUE)

  expect_named(
    fit,
    c("xyz.est", "no.X", "no.Y", "n", "m", "h", "extend", "sp", "b.box")
  )
  expect_equal(length(fit$xyz.est$x), 12)
  expect_equal(length(fit$xyz.est$y), 10)
  expect_equal(dim(fit$xyz.est$z), c(10L, 12L))
  expect_true(all(is.finite(fit$xyz.est$z)))
  expect_equal(fit$no.X, 12)
  expect_equal(fit$no.Y, 10)
  expect_true(fit$extend)
})

test_that("mba.points returns estimates at requested locations", {
  set.seed(2)
  xyz <- cbind(
    x = runif(50),
    y = runif(50),
    z = rnorm(50)
  )
  xy.est <- xyz[1:8, c("x", "y")]

  fit <- mba.points(xyz, xy.est, h = 3, verbose = FALSE)

  expect_named(fit, "xyz.est")
  expect_equal(dim(fit$xyz.est), c(8L, 3L))
  expect_equal(fit$xyz.est[, c("x", "y")], xy.est, tolerance = 1e-12)
  expect_true(all(is.finite(fit$xyz.est[, "z"])))
})

test_that("input validation catches malformed calls", {
  expect_error(mba.surf(matrix(1:6, ncol = 2), 5, 5), "3 columns")
  expect_error(mba.surf(matrix(1:9, ncol = 3), 1, 5), "no.X")
  expect_error(mba.points(matrix(1:9, ncol = 3), matrix(1:9, ncol = 3)), "2 columns")
})
