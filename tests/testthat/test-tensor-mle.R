test_that("tensor_mle models retain draw-shaped array parameters", {
  set.seed(20260722)

  dims <- c(2, 2)
  mu <- array(c(-0.5, 0, 0.5, 1), dim = dims)
  skew <- array(c(0.2, -0.1, 0.15, 0.1), dim = dims)
  sigmas <- list(
    matrix(c(1, 0.2, 0.2, 1), nrow = 2),
    matrix(c(1, 0.1, 0.1, 1), nrow = 2)
  )

  draws <- list(
    normal = rtnorm(40, mu = mu, sigmas = sigmas),
    skewt = rtskewt(40, mu = mu, sigmas = sigmas, skew = skew, nu = 8),
    vargamma = rtvargamma(
      40, mu = mu, sigmas = sigmas, skew = skew, scale = 4
    ),
    invgauss = rtinvgauss(
      40, mu = mu, sigmas = sigmas, skew = skew, kappa = 2
    ),
    genhyper = rtgenhyper(
      40, mu = mu, sigmas = sigmas, skew = skew,
      lambda = 1, omega = 4
    )
  )

  fits <- lapply(names(draws), function(model) {
    expect_no_error(
      suppressMessages(
        tensor_mle(draws[[model]], model = model, max_iter = 3)
      )
    )
  })
  names(fits) <- names(draws)

  for (model in names(fits)) {
    fit <- fits[[model]]
    expect_equal(dim(fit$mu), dims, info = model)
    expect_false(inherits(fit$mu, "tensor"), info = model)
    expect_equal(length(fit$sigmas), length(dims), info = model)
    expect_true(is.finite(fit$loglik), info = model)

    if (model != "normal") {
      expect_equal(dim(fit$skew), dims, info = model)
      expect_false(inherits(fit$skew, "tensor"), info = model)
    }
  }
})

test_that("weighted tensor draw sums preserve array shape and values", {
  draws <- tensor(array(seq_len(24), dim = c(3, 2, 4)), obs = 1L)
  weights <- c(0.5, -1, 2)

  observed <- tensortools:::.tensor_weighted_draw_sum(draws, weights)
  expected <- Reduce(
    `+`,
    lapply(seq_len(n_draws(draws)), function(i) {
      weights[i] * tensortools:::.tensor_single_draw_array(
        pull_draw(draws, i)
      )
    })
  )

  expect_equal(observed, expected)
  expect_equal(dim(observed), draw_shape(draws))
  expect_false(inherits(observed, "tensor"))
})

test_that("tensor_mle models handle vector-valued draws", {
  set.seed(20260723)

  dims <- 2L
  mu <- array(c(-0.25, 0.5), dim = dims)
  skew <- array(c(0.15, -0.1), dim = dims)
  sigmas <- list(matrix(c(1, 0.2, 0.2, 1), nrow = 2))

  draws <- list(
    normal = rtnorm(40, mu = mu, sigmas = sigmas),
    skewt = rtskewt(40, mu = mu, sigmas = sigmas, skew = skew, nu = 8),
    vargamma = rtvargamma(
      40, mu = mu, sigmas = sigmas, skew = skew, scale = 4
    ),
    invgauss = rtinvgauss(
      40, mu = mu, sigmas = sigmas, skew = skew, kappa = 2
    ),
    genhyper = rtgenhyper(
      40, mu = mu, sigmas = sigmas, skew = skew,
      lambda = 1, omega = 4
    )
  )

  for (model in names(draws)) {
    fit <- expect_no_error(
      suppressMessages(
        tensor_mle(draws[[model]], model = model, max_iter = 3)
      )
    )

    expect_equal(dim(fit$mu), dims, info = model)
    expect_equal(dim(fit$sigmas[[1L]]), c(dims, dims), info = model)
    expect_true(is.finite(fit$loglik), info = model)

    if (model != "normal") {
      expect_equal(dim(fit$skew), dims, info = model)
    }
  }
})

test_that("tensor_mle models handle scalar draws consistently", {
  set.seed(20260724)
  dims <- 1L
  mu <- array(0.25, dim = dims)
  skew <- array(0.1, dim = dims)
  sigmas <- list(matrix(1.5))

  draws <- list(
    normal = rtnorm(40, mu = mu, sigmas = sigmas),
    skewt = rtskewt(40, mu = mu, sigmas = sigmas, skew = skew, nu = 8),
    vargamma = rtvargamma(
      40, mu = mu, sigmas = sigmas, skew = skew, scale = 4
    ),
    invgauss = rtinvgauss(
      40, mu = mu, sigmas = sigmas, skew = skew, kappa = 2
    ),
    genhyper = rtgenhyper(
      40, mu = mu, sigmas = sigmas, skew = skew,
      lambda = 1, omega = 4
    )
  )

  for (model in names(draws)) {
    fit <- expect_no_error(
      suppressMessages(
        tensor_mle(draws[[model]], model = model, max_iter = 3)
      )
    )

    expect_equal(dim(fit$mu), dims, info = model)
    expect_equal(dim(fit$sigmas[[1L]]), c(dims, dims), info = model)
    expect_true(is.finite(fit$loglik), info = model)

    if (model != "normal") {
      expect_equal(dim(fit$skew), dims, info = model)
    }
  }
})

test_that("tensor_mle accepts fitted objects as warm starts", {
  set.seed(20260823)

  dims <- c(2, 2)
  mu <- array(c(-0.25, 0, 0.25, 0.5), dim = dims)
  skew <- array(c(0.1, -0.05, 0.05, 0.1), dim = dims)
  sigmas <- list(diag(2), diag(2))
  draws <- list(
    normal = rtnorm(30, mu = mu, sigmas = sigmas),
    skewt = rtskewt(
      30, mu = mu, sigmas = sigmas, skew = skew, nu = 8
    ),
    vargamma = rtvargamma(
      30, mu = mu, sigmas = sigmas, skew = skew, scale = 4
    ),
    invgauss = rtinvgauss(
      30, mu = mu, sigmas = sigmas, skew = skew, kappa = 2
    ),
    genhyper = rtgenhyper(
      30, mu = mu, sigmas = sigmas, skew = skew,
      lambda = 1, omega = 4
    )
  )

  for (model in names(draws)) {
    fit <- suppressMessages(
      tensor_mle(draws[[model]], model = model, max_iter = 2)
    )

    warm_fit <- expect_no_error(
      suppressMessages(
        tensor_mle(
          draws[[model]], model = model, max_iter = 1, init = fit
        )
      )
    )

    expect_equal(dim(warm_fit$mu), dims, info = model)
    expect_true(is.finite(warm_fit$loglik), info = model)
  }
})

test_that("tensor_mle validates warm-start parameters", {
  draws <- rtnorm(
    10,
    mu = array(0, dim = c(2, 2)),
    sigmas = list(diag(2), diag(2))
  )

  expect_error(
    tensor_mle(draws, model = "normal", init = 1),
    "named list",
    fixed = TRUE
  )
  expect_error(
    tensor_mle(draws, model = "normal", init = list(loglik = 0)),
    "does not contain any parameters",
    fixed = TRUE
  )
  expect_error(
    tensor_mle(
      draws, model = "normal",
      init = list(sigmas = list(diag(2), matrix(c(1, 2, 2, 1), 2)))
    ),
    "positive-definite",
    fixed = TRUE
  )
  expect_error(
    tensor_mle(
      draws, model = "skewt",
      init = list(nu = -1)
    ),
    "positive finite",
    fixed = TRUE
  )
})
