#' tensor_mle
#'
#' Estimates the mean array and covariance matrices from tensor-valued
#' observations. If tensor variate normal, uses the flop-flop algorithm. Otherwise,
#' uses the expectation-maximization algorithm.
#'
#' @param draws A `tensor` object containing IID observations.
#' @param max_iter A max number of iterations to try to get covariance matrices
#'   that converge.
#' @param tol A tolerance level to define the convergence of matrices.
#' @param model Specify the distribution of draws. Must be one of normal, skewt,
#'    vargamma, invgauss, or genhyper
#' @param init Optional named list of starting parameter values, or a fitted
#'   parameter list returned by `tensor_mle()`. Recognized components depend on
#'   `model`; output-only components such as `loglik`, `AIC`, and `BIC` are
#'   ignored. Missing components use the model's default initialization.
#' @return A list containing the estimated mean array and the list of covariance
#'   matrices.
#'
#' @export
tensor_mle <- function(draws, max_iter = 1000, tol = 1e-4, quiet = TRUE,
                       model, restrict = NULL, init = NULL) {
  if (!inherits(draws, "tensor")) {
    stop(
      "`draws` must be a `tensor` object. Use `tensor()` first.",
      call. = FALSE
    )
  }

  if (!model %in% c("normal", "skewt", "vargamma", "invgauss", "genhyper")) {
    stop("Not a valid model. Must be normal, skewt, vargamma, invgauss, or genhyper")
  } # check there is a real model

  n <- n_draws(draws)

  if (n == 0L) {
    stop("`draws` must contain at least one observation.", call. = FALSE)
  }

  .tensor_mle_validate_init(init, model, draw_shape(draws))

  # call the correct MLE function based on model
  switch(
    model,
    normal = tensor_mle_normal(
      draws, max_iter = max_iter, tol = tol, quiet = quiet,
      restrict = restrict, init = init
    ),
    skewt = tensor_mle_skewt(
      draws, max_iter = max_iter, tol = tol, quiet = quiet,
      restrict = restrict, init = init
    ),
    vargamma = tensor_mle_vargamma(
      draws, max_iter = max_iter, tol = tol, quiet = quiet,
      restrict = restrict, init = init
    ),
    invgauss = tensor_mle_invgauss(
      draws, max_iter = max_iter, tol = tol, quiet = quiet,
      restrict = restrict, init = init
    ),
    genhyper = tensor_mle_genhyper(
      draws, max_iter = max_iter, tol = tol, quiet = quiet,
      restrict = restrict, init = init
    )
  )
}

.tensor_mle_parameter_names <- function(model) {
  switch(
    model,
    normal = c("mu", "sigmas"),
    skewt = c("mu", "sigmas", "skew", "nu"),
    vargamma = c("mu", "sigmas", "skew", "gamma"),
    invgauss = c("mu", "sigmas", "skew", "kappa"),
    genhyper = c("mu", "sigmas", "skew", "lambda", "omega")
  )
}

.tensor_mle_initial_value <- function(init, name, default) {
  if (!is.null(init) && !is.null(init[[name]])) init[[name]] else default
}

.tensor_mle_validate_init <- function(init, model, dims) {
  if (is.null(init)) return(invisible(TRUE))

  if (!is.list(init) || is.null(names(init)) || any(names(init) == "")) {
    stop("`init` must be a named list.", call. = FALSE)
  }

  permitted <- .tensor_mle_parameter_names(model)
  supplied <- intersect(names(init), permitted)
  if (length(supplied) == 0L) {
    stop(
      "`init` does not contain any parameters recognized for model `",
      model, "`.", call. = FALSE
    )
  }

  validate_array <- function(value, name) {
    if (!is.array(value) ||
        !identical(as.integer(dim(value)), as.integer(dims)) ||
        !is.numeric(value) || any(!is.finite(value))) {
      stop(
        "`init$", name,
        "` must be a finite numeric array with dimensions ",
        paste(dims, collapse = " x "), ".", call. = FALSE
      )
    }
  }

  if (!is.null(init$mu)) validate_array(init$mu, "mu")
  if (!is.null(init$skew)) validate_array(init$skew, "skew")

  if (!is.null(init$sigmas)) {
    if (!is.list(init$sigmas) || length(init$sigmas) != length(dims)) {
      stop(
        "`init$sigmas` must contain one covariance matrix per tensor mode.",
        call. = FALSE
      )
    }

    for (mode in seq_along(dims)) {
      sigma <- init$sigmas[[mode]]
      expected_dim <- rep(as.integer(dims[mode]), 2L)
      positive_definite <- FALSE
      if (is.matrix(sigma) && is.numeric(sigma) &&
          identical(dim(sigma), expected_dim) &&
          all(is.finite(sigma)) && isSymmetric(sigma)) {
        positive_definite <- !inherits(
          try(chol(sigma), silent = TRUE), "try-error"
        )
      }
      if (!positive_definite) {
        stop(
          "`init$sigmas[[", mode,
          "]]` must be a finite symmetric positive-definite ",
          dims[mode], " by ", dims[mode], " matrix.", call. = FALSE
        )
      }
    }
  }

  for (name in intersect(c("nu", "gamma", "kappa", "omega"), supplied)) {
    value <- init[[name]]
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop(
        "`init$", name, "` must be one positive finite number.",
        call. = FALSE
      )
    }
  }

  if (!is.null(init$lambda) &&
      (!is.numeric(init$lambda) || length(init$lambda) != 1L ||
       !is.finite(init$lambda))) {
    stop("`init$lambda` must be one finite number.", call. = FALSE)
  }

  invisible(TRUE)
}

.tensor_weighted_draw_sum <- function(data, weights) {
  n <- n_draws(data)
  dims <- draw_shape(data)

  if (!is.numeric(weights) || length(weights) != n || anyNA(weights)) {
    stop("`weights` must contain one non-missing numeric value per draw.",
         call. = FALSE)
  }

  data_mat <- matrix(unclass(data), nrow = n, ncol = prod(dims))
  array(drop(crossprod(weights, data_mat)), dim = dims)
}
