#' dtvargamma
#'
#' Computes the density function for the tensor variate gamma.
#'
#' @param x A numeric vector/array, a `tensor`, or a `tensor` object.
#' @param mu An array of the mean.
#' @param skew An array of the skew.
#' @param sigmas A list of the sigma covariance matrices.
#' @param scale A parameter that describes the shape and rate of the gamma.
#' @param log Defaults to FALSE. If TRUE, return the log of the density.
#'
#' @return A single density for one tensor, or a numeric vector of densities
#'   for a tensor object.
#'
#' @examples
#' dtvargamma(array(1), mu = 0, skew = array(1), sigmas = list(matrix(1)), scale = 4)
#' @export
dtvargamma <- function(x, mu, skew, sigmas, scale, log = FALSE) {
  x_is_tensor <- inherits(x, "tensor")

  if (is.list(x)) {
    stop(
      "`x` must be a tensor/array or a `tensor` object; lists are not supported.",
      call. = FALSE
    )
  }

  if (x_is_tensor) {
    if (n_draws(x) == 0L) {
      stop("x must contain at least one tensor observation.", call. = FALSE)
    }

    dims <- draw_shape(x)
  } else {
    if (!is.numeric(x)) {
      stop(
        "`x` must be numeric, an array, a `tensor`, or a `tensor` object.",
        call. = FALSE
      )
    }

    dims <- .tensor_dims(x)
  }

  o <- length(dims)
  n_star <- prod(dims)

  .validate_same_dims(mu, dims, "mu")
  .validate_same_dims(skew, dims, "skew")
  sigmas <- .prepare_sigmas(sigmas, dims)

  mu <- array(mu, dim = dims)
  skew <- array(skew, dim = dims)

  all_det <- 0
  skew_tmp <- skew

  ve_skew <- matrix(c(skew), nrow = n_star)
  inv_sigmas <- lapply(sigmas, invert_safe)

  for(k in length(sigmas):1) {
    sigd <- sigmas[[k]]

    skew_tmp <- n_prod(skew_tmp, k, inv_sigmas[[k]])

    all_det <- all_det + (n_star/(2 * nrow(sigd))) * .logdet_safe(sigd)
  }

  skew_tmp <- as.numeric(skew_tmp)

  rho <- sum(skew_tmp * ve_skew)

  eval_one <- function(x_curr) {
    .validate_same_dims(x_curr, dims, "x", reference = "mu")

    if (inherits(x_curr, "tensor")) {
      x_curr <- .tensor_single_draw_array(x_curr)
    } else {
      x_curr <- array(x_curr, dim = dims)
    }

    xm <- x_curr - mu
    xm_tmp <- xm
    ve_xm <- matrix(c(xm), nrow = n_star)

    for(k in length(sigmas):1) {
      xm_tmp <- n_prod(xm_tmp, k, inv_sigmas[[k]])
    }

    xm_tmp <- as.numeric(xm_tmp)

    delta <- sum(xm_tmp * ve_xm)
    xm_skew <- sum(xm_tmp * ve_skew)
    y <- sqrt((rho + 2 * scale) * delta)
    log_bessel <- .log_besselK_asym(y, scale - n_star/2)

    loglik <-
      log(2) + scale * log(scale) - (n_star)/2 * log(2 * pi) - all_det -
      lgamma(scale) + xm_skew +
      (scale - (n_star/2))/2 * log(delta/(rho + 2 * scale)) +
      log_bessel - y

    if(log == FALSE) loglik <- exp(loglik)

    loglik
  }

  if (!x_is_tensor) {
    return(eval_one(x))
  }

  vals <- numeric(n_draws(x))

  for (i in seq_len(n_draws(x))) {
    vals[i] <- eval_one(pull_draw(x, i))
  }

  vals
}
