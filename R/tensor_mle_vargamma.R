#' tensor_mle_vargamma
#'
#' Computes the EM estimation algorithm for the vargamma distribution.
#' @return A list containing the estimated mean array, skew array, list of
#'   covariance matrices, and gamma.
#'
#' @noRd
tensor_mle_vargamma <- function(data, max_iter = 1000, tol = 1e-6,
                                quiet = TRUE, restrict = NULL, init = NULL) {

  # get dim of input
  n <- n_draws(data)
  dims <- draw_shape(data)
  o <- length(dims)
  n_star <- prod(dims)

  # Step 1: Initialize vals

  # different params based on model
  gamma <- .tensor_mle_initial_value(init, "gamma", 6)

  has_initial_skew <- !is.null(init) && !is.null(init$skew)

  mu <- .tensor_mle_initial_value(
    init,
    "mu",
    apply(unclass(data), 2:(o + 1), mean)
  )

  if (!is.null(init) && !is.null(init$sigmas)) {
    sigmas <- init$sigmas
  } else {
    normal_start <- tensor_mle_normal(
      data,
      max_iter = max_iter,
      tol = tol,
      quiet = TRUE,
      restrict = NULL
    )
    sigmas <- normal_start$sigmas
  }

  inv_sigma_start <- lapply(sigmas, invert_safe)

  precision_resids_mat <- if (has_initial_skew) NULL else
    matrix(0, nrow = n_star, ncol = n)
  if (!has_initial_skew) for (i in seq_len(n)) {
    centered <- .tensor_single_draw_array(pull_draw(data, i)) - mu

    for (d in seq_along(inv_sigma_start)) {
      centered <- n_prod(centered, d, inv_sigma_start[[d]])
    }

    precision_resids_mat[, i] <- as.numeric(centered)
  }
  precision_resids <- if (has_initial_skew) NULL else
    array(precision_resids_mat, dim = c(dims, n))

  # Use the median shift of precision-weighted residuals to seed skewness.
  skew <- if (has_initial_skew) init$skew else
    -apply(precision_resids, 1:o, median)

  logliks <- rep(0, max_iter)

  for (t in 1:max_iter) {
    # Step 2: Update a, b, c depending on expected values
    skew_compute <- skew
    inv_sigma <- lapply(sigmas, invert_safe)

    for (d in seq_along(sigmas)) {
      skew_compute <- n_prod(skew_compute, d, inv_sigma[[d]])
    }

    rho <- sum(skew * skew_compute)

    delta_vals <- rep(0, n)

    # mu_array <- replicate(n, mu, simplify = "array") |>
    #   aperm(c(o + 1, (1:(o))))

    #centered <- draws - mu_array

    for (i in 1:n) {
      center_draw <- .tensor_single_draw_array(pull_draw(data, i)) - mu

      centered_compute <- center_draw

      for (d in seq_along(sigmas)) {
        centered_compute <- n_prod(centered_compute, d, inv_sigma[[d]])
      }
      delta_vals[i] <- sum(center_draw * centered_compute)
    }

    rho <- rho + 2 * gamma
    param_vals <- gamma - n_star / 2

    bessel_arg <- sqrt(rho * delta_vals)
    k_ratio <- .besselK_asym_ratio(bessel_arg, param_vals + 1, param_vals)

    a <- sqrt(delta_vals / rho) * k_ratio

    b <- sqrt(rho / delta_vals) *
      k_ratio -
      (2 * param_vals) / delta_vals

    eps <- 1e-5

    c <- 1/2 * log(delta_vals / rho) +
      .dlog_besselK_asym_dnu(bessel_arg, param_vals, eps)

    # replace NaN or Inf values
    b[!is.finite(b)] <- 0
    c[!is.finite(c)] <- 0

    # Step 3: Update mu, skew, params
    weight_mean <- mean(a) * b - 1
    weight_skew <- mean(b) - b

    num_mean <- .tensor_weighted_draw_sum(data, weight_mean)
    num_skew <- .tensor_weighted_draw_sum(data, weight_skew)

    den_mean <- sum(mean(a) * b) - n

    new_mu <- num_mean / den_mean

    den_skew <- sum(a * mean(b)) - n
    new_skew <- num_skew / den_skew

    new_sigmas <- sigmas

    # update params based on model
    update_nu <- function(nu, b, c, n) {
      log(nu / 2) + 1 - digamma(nu / 2) - 1 / n * sum(b + c)
    }
    update_gamma <- function(gamma, a, c) {
      log(gamma) + 1 - digamma(gamma) + mean(c) - mean(a)
    }

    new_gamma <- uniroot(
      update_gamma,
      interval = c(1e-3, 1e3),
      a = a,
      c = c)$root

    for (j in 1:o) {
      inv_new_sigma <- lapply(new_sigmas, invert_safe) # compute sigmas

      n_d <- dims[j]
      other_modes <- (1:o)[-j]

      first <- matrix(0, n_d, n_d)
      x_sum  <- matrix(0, n_d, prod(dims[-j]))

      skew_tmp <- new_skew # compute skew once for each dim

      for (d in other_modes) { # multiply by inverse covars
        skew_tmp <- n_prod(skew_tmp, d, inv_new_sigma[[d]])
      }

      flat_skew_tmp <- matricization(skew_tmp, j) # flatten the skews
      flat_skew <- matricization(new_skew, j)

      for (i in 1:n) {
        xm <- .tensor_single_draw_array(pull_draw(data, i)) - new_mu # take centered draw
        flat_xm <- matricization(xm, j)

        x_sum <- x_sum + flat_xm

        xm_tmp <- xm
        for (d in other_modes) { # multiply by inverse covars
          xm_tmp <- n_prod(xm_tmp, d, inv_new_sigma[[d]])
        }

        flat_xm_tmp <- matricization(xm_tmp, j)

        first <- first + b[i] * (flat_xm_tmp %*% t(flat_xm))
      }

      second <- flat_skew_tmp %*% t(x_sum)
      third  <- x_sum %*% t(flat_skew_tmp)
      fourth <- sum(a) * (flat_skew_tmp %*% t(flat_skew))

      sigma_j <- n_d / (n * n_star) * (first - second - third + fourth)

      new_sigmas[[j]] <- sigma_j
    }

    scale_prod <- 1

    if (o > 1L) {
      for (j in seq_len(o - 1L)) { # normalize all but the last covariance
        curr_sigma <- new_sigmas[[j]]
        n_d <- dims[j]

        tr_j <- sum(diag(curr_sigma))

        scale_curr <- tr_j / n_d
        scale_prod <- scale_prod * scale_curr

        curr_sigma <- curr_sigma / scale_curr
        new_sigmas[[j]] <- curr_sigma
      }
    }

    new_sigmas[[o]] <- new_sigmas[[o]] * scale_prod

    # update all parameters
    mu <- new_mu
    skew <- new_skew
    sigmas <- new_sigmas
    gamma <- new_gamma

    # Step 5: Check convergence

    logliks[t] <- sum(dtvargamma(data, mu = mu, skew = skew, sigmas = sigmas,
                                 scale = gamma, log = TRUE))

    if(t >= 3) {
      ll_rel <- abs(logliks[t] - logliks[t - 1]) /
                (abs(logliks[t - 1]) + 1e-8)

      if (ll_rel < tol) {
        if (!quiet) message("Converged at iteration ", t)
        break
      }
    }

    if (t %% 50 == 0 & !quiet) {
      cat(sprintf("Iteration %d: criterion = %.3e\n", t, ll_rel))
    }
  }

  if(t == max_iter) message("Reached max iter ", max_iter)

  k <- 2 * n_star + sum((dims * (dims+1))/2) - (o - 1) + 1

  list(mu = mu, skew = skew, sigmas = sigmas, gamma = gamma,
       Ew = a, Einvw = b, Elogw = c,
       loglik = logliks[t],
       k = k,
       AIC = 2 * k - 2 * logliks[t],
       BIC = k * log(n) - 2 * logliks[t])
}
