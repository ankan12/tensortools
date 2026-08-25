#' tensor_mle_genhyper
#'
#' Computes the EM estimation algorithm for the genhyper distribution.
#' @return A list containing the estimated mean array, skew array, list of
#'   covariance matrices, lambda, and omega.
#'
#' @noRd
tensor_mle_genhyper <- function(data, max_iter = 1000, tol = 1e-6,
                                quiet = TRUE, restrict = NULL, init = NULL) {

  # get dim of input
  n <- n_draws(data)
  dims <- draw_shape(data)
  o <- length(dims)
  n_star <- prod(dims)

  # Step 1: Initialize vals

  # different params based on model
  lambda <- .tensor_mle_initial_value(init, "lambda", 10)
  omega <- .tensor_mle_initial_value(init, "omega", 10)

  flat_draws <- aperm(unclass(data), c(seq.int(2L, o + 1L), 1L))

  mean_draws <- apply(flat_draws, 1:o, mean)
  median_draws <- apply(flat_draws, 1:o, median)
  mu <- .tensor_mle_initial_value(init, "mu", median_draws)

  # E[X] = M + K_{lambda+1}(omega)/K_{lambda}(omega) * skew

  R_lambda <- .besselK_asym_ratio(omega, lambda + 1, lambda)

  skew <- if (!is.null(init) && !is.null(init$skew)) init$skew else
    1/R_lambda * (mean_draws - median_draws)

  # res_normal <- tensor_mle_normal(
  #   data,
  #   max_iter = max_iter,
  #   tol = tol,
  #   quiet = TRUE,
  #   restrict = NULL
  # )
  #
  # mu <- res_normal$mu
  # sigmas <- res_normal$sigmas
  #
  # inv_sigma_start <- lapply(sigmas, invert_safe)
  #
  # precision_resids <- lapply(data, function(x) {
  #   centered <- x - mu
  #
  #   for (d in seq_along(inv_sigma_start)) {
  #     centered <- n_prod(centered, d, inv_sigma_start[[d]])
  #   }
  #
  #   centered
  # })
  #
  # precision_resids <- simplify2array(precision_resids)

  # Use the median shift of precision-weighted residuals to seed skewness.
  #skew <- -apply(precision_resids, 1:o, median)

  logliks <- rep(0, max_iter)
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

  for (t in 1:max_iter) {
    # Step 2: Update a, b, c depending on expected values
    skew_compute <- skew
    inv_sigma <- lapply(sigmas, invert_safe)

    for (d in seq_along(sigmas)) {
      skew_compute <- n_prod(skew_compute, d, inv_sigma[[d]])
    }

    rho <- sum(skew * skew_compute)

    delta_vals <- rep(0, n)

    for (i in 1:n) {
      center_draw <- .tensor_single_draw_array(pull_draw(data, i)) - mu

      centered_compute <- center_draw

      for (d in seq_along(sigmas)) {
        centered_compute <- n_prod(centered_compute, d, inv_sigma[[d]])
      }

      delta_vals[i] <- sum(center_draw * centered_compute)
    }

    rho <- rho + omega
    delta_vals <- delta_vals + omega
    param_vals <- lambda - n_star / 2

    bessel_arg <- sqrt(rho * delta_vals)
    k_ratio <- .besselK_asym_ratio(bessel_arg, param_vals + 1, param_vals)

    a <- sqrt(delta_vals / rho) * k_ratio

    b <- sqrt(rho / delta_vals) *
      k_ratio -
      (2 * param_vals) / delta_vals

    eps <- 1e-5

    c <- 1/2 * log(delta_vals / rho) +
      .dlog_besselK_asym_dnu(bessel_arg, param_vals, eps)

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

    new_lambda <- mean(c) * lambda *
                  1 / .dlog_besselK_asym_dnu(omega, lambda, eps)

    R_lambda <- .besselK_asym_ratio(omega, new_lambda + 1, new_lambda)

    R_neg_lambda <- .besselK_asym_ratio(omega, -new_lambda + 1, -new_lambda)

    first_deriv <- 1 / 2 * (R_lambda + R_neg_lambda - (mean(a) + mean(b)))

    second_deriv <- 1 /2 * (R_lambda^2 - (1 + 2 * new_lambda) / omega * R_lambda -
                            1 + R_neg_lambda^2 -
                            (1 - 2 * new_lambda) / omega * R_neg_lambda - 1)

    new_omega <- omega - first_deriv / second_deriv

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

      if (j < o) { # force trace to be n_d for all sigmas except last
        sigma_j <- sigma_j / sum(diag(sigma_j)) * n_d
      }

      new_sigmas[[j]] <- sigma_j
    }

    # update all parameters
    mu <- new_mu
    skew <- new_skew
    sigmas <- new_sigmas

    lambda <- new_lambda
    omega <- new_omega

    # Step 5: Check convergence

    logliks[t] <- sum(dtgenhyper(data, mu = mu, skew = skew,
                                 sigmas = sigmas, lambda = lambda,
                                 omega = omega, log = TRUE))

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

  k <- 2 * n_star + sum((dims * (dims+1))/2) - (o - 1) + 2

  list(mu = mu, skew = skew, sigmas = sigmas,
       lambda = lambda, omega = omega,
       #Ew = a, Einvw = b, Elogw = c,
       loglik = logliks[t],
       k = k,
       AIC = 2 * k - 2 * logliks[t],
       BIC = k * log(n) - 2 * logliks[t])
}
