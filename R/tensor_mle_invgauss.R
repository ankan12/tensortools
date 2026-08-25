#' tensor_mle_invgauss
#'
#' Computes the EM estimation algorithm for the invgamma distribution.
#' @return A list containing the estimated mean array, skew array, list of
#'   covariance matrices, and kappa,
#'
#' @noRd
tensor_mle_invgauss <- function(data, max_iter, tol,
                                quiet = TRUE, restrict = NULL, init = NULL) {

  n <- n_draws(data)
  dims <- draw_shape(data)
  o <- length(dims)
  n_star <- prod(dims)

  # Step 1: Initialize vals

  # different params based on model
  kappa <- .tensor_mle_initial_value(init, "kappa", 2)
  flat_draws <- aperm(unclass(data), c(seq.int(2L, o + 1L), 1L))

  mean_draws <- apply(flat_draws, 1:o, mean)
  median_draws <- apply(flat_draws, 1:o, median)
  mu <- .tensor_mle_initial_value(init, "mu", median_draws)

  # E[X] = M + 1/kappa * skew
  skew <- if (!is.null(init) && !is.null(init$skew)) init$skew else
    kappa * (mean_draws - median_draws)

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

    rho <- rho + kappa^2
    delta_vals <- delta_vals + 1
    param_vals <- -(1 + n_star) / 2

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
    den_skew <- sum(a * mean(b)) - n

    new_mu <- num_mean / den_mean
    new_skew <- num_skew / den_skew

    new_sigmas <- sigmas

    # update params based on model
    new_kappa <- n / sum(a)

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

    #scale_prod <- 1
    #
    # for(j in 1:(o-1)) { # force trace to be n_d for all sigmas except last
    #   curr_sigma <- new_sigmas[[j]]
    #   n_d <- dims[o]
    #
    #   scale_curr <- sum(diag(curr_sigma)) * n_d
    #   scale_prod <- scale_prod * scale_curr
    #
    #   curr_sigma <- curr_sigma / scale_curr
    #
    #   new_sigmas[[j]] <- curr_sigma
    # }
    #
    # new_sigmas[[o]] <- new_sigmas[[o]] * scale_prod

    # update all parameters
    mu <- new_mu
    skew <- new_skew
    sigmas <- new_sigmas
    kappa <- new_kappa

    # Step 5: Check convergence

    logliks[t] <- sum(dtinvgauss(data, mu = mu, skew = skew,
                                 sigmas = sigmas, kappa = kappa,
                                 log = TRUE))

    if(t >= 3) {
      ll_rel <- abs(logliks[t] - logliks[t - 1]) / (abs(logliks[t - 1]) + 1e-8)

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

  list(mu = mu, skew = skew, sigmas = sigmas, kappa = kappa,
       #Ew = a, Einvw = b, Elogw = c,
       loglik = logliks[t],
       k = k,
       AIC = 2 * k - 2 * logliks[t],
       BIC = k * log(n) - 2 * logliks[t])
}
