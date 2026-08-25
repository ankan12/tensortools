#' tensor_mle_normal
#'
#' Estimates the mean array and covariance matrices from an array of tensor variate normal draws.
#'
#' @param data An array containing the draws, where the first mode represents each draw.
#' @param max_iter A max number of iterations to try to get covariance matrices that converge.
#' @param tol A tolerance level to define the convergence of matrices.
#' @param quiet Whether to suppress convergence messages.
#' @param restrict Optional covariance scale restrictions.
#' @param init Optional named list containing a starting value for `sigmas`.
#'   The normal mean is estimated directly as the sample mean.
#'
#' @return A list containing the estimated mean array and the list of covariance matrices.
#'
#' @examples
#' all_dims <- c(2, 10, 10, 2)
#' fourth_mu = array(rnorm(400), dim = all_dims)
#' fourth_sigma = lapply(seq_along(all_dims), function(k) diag(all_dims[k]))
#' fourth_tensor <- rtnorm(n = 1e3, mu = fourth_mu, sigmas = fourth_sigma)
#' est_fourth <- mle_est(fourth_tensor)
#' mean(fourth_mu - est_fourth$mu)^2
#' fourth_scaled <- lapply(fourth_sigma, function(S) S * (nrow(S) / sum(diag(S))))
#' (mse_each <- mapply(function(est, true) mean((est - true)^2), est_fourth$sigmas, fourth_scaled))
#' @export
tensor_mle_normal <- function(data, max_iter = 1000, tol = 1e-6,
                              quiet = TRUE, restrict = NULL, init = NULL) {
  #get dim of input
  n <- n_draws(data)
  dims <- draw_shape(data)
  o <- length(dims)
  n_star <- prod(dims)

  if (length(restrict) >= o) {
    stop("Invalid restriction: You must restrict at most number of dims - 1 scale parameters.")
  }

  if(o == 1 && dims == 1) { # univariate
    vector_data <- as.vector(unclass(data))
    mu_value <- mean(vector_data)
    mu <- array(mu_value, dim = dims)
    sigma <- matrix(mean((vector_data - mu_value)^2), nrow = 1L, ncol = 1L)
    total_loglik <- sum(dnorm(
      vector_data, mean = mu_value, sd = sqrt(sigma[1L, 1L]),
      log = TRUE
    ))
    k <- 2
    return(list(
      mu = mu,
      sigmas = list(sigma),
      loglik = total_loglik,
      k = k,
      AIC = 2 * k - 2 * total_loglik,
      BIC = k * log(n) - 2 * total_loglik
    ))
  }

  if(o == 1) { # multivariate
    array_data <- matrix(unclass(data), nrow = n, ncol = dims)
    mu <- array(colMeans(array_data), dim = dims)
    sigma <- cov(array_data) * (n - 1)/n
    total_loglik <- sum(dtnorm(data, mu = mu, sigmas = list(sigma),
                               log = TRUE))
    k <- dims + dims * (dims + 1)/2
    return(list(
      mu = mu,
      sigmas = list(sigma),
      loglik = total_loglik,
      k = k,
      AIC = 2 * k - 2 * total_loglik,
      BIC = k * log(n) - 2 * total_loglik
    ))
  }

  mu <- apply(unclass(data), 2:(o+1), mean)

  logliks <- rep(0, max_iter)

  if (!is.null(init) && !is.null(init$sigmas)) {
    sigmas <- init$sigmas
  } else {
    sigmas <- lapply(dims, diag)

    for(k in 1:o) {
      tot_sum <- 0
      for(i in 1:n) {
        curr_unfold <- matricization(
          .tensor_single_draw_array(pull_draw(data, i)) - mu,
          k
        )

        tot_sum <- tot_sum + tcrossprod(curr_unfold)
      }

      tot_sum <- tot_sum * dims[k]/(n * n_star)

      tot_sum <- tot_sum/(sum(diag(tot_sum))) * dims[k]

      sigmas[[k]] <- tot_sum
    }
  }

  for (t in 1:max_iter) {
    for (j in 1:o) {
      #inverses of all sigma
      inv_sigma <- lapply(sigmas, invert_safe)

      # dim for j-th mode
      n_d <- dims[j]
      other_modes <- (1:o)[-j]

      sigma_j <- matrix(0, nrow = n_d, ncol = n_d)

      # accumlate mode-j covar estimate
      for (draw in 1:n) {
        xm <- .tensor_single_draw_array(pull_draw(data, draw)) - mu
        xm_tmp <- xm
        flat_xm <- matricization(xm, j)

        for(d in other_modes) {
          xm_tmp <- n_prod(xm_tmp, d, inv_sigma[[d]])
        }

        sigma_j <- sigma_j + matricization(xm_tmp, j) %*% t(flat_xm)
      }

      sigma_j <- sigma_j * n_d / (n * n_star)

      if(j %in% restrict) {
        sigma_j <- diag(n_d)
      }
      else if(j < o) { # force trace to be n_d
        sigma_j <- sigma_j / sum(diag(sigma_j)) * n_d
      }

      sigmas[[j]] <- sigma_j

      # if (k %in% restrict) { # restrict to be product of identity
      #   s_k <- s_k / (n * d_negk)
      #   curr_sigma <- sum(diag(s_k)) / d_k
      #   sigmas[[k]] <- diag(d_k) * curr_sigma
      #   next
      # }
      # update sigma_k
      #sigmas[[k]] <- s_k / (n * d_negk)
    }

    # Step 5: Check convergence
    logliks[t] <- sum(dtnorm(data, mu = mu, sigmas = sigmas, log = TRUE))

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

  k <- n_star + sum((dims * (dims+1))/2) - (o - 1)

  for(r in restrict) {
    k <- k - (dims[r] * (dims[r] + 1))/2 + 1
  }

  list(mu = mu, sigmas = sigmas, loglik = logliks[t],
       k = k,
       AIC = 2 * k - 2 * logliks[t],
       BIC = k * log(n) - 2 * logliks[t])
}
