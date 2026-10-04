fit_converged <- function(fit) !is.null(fit) && isTRUE(fit$convergence == 0L)

reported_mixed_vector <- function(fit) {
  c(theta_report(fit$theta), xi0 = unname(fit$xi["xi0"]), xi1 = unname(fit$xi["xi1"]), alpha = fit$alpha)
}

truth_mixed_vector <- function(theta, xi0, xi1, alpha) {
  c(theta_report(theta), xi0 = unname(xi0), xi1 = unname(xi1), alpha = alpha)
}

recorded_se_vectors <- function(fit, p) {
  if (is.null(fit$vcov_internal)) return(list(parameter = NULL, beta = NULL))
  V <- fit$vcov_internal
  dt <- theta_length(p)
  f_report <- function(par) {
    th <- unpack_theta(par[seq_len(dt)], p)
    c(theta_report(th), xi0 = unname(par[dt + 1L]), xi1 = exp(unname(par[dt + 2L])))
  }
  f_beta <- function(par) theta_to_beta(unpack_theta(par[seq_len(dt)], p))
  Jr <- numeric_jacobian(f_report, fit$par)
  Jb <- numeric_jacobian(f_beta, fit$par)
  Vr <- Jr %*% V %*% t(Jr)
  Vb <- Jb %*% V %*% t(Jb)
  ser <- sqrt(pmax(diag(Vr), 0)); names(ser) <- rownames(Jr)
  seb <- sqrt(pmax(diag(Vb), 0)); names(seb) <- rownames(Jb)
  list(parameter = c(ser, alpha = sqrt(fit$alpha * (1 - fit$alpha) / NA_real_)), beta = seb)
}

recorded_se_vectors_with_n <- function(fit, p, n) {
  x <- recorded_se_vectors(fit, p)
  if (!is.null(x$parameter)) x$parameter["alpha"] <- sqrt(fit$alpha * (1 - fit$alpha) / n)
  x
}

classification_metrics <- function(fit, theta_true, test = NULL) {
  if (!fit_converged(fit)) return(c(error = NA_real_, excess = NA_real_, test_error = NA_real_))
  err <- conditional_classification_error(fit$beta, theta_true)
  bayes <- conditional_classification_error(theta_to_beta(theta_true), theta_true)
  te <- NA_real_
  if (!is.null(test)) {
    pred <- predict_class(fit$beta, test$y)
    te <- mean(pred != test$z)
  }
  c(error = err, excess = err - bayes, test_error = te)
}

composite_curve_rmse <- function(fit, theta_true, alpha_true, xi_true, n_eval = 5000L,
                                 entropy_floor = 1e-300) {
  ev <- simulate_gaussian2(n_eval, theta_true)$y
  en_true <- posterior_entropy(ev, theta_true, entropy_floor)
  q_true <- q_from_entropy(en_true, xi_true["xi0"], xi_true["xi1"])
  r_true <- alpha_true + (1 - alpha_true) * q_true

  en_hat <- posterior_entropy(ev, fit$theta, entropy_floor)
  q_hat <- q_from_entropy(en_hat, fit$xi["xi0"], fit$xi["xi1"])
  r_hat <- fit$alpha + (1 - fit$alpha) * q_hat
  sqrt(mean((r_hat - r_true)^2))
}

start_instability <- function(fit) {
  ss <- fit$start_summary
  if (is.null(ss) || nrow(ss) < 2L) return(c(loglik_spread = 0, alpha_spread = NA_real_))
  ll_spread <- diff(range(-ss$objective, finite = TRUE))
  a_spread <- if ("alpha_start_solution" %in% names(ss)) diff(range(ss$alpha_start_solution, finite = TRUE)) else NA_real_
  c(loglik_spread = ll_spread, alpha_spread = a_spread)
}

summarise_scalar <- function(x) {
  c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    median = median(x, na.rm = TRUE),
    q975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE),
    n = sum(is.finite(x)))
}

bind_named_rows <- function(lst) {
  nm <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(x) {
    out <- setNames(rep(NA_real_, length(nm)), nm)
    out[names(x)] <- x
    out
  }))
}
