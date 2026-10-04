# Optional numerical leading-order ARE calculation.
# This uses a large Monte Carlo sample to approximate expected observed information,
# then applies the delta method to beta and a numerical Hessian of classification risk.
# It is not called by default by the simulation runners.

risk_from_beta <- function(beta, theta_true) conditional_classification_error(beta, theta_true)

risk_hessian_beta <- function(beta_true, theta_true) {
  optimHess(beta_true, function(b) risk_from_beta(b, theta_true))
}

expected_info_cc_mc <- function(theta_true, n_mc = 20000L, seed = 1L) {
  set.seed(seed)
  d <- simulate_gaussian2(n_mc, theta_true)
  p <- ncol(d$y)
  par0 <- pack_theta(theta_true)
  fn <- function(par) nll_cc_internal(par, d$y, d$z, p) / n_mc
  H <- optimHess(par0, fn)
  safe_inverse(H)
}

expected_info_recorded_mc <- function(theta_true, xi0, xi1, alpha, n_mc = 20000L,
                                      seed = 1L, entropy_floor = 1e-300) {
  set.seed(seed)
  d <- make_partially_labelled_data(n_mc, theta_true, alpha, xi0, xi1, entropy_floor)
  p <- ncol(d$y)
  par0 <- c(pack_theta(theta_true), xi0 = xi0, eta_xi = log(xi1))
  fn <- function(par) nll_recorded_mixed(par, d, p, alpha, entropy_floor) / n_mc
  H <- optimHess(par0, fn)
  safe_inverse(H)
}

beta_vcov_from_internal <- function(V, par0, p, has_xi = FALSE) {
  if (is.null(V)) return(NULL)
  dt <- theta_length(p)
  f <- function(par) theta_to_beta(unpack_theta(par[seq_len(dt)], p))
  J <- numeric_jacobian(f, par0)
  J %*% V %*% t(J)
}

estimate_leading_are <- function(theta_true, xi0, xi1, alpha, n_mc = 20000L,
                                 seed = 1L, entropy_floor = 1e-300) {
  p <- length(theta_true$mu1)
  btrue <- theta_to_beta(theta_true)
  Hrisk <- risk_hessian_beta(btrue, theta_true)

  Vcc_int <- expected_info_cc_mc(theta_true, n_mc = n_mc, seed = seed)
  Vpc_int <- expected_info_recorded_mc(theta_true, xi0, xi1, alpha,
                                       n_mc = n_mc, seed = seed + 1L,
                                       entropy_floor = entropy_floor)
  if (is.null(Vcc_int) || is.null(Vpc_int)) return(NA_real_)

  pcc <- pack_theta(theta_true)
  ppc <- c(pack_theta(theta_true), xi0 = xi0, eta_xi = log(xi1))
  Vcc_beta <- beta_vcov_from_internal(Vcc_int, pcc, p)
  Vpc_beta <- beta_vcov_from_internal(Vpc_int, ppc, p, has_xi = TRUE)
  num <- sum(diag(Hrisk %*% Vcc_beta))
  den <- sum(diag(Hrisk %*% Vpc_beta))
  num / den
}
