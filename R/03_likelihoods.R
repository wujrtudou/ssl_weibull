known_label_logjoint <- function(data, lc) {
  z <- data$z
  out <- rep(NA_real_, length(z))
  i1 <- !is.na(z) & z == 1L
  i2 <- !is.na(z) & z == 2L
  out[i1] <- lc$log_joint1[i1]
  out[i2] <- lc$log_joint2[i2]
  out
}

nll_pc <- function(theta_par, data, p, entropy_floor = 1e-300) {
  theta <- unpack_theta(theta_par, p)
  lc <- model_log_components(data$y, theta)
  lj <- known_label_logjoint(data, lc)
  obs <- !data$m
  val <- sum(lj[obs]) + sum(lc$log_marginal[data$m])
  if (!is.finite(val)) return(1e100)
  -val
}

nll_recorded_mixed <- function(par, data, p, alpha, entropy_floor = 1e-300) {
  dt <- theta_length(p)
  theta <- unpack_theta(par[seq_len(dt)], p)
  xi0 <- par[dt + 1L]
  xi1 <- exp(par[dt + 2L])

  lc <- model_log_components(data$y, theta)
  lj <- known_label_logjoint(data, lc)
  en <- posterior_entropy(data$y, theta, entropy_floor)
  eta <- xi0 + xi1 * log(en)
  lq <- log_sigmoid(eta)
  l1q <- log1m_sigmoid(eta)

  i0 <- !data$m1 & !data$m2
  i1 <- data$m1
  i2 <- data$m2

  val <- 0
  if (any(i0)) val <- val + sum(lj[i0] + log1p(-alpha) + l1q[i0])
  if (any(i1)) {
    if (alpha <= 0) return(1e100)
    val <- val + sum(lc$log_marginal[i1] + log(alpha))
  }
  if (any(i2)) val <- val + sum(lc$log_marginal[i2] + log1p(-alpha) + lq[i2])
  if (!is.finite(val)) return(1e100)
  -val
}

nll_mar <- function(par, data, p, entropy_floor = 1e-300) {
  dt <- theta_length(p)
  theta <- unpack_theta(par[seq_len(dt)], p)
  xi0 <- par[dt + 1L]
  xi1 <- exp(par[dt + 2L])
  lc <- model_log_components(data$y, theta)
  lj <- known_label_logjoint(data, lc)
  en <- posterior_entropy(data$y, theta, entropy_floor)
  eta <- xi0 + xi1 * log(en)
  lq <- log_sigmoid(eta)
  l1q <- log1m_sigmoid(eta)
  obs <- !data$m
  mis <- data$m
  val <- sum(lj[obs] + l1q[obs]) + sum(lc$log_marginal[mis] + lq[mis])
  if (!is.finite(val)) return(1e100)
  -val
}

nll_latent_mixed <- function(par, data, p, entropy_floor = 1e-300) {
  dt <- theta_length(p)
  theta <- unpack_theta(par[seq_len(dt)], p)
  xi0 <- par[dt + 1L]
  xi1 <- exp(par[dt + 2L])
  alpha <- par[dt + 3L]
  if (alpha < 0 || alpha >= 1) return(1e100)

  lc <- model_log_components(data$y, theta)
  lj <- known_label_logjoint(data, lc)
  en <- posterior_entropy(data$y, theta, entropy_floor)
  eta <- xi0 + xi1 * log(en)
  lq <- log_sigmoid(eta)
  l1q <- log1m_sigmoid(eta)

  obs <- !data$m
  mis <- data$m
  log1ma <- log1p(-alpha)
  lmiss <- if (alpha == 0) {
    lq
  } else {
    logsumexp2(rep(log(alpha), length(lq)), log1ma + lq)
  }

  val <- sum(lj[obs] + log1ma + l1q[obs]) +
    sum(lc$log_marginal[mis] + lmiss[mis])
  if (!is.finite(val)) return(1e100)
  -val
}

nll_oracle_latent <- function(theta_par, data, p, alpha, xi0, xi1,
                              entropy_floor = 1e-300) {
  theta <- unpack_theta(theta_par, p)
  lc <- model_log_components(data$y, theta)
  lj <- known_label_logjoint(data, lc)
  en <- posterior_entropy(data$y, theta, entropy_floor)
  eta <- xi0 + xi1 * log(en)
  lq <- log_sigmoid(eta)
  l1q <- log1m_sigmoid(eta)
  obs <- !data$m
  mis <- data$m
  log1ma <- log1p(-alpha)
  lmiss <- if (alpha == 0) lq else logsumexp2(rep(log(alpha), length(lq)), log1ma + lq)
  val <- sum(lj[obs] + log1ma + l1q[obs]) +
    sum(lc$log_marginal[mis] + lmiss[mis])
  if (!is.finite(val)) return(1e100)
  -val
}

nll_cc_internal <- function(theta_par, y, z, p) {
  theta <- unpack_theta(theta_par, p)
  lc <- model_log_components(y, theta)
  val <- sum(ifelse(z == 1L, lc$log_joint1, lc$log_joint2))
  if (!is.finite(val)) return(1e100)
  -val
}
