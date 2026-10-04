n_chol_par <- function(p) p * (p + 1L) / 2L

chol_positions <- function(p) {
  out <- vector("list", n_chol_par(p))
  k <- 1L
  for (i in seq_len(p)) {
    for (j in seq_len(i)) {
      out[[k]] <- c(i, j)
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

theta_length <- function(p) 1L + 2L * p + n_chol_par(p)

pack_theta <- function(theta) {
  p <- length(theta$mu1)
  L <- t(chol(theta$Sigma))
  cp <- numeric(n_chol_par(p))
  pos <- chol_positions(p)
  for (k in seq_len(nrow(pos))) {
    i <- pos[k, 1]; j <- pos[k, 2]
    cp[k] <- if (i == j) log(L[i, j]) else L[i, j]
  }
  out <- c(qlogis(theta$pi1), theta$mu1, theta$mu2, cp)
  names(out) <- c("eta_pi",
                  paste0("mu1_", seq_len(p)),
                  paste0("mu2_", seq_len(p)),
                  paste0("chol_", seq_along(cp)))
  out
}

unpack_theta <- function(par, p) {
  stopifnot(length(par) == theta_length(p))
  k <- 1L
  eta_pi <- par[k]; k <- k + 1L
  mu1 <- par[k:(k + p - 1L)]; k <- k + p
  mu2 <- par[k:(k + p - 1L)]; k <- k + p
  cp <- par[k:length(par)]
  L <- matrix(0, p, p)
  pos <- chol_positions(p)
  for (h in seq_len(nrow(pos))) {
    i <- pos[h, 1]; j <- pos[h, 2]
    L[i, j] <- if (i == j) exp(cp[h]) else cp[h]
  }
  list(pi1 = plogis(eta_pi), mu1 = as.numeric(mu1), mu2 = as.numeric(mu2),
       Sigma = L %*% t(L), L = L)
}

theta_bounds <- function(p) {
  d <- theta_length(p)
  lo <- rep(-Inf, d); hi <- rep(Inf, d)
  lo[1] <- -8; hi[1] <- 8
  pos <- chol_positions(p)
  start_cp <- 1L + 2L * p + 1L
  for (h in seq_len(nrow(pos))) {
    idx <- start_cp + h - 1L
    if (pos[h, 1] == pos[h, 2]) {
      lo[idx] <- -5
      hi[idx] <- 5
    } else {
      lo[idx] <- -20
      hi[idx] <- 20
    }
  }
  names(lo) <- names(hi) <- names(pack_theta(list(pi1=.5, mu1=rep(0,p), mu2=rep(0,p), Sigma=diag(p))))
  list(lower = lo, upper = hi)
}

dmvn_log <- function(y, mu, Sigma) {
  y <- as.matrix(y)
  p <- ncol(y)
  R <- try(chol(Sigma), silent = TRUE)
  if (inherits(R, "try-error")) return(rep(-Inf, nrow(y)))
  z <- sweep(y, 2L, mu, "-")
  zz <- backsolve(R, t(z), transpose = TRUE)
  quad <- colSums(zz^2)
  logdet <- 2 * sum(log(diag(R)))
  -0.5 * (p * log(2 * pi) + logdet + quad)
}

model_log_components <- function(y, theta) {
  l1 <- log(theta$pi1) + dmvn_log(y, theta$mu1, theta$Sigma)
  l2 <- log1p(-theta$pi1) + dmvn_log(y, theta$mu2, theta$Sigma)
  lf <- logsumexp2(l1, l2)
  list(log_joint1 = l1, log_joint2 = l2, log_marginal = lf)
}

posterior_entropy <- function(y, theta, entropy_floor = 1e-300) {
  lc <- model_log_components(y, theta)
  tau1 <- plogis(lc$log_joint1 - lc$log_joint2)
  tau2 <- 1 - tau1
  xlogx <- function(x) {
    out <- numeric(length(x))
    nz <- x > 0
    out[nz] <- x[nz] * log(x[nz])
    out
  }
  en <- -(xlogx(tau1) + xlogx(tau2))
  pmax(en, entropy_floor)
}

theta_to_beta <- function(theta) {
  b1 <- as.numeric(solve(theta$Sigma, theta$mu1 - theta$mu2))
  b0 <- -0.5 * sum((theta$mu1 + theta$mu2) * b1) + log(theta$pi1 / (1 - theta$pi1))
  out <- c(beta0 = b0, setNames(b1, paste0("beta", seq_along(b1))))
  out
}

theta_report <- function(theta) {
  p <- length(theta$mu1)
  Sidx <- which(lower.tri(theta$Sigma, diag = TRUE), arr.ind = TRUE)
  svals <- theta$Sigma[Sidx]
  snames <- paste0("Sigma_", Sidx[,1], Sidx[,2])
  c(pi1 = theta$pi1,
    setNames(theta$mu1, paste0("mu1_", seq_len(p))),
    setNames(theta$mu2, paste0("mu2_", seq_len(p))),
    setNames(svals, snames))
}

make_theta <- function(p, Delta, Sigma = diag(p), pi1 = 0.5) {
  mu1 <- rep(0, p); mu2 <- rep(0, p)
  mu1[1] <- Delta / 2
  mu2[1] <- -Delta / 2
  list(pi1 = pi1, mu1 = mu1, mu2 = mu2, Sigma = as.matrix(Sigma))
}

simulate_gaussian2 <- function(n, theta) {
  p <- length(theta$mu1)
  z <- ifelse(runif(n) < theta$pi1, 1L, 2L)
  L <- t(chol(theta$Sigma))
  e <- matrix(rnorm(n * p), n, p) %*% t(L)
  y <- e
  i1 <- z == 1L; i2 <- !i1
  y[i1, ] <- sweep(y[i1, , drop = FALSE], 2L, theta$mu1, "+")
  y[i2, ] <- sweep(y[i2, , drop = FALSE], 2L, theta$mu2, "+")
  list(y = y, z = z)
}

predict_class <- function(beta, y) {
  d <- beta[1] + as.numeric(as.matrix(y) %*% beta[-1])
  ifelse(d > 0, 1L, 2L)
}

conditional_classification_error <- function(beta_hat, theta_true) {
  b0 <- beta_hat[1]; b <- beta_hat[-1]
  v <- as.numeric(t(b) %*% theta_true$Sigma %*% b)
  if (!is.finite(v) || v <= 1e-14) return(NA_real_)
  sdv <- sqrt(v)
  m1 <- b0 + sum(b * theta_true$mu1)
  m2 <- b0 + sum(b * theta_true$mu2)
  theta_true$pi1 * pnorm(0, mean = m1, sd = sdv) +
    (1 - theta_true$pi1) * (1 - pnorm(0, mean = m2, sd = sdv))
}

closed_form_cc <- function(y, z, ridge = 1e-8) {
  y <- as.matrix(y); p <- ncol(y); n <- nrow(y)
  i1 <- z == 1L; i2 <- z == 2L
  if (sum(i1) < 2L || sum(i2) < 2L) stop("Insufficient labelled observations for CC fit.")
  pi1 <- mean(i1)
  mu1 <- colMeans(y[i1, , drop = FALSE])
  mu2 <- colMeans(y[i2, , drop = FALSE])
  r1 <- sweep(y[i1, , drop = FALSE], 2L, mu1, "-")
  r2 <- sweep(y[i2, , drop = FALSE], 2L, mu2, "-")
  Sigma <- (crossprod(r1) + crossprod(r2)) / n
  Sigma <- Sigma + diag(ridge, p)
  theta <- list(pi1 = pi1, mu1 = mu1, mu2 = mu2, Sigma = Sigma)
  list(method = "CC", theta = theta, par = pack_theta(theta), beta = theta_to_beta(theta),
       convergence = 0L, loglik = NA_real_)
}
