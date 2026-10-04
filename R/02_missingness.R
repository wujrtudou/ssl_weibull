q_from_entropy <- function(entropy, xi0, xi1) {
  plogis(xi0 + xi1 * log(pmax(entropy, 1e-300)))
}

missing_eta <- function(y, theta, xi0, xi1, entropy_floor = 1e-300) {
  en <- posterior_entropy(y, theta, entropy_floor = entropy_floor)
  xi0 + xi1 * log(en)
}

binary_entropy_from_score <- function(d, entropy_floor = 1e-300) {
  p1 <- plogis(d)
  p2 <- 1 - p1
  xlogx <- function(x) {
    out <- numeric(length(x))
    nz <- x > 0
    out[nz] <- x[nz] * log(x[nz])
    out
  }
  pmax(-(xlogx(p1) + xlogx(p2)), entropy_floor)
}

calibrate_xi0 <- function(theta, xi1, gamma, tol = 1e-10, entropy_floor = 1e-300) {
  stopifnot(gamma > 0, gamma < 1, xi1 > 0)
  beta <- theta_to_beta(theta)
  b0 <- beta[1]; b <- beta[-1]
  v <- as.numeric(t(b) %*% theta$Sigma %*% b)
  sdv <- sqrt(v)
  md1 <- b0 + sum(b * theta$mu1)
  md2 <- b0 + sum(b * theta$mu2)

  eq <- function(xi0) {
    f <- function(d, md) {
      en <- binary_entropy_from_score(d, entropy_floor)
      plogis(xi0 + xi1 * log(en)) * dnorm(d, mean = md, sd = sdv)
    }
    e1 <- integrate(f, lower = -Inf, upper = Inf, md = md1,
                    rel.tol = 1e-9, subdivisions = 300L)$value
    e2 <- integrate(f, lower = -Inf, upper = Inf, md = md2,
                    rel.tol = 1e-9, subdivisions = 300L)$value
    theta$pi1 * e1 + (1 - theta$pi1) * e2
  }

  lo <- -40; hi <- 40
  flo <- eq(lo) - gamma; fhi <- eq(hi) - gamma
  while (flo > 0) { lo <- lo - 20; flo <- eq(lo) - gamma }
  while (fhi < 0) { hi <- hi + 20; fhi <- eq(hi) - gamma }
  uniroot(function(x) eq(x) - gamma, interval = c(lo, hi), tol = tol)$root
}

generate_missingness <- function(y, theta, alpha, xi0, xi1, entropy_floor = 1e-300) {
  eta <- missing_eta(y, theta, xi0, xi1, entropy_floor)
  q <- plogis(eta)
  n <- nrow(as.matrix(y))
  m1 <- runif(n) < alpha
  m2 <- rep(FALSE, n)
  eligible <- !m1
  m2[eligible] <- runif(sum(eligible)) < q[eligible]
  m <- m1 | m2
  list(m1 = m1, m2 = m2, m = m, q = q, entropy_eta = eta)
}

make_partially_labelled_data <- function(n, theta, alpha, xi0, xi1, entropy_floor = 1e-300) {
  d <- simulate_gaussian2(n, theta)
  mm <- generate_missingness(d$y, theta, alpha, xi0, xi1, entropy_floor)
  z_obs <- d$z
  z_obs[mm$m] <- NA_integer_
  list(y = d$y, z_full = d$z, z = z_obs,
       m1 = mm$m1, m2 = mm$m2, m = mm$m,
       q_true = mm$q, theta_true = theta,
       alpha_true = alpha, xi_true = c(xi0 = xi0, xi1 = xi1))
}
