initial_theta_observed <- function(data, ridge = 1e-4) {
  y <- as.matrix(data$y); p <- ncol(y)
  obs <- !is.na(data$z)
  z <- data$z[obs]; yo <- y[obs, , drop = FALSE]
  if (sum(z == 1L) >= 2L && sum(z == 2L) >= 2L) {
    pi1 <- min(max(mean(z == 1L), 0.05), 0.95)
    mu1 <- colMeans(yo[z == 1L, , drop = FALSE])
    mu2 <- colMeans(yo[z == 2L, , drop = FALSE])
    r1 <- sweep(yo[z == 1L, , drop = FALSE], 2L, mu1, "-")
    r2 <- sweep(yo[z == 2L, , drop = FALSE], 2L, mu2, "-")
    S <- (crossprod(r1) + crossprod(r2)) / max(1, nrow(yo))
  } else {
    km <- kmeans(y, centers = 2, nstart = 5)
    cc <- km$centers
    # orient by first coordinate for a deterministic start
    ord <- order(cc[, 1], decreasing = TRUE)
    mu1 <- cc[ord[1], ]; mu2 <- cc[ord[2], ]
    pi1 <- 0.5
    S <- cov(y)
  }
  S <- as.matrix(S) + diag(ridge, p)
  list(pi1 = pi1, mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), Sigma = S)
}

perturb_start <- function(x, scale = 0.05) x + rnorm(length(x), sd = scale)

fit_pc <- function(data, n_starts = 3L, control = list(), entropy_floor = 1e-300,
                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- ncol(as.matrix(data$y))
  th0 <- initial_theta_observed(data)
  p0 <- pack_theta(th0)
  b <- theta_bounds(p)
  starts <- c(list(p0), replicate(max(0L, n_starts - 1L), perturb_start(p0), simplify = FALSE))
  ms <- run_multistart_nlminb(starts, nll_pc, lower = b$lower, upper = b$upper,
                              control = control, data = data, p = p,
                              entropy_floor = entropy_floor)
  best <- ms$best
  theta <- unpack_theta(best$par, p)
  list(method = "MCAR", par = best$par, theta = theta, beta = theta_to_beta(theta),
       loglik = -best$objective, convergence = best$convergence,
       start_summary = ms$summary, missing_rate = mean(data$m))
}

initial_xi <- function(data, theta, response, eligible = rep(TRUE, nrow(as.matrix(data$y))),
                       entropy_floor = 1e-300) {
  en <- posterior_entropy(data$y, theta, entropy_floor)
  x <- log(en[eligible]); yy <- as.integer(response[eligible])
  if (length(unique(yy)) < 2L || length(yy) < 10L) {
    pr <- min(max(mean(yy), 1e-3), 1 - 1e-3)
    return(c(xi0 = qlogis(pr), eta_xi = 0))
  }
  g <- try(glm(yy ~ x, family = binomial()), silent = TRUE)
  if (inherits(g, "try-error") || any(!is.finite(coef(g)))) {
    pr <- min(max(mean(yy), 1e-3), 1 - 1e-3)
    return(c(xi0 = qlogis(pr), eta_xi = 0))
  }
  a <- unname(coef(g)[1]); slope <- unname(coef(g)[2])
  if (!is.finite(slope) || slope <= 0) slope <- 1
  c(xi0 = a, eta_xi = log(max(slope, 1e-3)))
}

mixed_bounds <- function(p, include_alpha = FALSE, alpha_upper = 1 - 1e-8) {
  tb <- theta_bounds(p)
  lo <- c(tb$lower, xi0 = -40, eta_xi = -7)
  hi <- c(tb$upper, xi0 = 40, eta_xi = 7)
  if (include_alpha) {
    lo <- c(lo, alpha = 0)
    hi <- c(hi, alpha = alpha_upper)
  }
  list(lower = lo, upper = hi)
}

fit_recorded_mixed <- function(data, n_starts = 3L, control = list(), entropy_floor = 1e-300,
                               pc_fit = NULL, compute_hessian = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- ncol(as.matrix(data$y))
  pc_fit <- pc_fit %||% fit_pc(data, n_starts = max(2L, n_starts), control = control,
                               entropy_floor = entropy_floor)
  alpha_hat <- mean(data$m1)
  xi0 <- initial_xi(data, pc_fit$theta, response = data$m2, eligible = !data$m1,
                    entropy_floor = entropy_floor)
  base <- c(pc_fit$par, xi0)
  bd <- mixed_bounds(p, include_alpha = FALSE)
  starts <- c(list(base), replicate(max(0L, n_starts - 1L), perturb_start(base), simplify = FALSE))
  ms <- run_multistart_nlminb(starts, nll_recorded_mixed, lower = bd$lower, upper = bd$upper,
                              control = control, data = data, p = p, alpha = alpha_hat,
                              entropy_floor = entropy_floor)
  best <- ms$best
  dt <- theta_length(p)
  theta <- unpack_theta(best$par[seq_len(dt)], p)
  xi <- c(xi0 = best$par[dt + 1L], xi1 = exp(best$par[dt + 2L]))
  H <- V <- NULL
  if (compute_hessian && best$convergence == 0L) {
    H <- try(optimHess(best$par, nll_recorded_mixed, data = data, p = p,
                       alpha = alpha_hat, entropy_floor = entropy_floor), silent = TRUE)
    if (!inherits(H, "try-error")) V <- safe_inverse(H) else H <- NULL
  }
  list(method = "MIXED_RECORDED", par = best$par, theta = theta, xi = xi, alpha = alpha_hat,
       beta = theta_to_beta(theta), loglik = -best$objective,
       convergence = best$convergence, start_summary = ms$summary,
       hessian = H, vcov_internal = V)
}

fit_mar <- function(data, n_starts = 3L, control = list(), entropy_floor = 1e-300,
                    pc_fit = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- ncol(as.matrix(data$y))
  pc_fit <- pc_fit %||% fit_pc(data, n_starts = max(2L, n_starts), control = control,
                               entropy_floor = entropy_floor)
  xi0 <- initial_xi(data, pc_fit$theta, response = data$m, entropy_floor = entropy_floor)
  base <- c(pc_fit$par, xi0)
  bd <- mixed_bounds(p, include_alpha = FALSE)
  starts <- c(list(base), replicate(max(0L, n_starts - 1L), perturb_start(base), simplify = FALSE))
  ms <- run_multistart_nlminb(starts, nll_mar, lower = bd$lower, upper = bd$upper,
                              control = control, data = data, p = p,
                              entropy_floor = entropy_floor)
  best <- ms$best
  dt <- theta_length(p)
  theta <- unpack_theta(best$par[seq_len(dt)], p)
  xi <- c(xi0 = best$par[dt + 1L], xi1 = exp(best$par[dt + 2L]))
  list(method = "MAR", par = best$par, theta = theta, xi = xi, alpha = 0,
       beta = theta_to_beta(theta), loglik = -best$objective,
       convergence = best$convergence, start_summary = ms$summary)
}

fit_latent_mixed <- function(data, alpha_starts = c(0, 0.05, 0.15, 0.30), control = list(),
                             entropy_floor = 1e-300, alpha_upper = 1 - 1e-8,
                             pc_fit = NULL, compute_hessian = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- ncol(as.matrix(data$y))
  pc_fit <- pc_fit %||% fit_pc(data, n_starts = 3L, control = control,
                               entropy_floor = entropy_floor)
  xi0 <- initial_xi(data, pc_fit$theta, response = data$m, entropy_floor = entropy_floor)
  bd <- mixed_bounds(p, include_alpha = TRUE, alpha_upper = alpha_upper)
  alpha_starts <- unique(pmin(pmax(alpha_starts, 0), alpha_upper))
  starts <- lapply(alpha_starts, function(a) c(pc_fit$par, xi0, alpha = a))
  # add one mild perturbation when only a few starts are supplied
  if (length(starts) < 3L) starts <- c(starts, list(perturb_start(starts[[1]], 0.03)))
  ms <- run_multistart_nlminb(starts, nll_latent_mixed, lower = bd$lower, upper = bd$upper,
                              control = control, data = data, p = p,
                              entropy_floor = entropy_floor)
  best <- ms$best
  dt <- theta_length(p)
  theta <- unpack_theta(best$par[seq_len(dt)], p)
  xi <- c(xi0 = best$par[dt + 1L], xi1 = exp(best$par[dt + 2L]))
  alpha <- best$par[dt + 3L]
  H <- V <- NULL
  if (compute_hessian && best$convergence == 0L && alpha > 1e-8 && alpha < alpha_upper - 1e-8) {
    H <- try(optimHess(best$par, nll_latent_mixed, data = data, p = p,
                       entropy_floor = entropy_floor), silent = TRUE)
    if (!inherits(H, "try-error")) V <- safe_inverse(H) else H <- NULL
  }
  ss <- ms$summary
  ss$alpha_start_solution <- vapply(ms$all, function(z) z$par[dt + 3L], numeric(1))
  list(method = "MIXED_LATENT", par = best$par, theta = theta, xi = xi, alpha = alpha,
       beta = theta_to_beta(theta), loglik = -best$objective,
       convergence = best$convergence, start_summary = ss,
       hessian = H, vcov_internal = V)
}

fit_oracle_latent <- function(data, alpha, xi, n_starts = 3L, control = list(),
                              entropy_floor = 1e-300, pc_fit = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- ncol(as.matrix(data$y))
  pc_fit <- pc_fit %||% fit_pc(data, n_starts = n_starts, control = control,
                               entropy_floor = entropy_floor)
  tb <- theta_bounds(p)
  starts <- c(list(pc_fit$par), replicate(max(0L, n_starts - 1L), perturb_start(pc_fit$par), simplify = FALSE))
  ms <- run_multistart_nlminb(starts, nll_oracle_latent, lower = tb$lower, upper = tb$upper,
                              control = control, data = data, p = p, alpha = alpha,
                              xi0 = xi["xi0"], xi1 = xi["xi1"], entropy_floor = entropy_floor)
  best <- ms$best
  theta <- unpack_theta(best$par, p)
  list(method = "ORACLE_LATENT", par = best$par, theta = theta, xi = xi, alpha = alpha,
       beta = theta_to_beta(theta), loglik = -best$objective,
       convergence = best$convergence, start_summary = ms$summary)
}

fit_cc <- function(data) closed_form_cc(data$y, data$z_full)

fit_mcar <- function(data, ...) fit_pc(data, ...)
