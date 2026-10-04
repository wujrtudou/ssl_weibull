make_ar1_sigma <- function(p, rho = 0.3) {
  outer(seq_len(p), seq_len(p), function(i, j) rho^abs(i - j))
}

parameter_recovery_configs <- function() {
  base <- list(n = 1000L, Delta = 2, gamma = 0.30, xi1 = 3, alpha = 0.10, p = 3L)
  cfgs <- list(base)
  cfgs <- c(cfgs, list(modifyList(base, list(n = 500L))))
  for (v in c(1, 3)) cfgs <- c(cfgs, list(modifyList(base, list(Delta = v))))
  for (v in c(0.15, 0.45)) cfgs <- c(cfgs, list(modifyList(base, list(gamma = v))))
  for (v in c(1, 5)) cfgs <- c(cfgs, list(modifyList(base, list(xi1 = v))))
  for (v in c(0, 0.30)) cfgs <- c(cfgs, list(modifyList(base, list(alpha = v))))
  for (i in seq_along(cfgs)) {
    cfgs[[i]]$id_num <- i
    cfgs[[i]]$id <- make_id("recovery", paste0("n", cfgs[[i]]$n), paste0("D", cfgs[[i]]$Delta),
                            paste0("g", cfgs[[i]]$gamma), paste0("x", cfgs[[i]]$xi1),
                            paste0("a", cfgs[[i]]$alpha))
  }
  cfgs
}

efficiency_configs <- function() {
  g <- expand.grid(n = 1000L, p = c(1L, 5L), Delta = c(1,2,3),
                   alpha = c(0,0.10,0.20,0.30,0.40), KEEP.OUT.ATTRS = FALSE)
  g2 <- expand.grid(n = 500L, p = c(1L,5L), Delta = 2,
                    alpha = c(0,0.10,0.20,0.30,0.40), KEEP.OUT.ATTRS = FALSE)
  g <- rbind(g, g2)
  g$gamma <- 0.30; g$xi1 <- 3
  lapply(seq_len(nrow(g)), function(i) {
    x <- as.list(g[i, ]); x$id_num <- i
    x$id <- make_id("eff", paste0("n",x$n), paste0("p",x$p), paste0("D",x$Delta), paste0("a",x$alpha))
    x
  })
}

fixed_missing_configs <- function() {
  g <- expand.grid(Gamma = c(0.30,0.50), p = c(1L,5L), alpha = c(0,0.10,0.20),
                   KEEP.OUT.ATTRS = FALSE)
  g$n <- 1000L; g$Delta <- 2; g$xi1 <- 3
  lapply(seq_len(nrow(g)), function(i) {
    x <- as.list(g[i, ])
    x$gamma <- (x$Gamma - x$alpha) / (1 - x$alpha)
    x$id_num <- i
    x$id <- make_id("fixed", paste0("G",x$Gamma), paste0("p",x$p), paste0("a",x$alpha))
    x
  })
}

latent_configs <- function() {
  g <- expand.grid(n = c(500L,1000L), Delta = c(1,2), alpha = c(0,0.10,0.20,0.30),
                   KEEP.OUT.ATTRS = FALSE)
  g$p <- 5L; g$gamma <- 0.30; g$xi1 <- 3
  lapply(seq_len(nrow(g)), function(i) {
    x <- as.list(g[i, ]); x$id_num <- i
    x$id <- make_id("latent", paste0("n",x$n), paste0("D",x$Delta), paste0("a",x$alpha))
    x
  })
}

prepare_config_truth <- function(cfg, recovery = FALSE) {
  Sigma <- if (recovery) make_ar1_sigma(cfg$p, 0.3) else diag(cfg$p)
  theta <- make_theta(cfg$p, cfg$Delta, Sigma = Sigma, pi1 = 0.5)
  xi0 <- calibrate_xi0(theta, cfg$xi1, cfg$gamma)
  cfg$theta <- theta
  cfg$xi0 <- xi0
  cfg
}

bootstrap_ratio <- function(num, den, R = 2000L, seed = 1L) {
  ok <- is.finite(num) & is.finite(den)
  num <- num[ok]; den <- den[ok]
  if (length(num) < 5L) return(c(se = NA, lo = NA, hi = NA))
  set.seed(seed)
  vals <- replicate(R, {
    ii <- sample.int(length(num), replace = TRUE)
    mean(num[ii]) / mean(den[ii])
  })
  c(se = sd(vals), lo = quantile(vals, 0.025, names = FALSE), hi = quantile(vals, 0.975, names = FALSE))
}

as_config_df <- function(cfg) {
  keep <- setdiff(names(cfg), c("theta"))
  as.data.frame(cfg[keep], stringsAsFactors = FALSE)
}

summarise_recovery_config <- function(res, cfg) {
  ok <- vapply(res, function(x) is.null(x$error) && isTRUE(x$mixed_convergence == 0L), logical(1))
  rr <- res[ok]
  if (!length(rr)) return(data.frame(config = cfg$id, quantity = NA, parameter = NA, n_valid = 0))

  make_rows <- function(truth_name, est_name, se_name, quantity) {
    truth <- rr[[1]][[truth_name]]
    pars <- names(truth)
    do.call(rbind, lapply(pars, function(nm) {
      est <- vapply(rr, function(x) x[[est_name]][nm], numeric(1))
      se <- vapply(rr, function(x) {
        sv <- x[[se_name]]
        if (is.null(sv) || !nm %in% names(sv)) NA_real_ else unname(sv[nm])
      }, numeric(1))
      tr <- unname(truth[nm])
      coverage <- if (all(!is.finite(se))) NA_real_ else mean(abs(est-tr) <= 1.96*se, na.rm=TRUE)
      data.frame(config=cfg$id, quantity=quantity, parameter=nm, truth=tr,
                 bias=mean(est-tr, na.rm=TRUE), rmse=rmse(est,tr),
                 emp_sd=sd(est, na.rm=TRUE), avg_se=mean(se, na.rm=TRUE),
                 coverage95=coverage, n_valid=sum(is.finite(est)), stringsAsFactors=FALSE)
    }))
  }

  r1 <- make_rows("truth_param", "est_param", "se_param", "model_parameter")
  r2 <- make_rows("beta_truth", "beta_est", "beta_se", "discriminant")
  out <- rbind(r1, r2)
  out$coverage95[out$parameter == "alpha" & out$truth == 0] <- NA_real_
  out
}

summarise_efficiency_config <- function(res, cfg, bootstrap_R = 2000L) {
  ok <- vapply(res, function(x) is.null(x$error), logical(1))
  rr <- res[ok]
  cc <- vapply(rr, function(x) x$cc_excess, numeric(1))
  pc <- vapply(rr, function(x) x$mixed_excess, numeric(1))
  br <- bootstrap_ratio(cc, pc, R = bootstrap_R, seed = 5000L + cfg$id_num)
  data.frame(config = cfg$id, n = cfg$n, p = cfg$p, Delta = cfg$Delta,
             alpha = cfg$alpha, gamma = cfg$gamma, xi1 = cfg$xi1,
             mean_cc_excess = mean(cc, na.rm = TRUE),
             mean_mixed_excess = mean(pc, na.rm = TRUE),
             simulated_RE = mean(cc, na.rm = TRUE) / mean(pc, na.rm = TRUE),
             RE_boot_se = br["se"], RE_lo = br["lo"], RE_hi = br["hi"],
             n_valid = sum(is.finite(cc) & is.finite(pc)), stringsAsFactors = FALSE)
}

summarise_latent_config <- function(res, cfg) {
  ok <- vapply(res, function(x) is.null(x$error), logical(1))
  rr <- res[ok]
  get <- function(nm) vapply(rr, function(x) x[[nm]] %||% NA_real_, numeric(1))
  out <- data.frame(
    config = cfg$id, n = cfg$n, p = cfg$p, Delta = cfg$Delta, alpha = cfg$alpha,
    alpha_hat_mean = mean(get("latent_alpha"), na.rm = TRUE),
    alpha_hat_bias = mean(get("latent_alpha") - cfg$alpha, na.rm = TRUE),
    alpha_hat_rmse = rmse(get("latent_alpha"), cfg$alpha),
    xi0_bias = mean(get("latent_xi0") - cfg$xi0, na.rm = TRUE),
    xi0_rmse = rmse(get("latent_xi0"), cfg$xi0),
    xi1_bias = mean(get("latent_xi1") - cfg$xi1, na.rm = TRUE),
    xi1_rmse = rmse(get("latent_xi1"), cfg$xi1),
    alpha_boundary_frequency = mean(get("latent_alpha") <= 1e-10, na.rm = TRUE),
    Dr_mean = mean(get("Dr"), na.rm = TRUE),
    recorded_error = mean(get("recorded_error"), na.rm = TRUE),
    oracle_error = mean(get("oracle_error"), na.rm = TRUE),
    latent_error = mean(get("latent_error"), na.rm = TRUE),
    mar_error = mean(get("mar_error"), na.rm = TRUE),
    mcar_error = mean(get("mcar_error"), na.rm = TRUE),
    cc_error = mean(get("cc_error"), na.rm = TRUE),
    recorded_excess = mean(get("recorded_excess"), na.rm = TRUE),
    oracle_excess = mean(get("oracle_excess"), na.rm = TRUE),
    latent_excess = mean(get("latent_excess"), na.rm = TRUE),
    mar_excess = mean(get("mar_excess"), na.rm = TRUE),
    mcar_excess = mean(get("mcar_excess"), na.rm = TRUE),
    cc_excess = mean(get("cc_excess"), na.rm = TRUE),
    latent_beta_rmse = mean(get("latent_beta_rmse"), na.rm = TRUE),
    loglik_start_spread = mean(get("loglik_start_spread"), na.rm = TRUE),
    alpha_start_spread = mean(get("alpha_start_spread"), na.rm = TRUE),
    n_valid = length(rr), stringsAsFactors = FALSE)
  out
}


summarise_latent_beta_config <- function(res, cfg) {
  ok <- vapply(res, function(x) is.null(x$error) && !is.null(x$latent_beta), logical(1))
  rr <- res[ok]
  if (!length(rr)) return(data.frame(config=cfg$id, parameter=NA, n_valid=0))
  truth <- rr[[1]]$beta_true
  do.call(rbind, lapply(names(truth), function(nm) {
    est <- vapply(rr, function(x) x$latent_beta[nm], numeric(1))
    tr <- unname(truth[nm])
    data.frame(config=cfg$id, parameter=nm, truth=tr,
               bias=mean(est-tr,na.rm=TRUE), rmse=rmse(est,tr),
               emp_sd=sd(est,na.rm=TRUE), n_valid=sum(is.finite(est)),
               stringsAsFactors=FALSE)
  }))
}

summarise_recovery_classification_config <- function(res, cfg) {
  ok <- vapply(res, function(x) is.null(x$error), logical(1))
  rr <- res[ok]
  get <- function(nm) vapply(rr, function(x) x[[nm]] %||% NA_real_, numeric(1))
  data.frame(
    config=cfg$id, n=cfg$n, p=cfg$p, Delta=cfg$Delta, gamma=cfg$gamma,
    xi1=cfg$xi1, alpha=cfg$alpha,
    mixed_test_error=mean(get("mixed_test_error"),na.rm=TRUE),
    pc_test_error=mean(get("pc_test_error"),na.rm=TRUE),
    cc_test_error=mean(get("cc_test_error"),na.rm=TRUE),
    mixed_conditional_error=mean(get("mixed_error"),na.rm=TRUE),
    pc_conditional_error=mean(get("pc_error"),na.rm=TRUE),
    cc_conditional_error=mean(get("cc_error"),na.rm=TRUE),
    convergence_frequency=mean(get("mixed_convergence")==0,na.rm=TRUE),
    valid_information_frequency=mean(vapply(rr,function(x)isTRUE(x$valid_information),logical(1))),
    n_valid=length(rr), stringsAsFactors=FALSE)
}
